module GenomeFragments

## This is a module for handling high throughput sequencing reads

using GenomicFeatures, DataStructures, BioAlignments, Mmap, UnalignedVectors, XAM, IterTools

include("buildfm.jl")
include("processing.jl")

export FragMatrix, FragMatrixPair, FragMatrixSingle, get_frags, chrom_lt, load_frag_matrix, totalfrags, streambam, qualityfilt, filterpairfirst, filterpairfirstpp, 
       barcode, ishighquality, isunique,
       inc_meta_frag!, inc_meta_cut!, inc_meta_mid!,
       inc_heat_frag!, inc_heat_cut!, inc_heat_mid!,
       inc_heat_atac_cut!, inc_heat_atac_frag!, inc_meta_atac_frag!,
       inc_meta_mid_width, inc_heat_mid_width,
       inc_heat_atac_cut_width, inc_meta_atac_cut_width

abstract type FragMatrix end

struct FragMatrixPair{T} <: FragMatrix

    name::String

    numregions::Int
    totalfrags::Int
    numfields::Int
    fields::Vector{String}
    totalecfrags::Vector{Int}

    fragstats::Dict{Symbol, Int}

    index::Dict{String, UnitRange{Int}}
    chromindex::Dict{Int, String}

    FM::T
end


struct FragMatrixSingle{T} <: FragMatrix

    name::String
    readlength::Int
    numregions::Int
    totalfrags::Int
    numfields::Int
    fields::Vector{String}
    totalecfrags::Vector{Int}

    fragstats::Dict{Symbol, Int}

    index::Dict{String, UnitRange{Int}}
    chromindex::Dict{Int, String}

    FM::T
end

function Base.show(io::IO, fm::FragMatrixPair)
    if get(io, :compact, false)
        ### compact
        print(io, "FM Paired Indexed: ", fm.name)
        print(io, "    ", fm.numregions)
        print(io, "    ", fm.totalfrags)
    else
        println(io, "Index FM      :   ", fm.name)
        #println(io, "paired        :   ", fm.paired)
        println(io, "numregions    :   ", fm.numregions)
        println(io, "totalfrags    :   ", fm.totalfrags)
        println(io, "numfields     :   ", fm.numfields)
        println(io, "fields        :   ", fm.fields)
        println(io, "totalecfrags  :   ", fm.totalecfrags)
        println(io, "index         :   ", fm.index)
        println(io, "chromindex    :   ", fm.chromindex)
        println(io, "FM            :   ", summary(fm.FM))
    end
end

function Base.show(io::IO, fm::FragMatrixSingle)
    if get(io, :compact, false)
        ### compact
        print(io, "FM Paired Indexed: ", fm.name)
        print(io, "    ", fm.numregions)
        print(io, "    ", fm.totalfrags)
    else
        println(io, "Index FM      :   ", fm.name)
        println(io, "Readlength    :   ", fm.readlength)
        #println(io, "paired        :   ", fm.paired)
        println(io, "numregions    :   ", fm.numregions)
        println(io, "totalfrags    :   ", fm.totalfrags)
        println(io, "numfields     :   ", fm.numfields)
        println(io, "fields        :   ", fm.fields)
        println(io, "totalecfrags  :   ", fm.totalecfrags)
        println(io, "index         :   ", fm.index)
        println(io, "chromindex    :   ", fm.chromindex)
        println(io, "FM            :   ", summary(fm.FM))
    end
end

##### Find indexed frags
##### Return empty view if chrom is missing
function get_frags(chrom, location::UnitRange{Int}, fm::FragMatrix)

    !haskey(fm.index, chrom) && return view(fm.FM, :, 1:0)

    ind = fm.index[chrom]
    starts = view(fm.FM, 1, ind)
    stops  = view(fm.FM, 2, ind)

    st = min(searchsortedfirst(starts, location.start), searchsortedfirst(stops, location.start))
    (st > size(starts, 1)) && return view(fm.FM, :, 1:0)

    en = min(max(searchsortedlast(starts, location.stop), searchsortedlast(stops, location.stop)), length(stops))

    view(fm.FM, :, ind.start .+ (st:en) .- 1)
end

function get_frags(chrom, fm::FragMatrix)
    !haskey(fm.index, chrom) && return view(fm.FM, :, 1:0)
    ind = fm.index[chrom]
    view(fm.FM, :, ind)
end

#### sorts chroms in numerical order then X, Y, M
function chrom_lt(chrA, chrB)
    chrA[4] == 'M' && return false
    chrB[4] == 'M' && return true

    if isnumber(chrA[4])
        !isnumber(chrB[4]) && return true
        return isless(parse(Int, chrA[4:end]), parse(Int, chrB[4:end]))
    else
        isnumber(chrB[4]) && return false
        return isless(chrA[4:end], chrB[4:end])
    end
end


function load_frag_matrix(file, paired=true)

     io = open(file)

    paired, numregions, totalfrags, numfields, fields, totalecfrags, datatype, index, chromindex, fragstats = read_header_footer(io, paired)

    #  @show paired, numregions, totalfrags, numfields, fields, totalecfrags, datatype, index, chromindex, fragstats

    a = Mmap.mmap(io, Vector{UInt8}, sizeof(datatype)*numfields*totalfrags, position(io))
    # Create an array of the desired eltype and size:
    ua = UnalignedVector{datatype}(a)

    F = reshape(ua, (numfields, totalfrags))

    ft = typeof(F)
    if paired
        FragMatrixPair{ft}(file, numregions, totalfrags, numfields, fields, totalecfrags, fragstats, index, chromindex, F)
    else
        ### estimate readlength
        readlength = estreadlength(F)
        FragMatrixSingle{ft}(file, readlength, numregions, totalfrags, numfields, fields, totalecfrags, fragstats, index, chromindex, F)
    end
end



### for single reads
function estreadlength(F, nr=min(500, size(F, 2)))
    readlengths = counter(Int)
    for i = 1:nr
        readlength = F[2, i] - F[1, i] + 1
        push!(readlengths, readlength)
    end
    rls = collect(keys(readlengths))
    c = [readlengths[k] for k in rls]
    si = sortperm(c, rev=true)

    if rls[si[1]] != maximum(rls)
        println("Warning: maxfraglength ≠ most common frag length: ", rls[si[1]], " ≠ ", maximum(rls))
    end
    return rls[si[1]]
end

function read_header_footer(io, paired=true)

    numregions = -1
    numfields = -1
    totalfrags = -1
    fmfields = String[]
    totalecfrags = Int[]
    dataencoding = Int32
    index = Dict{String, UnitRange{Int}}()
    chromindex = Dict{Int,  String}()


    fragstats = Dict{Symbol, Int}()

    bottom_header = false
    first_line = readline(io)

    if chomp(first_line) == "#bottom_header"
        headerpos = read(io, Int64)
        bottom_header = true
        seek(io, headerpos)
    else
        seekstart(io)
    end

    ### assumes index is final field
    while true
         fields = split(readline(io), '\t')
        (fields[1] == "#paired")     && (paired     = parse(Bool, chomp(fields[2]))::Bool)
        (fields[1] == "#numregions")   && (numregions   = parse(Int,  fields[2])::Int)
        (fields[1] == "#totalfrags") && (totalfrags = parse(Int,  fields[2])::Int)
        (fields[1] == "#numfields")  && (numfields  = parse(Int,  fields[2])::Int)
        (fields[1] == "#fields")     && (fmfields     = split(chomp(fields[2]), [',', ' ', '[', ']'], keepempty=false))
        (fields[1] == "#totalecfrags") && (totalecfrags = parse.([Int], split(chomp(fields[2]), [',', ' ', '[', ']'], keepempty=false))::Vector{Int})
        (fields[1] == "#dataencoding")   && (dataencoding = get_datatype(chomp(fields[2])))


        if (fields[1] == "#numregions") || (fields[1] == "#totalfrags") || (fields[1] == "#totaldupfrags")
            fragstats[Symbol(fields[1][2:end])] = parse(Int,  fields[2])::Int
        end

        if fields[1] == "#index"
            index_str = split(replace(chomp(fields[2]), r"\(|\"" => ""), r"\)[,]*", keepempty=false)
            index_vec = parse_chrom_range.(index_str)

            index = Dict(index_vec)
            chromindex = Dict(i => c[1] for (i, c) in enumerate(index_vec))
            break
        end
    end

    if isempty(totalecfrags)
        totalecfrags = [totalfrags]
    end

    if bottom_header
        seekstart(io)
        readline(io)
        read(io, Int64)
    end

    paired, numregions, totalfrags, numfields, fmfields, totalecfrags, dataencoding, index, chromindex, fragstats


end

function parse_chrom_range(s)
    chrom, rind = split(s, ",")
    sstart, sstop = split(rind, ":")
    chrom, parse(Int, sstart):parse(Int, sstop)
end


function get_datatype(s)
    (s == "Int32")  && return Int32
    (s == "UInt32") && return UInt32
    (s == "Int64")  && return Int64
    (s == "UInt64") && return UInt64

    error("Unsupported data type: $s")
end


####### Encode Decode of strand_chrom
set_strand_chrom_enc(strand, chromindex, ::Type{UInt32}, topbit=0x80000000) = (strand == STRAND_NEG) ?  UInt32(chromindex) | topbit : UInt32(chromindex)
get_strand_chrom_enc(enc::UInt32, topbit=0x80000000) = (enc & topbit) > 0 ? (STRAND_NEG, enc $ topbit) : (STRAND_POS, enc)
set_strand_chrom_enc(strand, chromindex, ::Type{UInt64}, topbit=0x8000000000000000) = (strand == STRAND_NEG) ?  UInt64(chromindex) | topbit : UInt64(chromindex)
get_strand_chrom_enc(enc::UInt64, topbit=0x8000000000000000) = (enc & topbit) > 0 ? (STRAND_NEG, enc $ topbit) : (STRAND_POS, enc)
set_strand_chrom_enc(strand, chromindex, ::Type{Int32}) = (strand == STRAND_NEG) ? -Int32(chromindex) : Int32(chromindex)
get_strand_chrom_enc(enc::Int32) = (enc < 0) ? (STRAND_NEG, -enc) : (STRAND_POS, enc)
set_strand_chrom_enc(strand, chromindex, ::Type{Int64}) = (strand == STRAND_NEG) ? -Int64(chromindex) : Int64(chromindex)
get_strand_chrom_enc(enc::Int64) = (enc < 0) ? (STRAND_NEG, -enc) : (STRAND_POS, enc)


isneg(enc::UInt32, topbit=0x80000000)         = (enc & topbit) > 0 ? true : false
isneg(enc::UInt64, topbit=0x8000000000000000) = (enc & topbit) > 0 ? true : false
isneg(enc::Int32)                             =  enc           < 0 ? true : false
isneg(enc::Int64)                             =  enc           < 0 ? true : false




# chromsome, strand = GenomeFrags2.get_strand_chrom_enc(v[3, i])
# STRAND_POS
# STRAND_NEG


#### Accessor
totalfrags(FM) = FM.totalfrags
totalfrags(FM, s::Symbol) = totalfrags(FM, string(s))
function totalfrags(FM, s::String)
    ind = findfirst(f -> f == s, FM.fields)
    (ind == 0) && error("Field: $s not found")
    FM.totalecfrags[s - 3]
end


end # module

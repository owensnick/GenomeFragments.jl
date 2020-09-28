
using XAM, GenomicFeatures
using IterTools

qualityfilt(τ) = r -> BAM.mappingquality(r) >= τ
@inline isunique(record) = BAM.mappingquality(record) == 255
#@inline isreverse(record)   = (BAM.flag(record) & SAM.FLAG_REVERSE) != 0
@inline ispair(record)      = (BAM.flag(record) & SAM.FLAG_PAIRED)  != 0
@inline firstinpair(record) = (BAM.flag(record) & SAM.FLAG_READ1)   != 0

@inline posstrand(record)  = (BAM.position(record), BAM.ispositivestrand(record))

@inline getstrand(record) = BAM.ispositivestrand(record) ? STRAND_POS : STRAND_NEG
@inline filterpairfirst(record) = BAM.ismapped(record) && (!ispair(record) ? true : (((BAM.position(record) == BAM.nextposition(record)) && BAM.ispositivestrand(record)) || (BAM.templength(record) ≥ 0)))

#### Function to get fragment coords based on flag (see classify_record.jl)
#### Works for single and paired end
#### For single end reads coords are calculated
#### For paired end the pos_pair_filter will guarantee that we only ever consider the leftmsot positive strand read
#### For paired Set the strand to the strand of first in pair

@inline function get_frag_coords(record)

    pos        = BAM.position(record)
    fraglength = abs(BAM.templength(record))
    readlength = BAM.seqlength(record)
    strand     = getstrand(record)

    #### Get start and stop coords and set strand to the strand of the first in pair
    if ispair(record)
        if strand == STRAND_POS
            start_pos = pos
            stop_pos  = pos + fraglength - 1
            !firstinpair(record) && (strand = STRAND_NEG)
        else
            start_pos = pos + readlength - fraglength
            stop_pos  = pos + readlength - 1
            !firstinpair(record) && (strand = STRAND_POS)
        end
    else
        start_pos = pos
        stop_pos  = pos + readlength - 1
    end

    start_pos, stop_pos, strand
end

function fragcoords(record)
    if ispair(record)
        if BAM.ispositivestrand(record)
            pos = BAM.position(record)
            return BAM.refname(record), pos, pos + BAM.templength(record) - 1, ifelse(firstinpair(record), STRAND_POS, STRAND_NEG)
        else
            epos = BAM.rightposition(record)
            return BAM.refname(record), epos - BAM.templength(record), epos, ifelse(firstinpair(record), STRAND_NEG, STRAND_POS)
        end
    else
        return BAM.refname(record), BAM.position(record), BAM.rightposition(record), getstrand(record)
    end
end


function test_bam(bf="test\\test.bam")

    haspaired = false
    bamreader = open(BAM.Reader, bf)
    chroms = [v["SN"] for v in findall(BAM.header(bamreader), "SQ")]
    equivs=[(isunique, seqname)]
    ng = 0
    for group in IterTools.groupby(posstrand, Iterators.filter(filterpairfirst, bamreader))
        sort!(group, by=BAM.templength)
        ng += 1
        println(length(group))
        # println(first(group))

        for record in group
            leftpos    = BAM.position(record)
            rightpos   = BAM.rightposition(record)
            fraglength = abs(BAM.templength(record))
            readlength = BAM.seqlength(record)
            strand     = getstrand(record)
            println((leftpos, rightpos, leftpos + readlength - 1, fraglength, readlength, strand, BAM.mappingquality(record)))
        end

        println("================")


        for tlg in IterTools.groupby(BAM.templength, group)
            haspaired |= ispair(first(tlg))
            chrom, start, stop, strand = fragcoords(first(tlg))
            fragcounts = [length(unique(ef, filter(ff, tlg))) for (ff, ef) in equivs]

            println((chrom, start, stop, strand, fragcounts))
            #all(iszero, fragcounts) && continue


            # ### Get read start, stop and strand
            # ### Count group sizes after splitting into equivalence classes, e.g. by read name or barcode
            # fragstart, fragstop, fragstrand = get_frag_coords(template_group[1])
            # frag_counts = [length(unique(ef, filter(ff, template_group))) for (ff, ef) in equiv_funs_filter]
            # !any(!iszero, frag_counts) && continue
            # chrom = BAM.refname(template_group[1])
            # chromind = chrom_to_ind[chrom]
            #
            # # display((fragstart, fragstop, set_strand_chrom_enc(fragstrand, chromind, UInt32), (fc for fc in frag_counts)...))
            # write(io, T(fragstart))         ### start
            # write(io, T(fragstop))          ### stop
            # write(io, set_strand_chrom_enc(fragstrand, chromind, T))          ### chrom strand enc
            #
            # #### first column total written frags
            # total_frags[chromind, 1] += 1
            # #### remaining cols, equiv_fun frags
            # for (i, fc) in enumerate(frag_counts)
            #     write(io, T(fc))
            #     total_frags[chromind, i+1] += fc
            # end
        end
    end
    @show ng
    println("=====================================")
    close(bamreader)



    chroms
end
test_bam()




function streambam(bamfile, outfile, equivs=[(isunique, seqname)], labels=["uniseq"], T=Int32)
    println("[CFB]\tCollecting frags and counting equivalences streaming....")
    println("[CFB]\tEquiv funs   : ", equivs)
    println("[CFB]\tEquiv labels : ", labels)
    println("[CFB]\tWriting to   : ", outfile)

    io = open(outfile, "w")
    ##### write header
    println(io, "#bottom_header")
    pos_buff = 0
    write(io, pos_buff)


    starttime = time()
    reader = open(BAM.Reader, bamfile)

    chroms, totalfrags, paired = stream_reads_equiv_class_filter(reader, io, equivs, T)
    @show sum(totalfrags)
    display(totalfrags)


    index_table = (chrom=chroms, frags=totalfrags[:, 1], frags_ec=[totalfrags[i, 2:end] for i = 1:size(totalfrags, 1)])
    display(index_table)

    display(sum(index_table.frags_ec))
    display(sum(totalfrags[:, 2:end], dims=1))
    # display(index_table)
    # println("[CFB]\tWriting header:")

    # pos_buff = position(io)

    # println(io, "#paired\t",   paired)
    # println(io, "#numregions\t",   length(chroms))
    # println(io, "#totalfrags\t", sum(index_table[:frags]))
    # println(io, "#totalecfrags\t", sum(index_table[:frags_ec]))
    # println(io, "#numfields\t", 3 + length(equiv_funs))


    println("[CFB]\tWriting header:")

    pos_buff = position(io)

    println(io, "#paired\t",   paired)
    println(io, "#numregions\t",   length(chroms))
    println(io, "#totalfrags\t", sum(totalfrags[:, 1]))
    println(io, "#totalecfrags\t", vec(sum(totalfrags[:, 2:end], dims=1)))
    println(io, "#numfields\t", 3 + length(equivs))
    println(io, "#fields\t[", join(["start" ; "stop" ; "strand_chrom_enc"; labels], ", "), "]")
    println(io, "#dataencoding\t", T)

    ### write index
    frags    = totalfrags[:, 1]
    cumfrags = cumsum(frags)
    fragind  = map((s, e) -> s:e, [0 ; cumfrags[1:end-1]] .+ 1, cumfrags)
    println(io, "#index\t", join(zip(chroms, fragind), ","))

    ##### rewrite header
    seekstart(io)
    println(io, "#bottom_header")
    write(io, pos_buff)

    close(io)
    println("[CFB]\tcomplete in ", time() - starttime, " seconds")

end





function stream_reads_equiv_class_filter(bamreader, io, equivs=[(isunique, seqname)], T=Int32)


    haspaired = false
    chroms = [v["SN"] for v in findall(BAM.header(bamreader), "SQ")]
    chromindex = Dict(chroms[i] => i for i = 1:length(chroms))
    totalfrags = zeros(Int, length(chroms), length(equivs) + 1)

    index = 0
    for group in IterTools.groupby(posstrand, Iterators.filter(filterpairfirst, bamreader))    
        sort!(group, by=BAM.templength)
        for template_group in IterTools.groupby(BAM.templength, group)
            
            ispair(template_group[1]) && (haspaired = true)

            ### Get read start, stop and strand
            ### Count group sizes after splitting into equivalence classes, e.g. by read name or barcode
            chrom, fragstart, fragstop, fragstrand = fragcoords(template_group[1])
            fragcounts = [length(unique(ef, filter(ff, template_group))) for (ff, ef) in equivs]
            !any(!iszero, fragcounts) && continue
            
            chromind = chromindex[chrom]

            # display((fragstart, fragstop, set_strand_chrom_enc(fragstrand, chromind, UInt32), (fc for fc in frag_counts)...))
            write(io, T(fragstart))         ### start
            write(io, T(fragstop))          ### stop
            write(io, set_strand_chrom_enc(fragstrand, chromind, T))          ### chrom strand enc

            #### first column total written frags
            totalfrags[chromind, 1] += 1
            #### remaining cols, equiv_fun frags
            for (i, fc) in enumerate(fragcounts)
                write(io, T(fc))
                totalfrags[chromind, i+1] += fc
            end
        end
    end
    chroms, totalfrags, haspaired
end
streambam("c:\\home\\julia\\test\\test.bam", "c:\\home\\julia\\test\\test.bin", [(qualityfilt(0), seqname), (qualityfilt(20), seqname)], ["q0", "q20"])


streambam("c:\\home\\projects\\islets\\endoc\\k9\\bw2_k9me2.sort.bam", "c:\\home\\projects\\islets\\endoc\\k9\\bw2_k9me2.q30.sort.bin", [(qualityfilt(30), seqname)], ["q30"])
streambam("c:\\home\\projects\\islets\\endoc\\k9\\bw2_k9me3.sort.bam", "c:\\home\\projects\\islets\\endoc\\k9\\bw2_k9me3.q30.sort.bin", [(qualityfilt(30), seqname)], ["q30"])
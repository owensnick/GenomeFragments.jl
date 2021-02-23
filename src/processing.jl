
############## param types

abstract type AccParams end

struct AccParamsSingle{V, W} <: AccParams
    fraglength::Int
    incfun::V
    dataentry::Int
    countfun::W
end

struct AccParamsPair{V, W} <: AccParams
    minfragsize::Int
    maxfragsize::Int
    incfun::V
    dataentry::Int
    countfun::W
end




### helper functions


@inline inc_meta_frag!(V, mpl, mpr, w) =  V[intersect(mpl:mpr, 1:size(V, 1))] .+= w
@inline inc_meta_cut!(V, mpl, mpr, w)  =  (inc_meta_index!(V, mpl, w) ; inc_meta_index!(V, mpr, w))
@inline inc_meta_mid!(V, mpl, mpr, w)  =  inc_meta_index!(V, div(mpl + mpr, 2), w)
@inline inc_meta_index!(V, i::Int, w)  =  (1 ≤ i ≤ size(V, 1)) && (V[i] += w)


@inline inc_heat_frag!(V, k, mpl, mpr, w) =  V[intersect(mpl:mpr, 1:size(V, 1)), k] .+= w
@inline inc_heat_cut!(V, k, mpl, mpr, w)  =  (inc_heat_index!(V, k, mpl, w) ; inc_heat_index!(V, k, mpr, w))
@inline inc_heat_left_cut!(V, k, mpl, mpr, w)  =  inc_heat_index!(V, k, mpl, w)
@inline inc_heat_mid!(V, k, mpl, mpr, w)  =  inc_heat_index!(V, k, div(mpl + mpr, 2), w)
@inline @inbounds inc_heat_index!(V, k, i::Int, w)  =  (1 ≤ i ≤ size(V, 1)) && (V[i, k] += w)



inc_meta_mid_width(bw) = (V, mpl, mpr, w) -> begin
    mid = div(mpl + mpr, 2)
    inc_meta_frag!(V, mid - bw, mid + bw, w)
end

inc_heat_mid_width(bw) = (V, k, mpl, mpr, w) -> begin
    mid = div(mpl + mpr, 2)
    bc = intersect((mid - bw):(mid + bw), 1:size(V, 1))
    V[bc, k] .+= w
end

inc_heat_atac_cut_width(bw) = (V, k, mpl, mpr, w, ao=4) -> begin
    inc_heat_frag!(V, k, mpl + ao - bw, mpl + ao + bw, w)
    inc_heat_frag!(V, k, mpr - ao - bw, mpr - ao + bw, w)
end

inc_meta_atac_cut_width(bw) = (V, mpl, mpr, w, ao=4) -> begin
    inc_meta_frag!(V, mpl + ao - bw, mpl + ao + bw, w)
    inc_meta_frag!(V, mpr - ao - bw, mpr - ao + bw, w)
end

@inline inc_heat_atac_cut!(V, k, mpl, mpr, w, ao=4)   =  (inc_heat_index!(V, k, mpl + ao, w) ; inc_heat_index!(V, k, mpr - ao, w))
@inline inc_heat_atac_frag!(V, k, mpl, mpr, w, ao=4)  =  inc_heat_frag!(V, k, mpl+ao, mpr-ao, w);
@inline inc_meta_atac_frag!(V, mpl, mpr, w, ao=4)     =  inc_meta_frag!(V, mpl+ao, mpr-ao, w);


###############################
function readsintersecting(chroms, locations, FM, dataentry=4, fn=identity)
    total = zeros(Int, length(chroms))
    for (k, (c, l)) in enumerate(zip(chroms, locations))
        
        v = get_frags(c, l, FM)
        for i = 1:size(v, 2)
            total[k] += fn(v[dataentry, i])
        end
    end
    total
end

function fragregion(chrom, location, strand, FM::FragMatrixSingle{T}, inc_fun=inc_meta_frag!, fraglength=120, dataentry=4, fn=identity; pos=true, neg=true) where {T}
    fw = length(location)
    P = zeros(fw)
    V = get_frags(chrom, location, FM)
    
    for i = 1:size(V, 2)
        rstrand, _ = GenomeFragments.get_strand_chrom_enc(V[3, i])
        !pos && (rstrand == STRAND_POS) && continue
        !neg && (rstrand == STRAND_NEG) && continue
        if rstrand == STRAND_POS
            s = V[1, i] - location.start + 1
            e = s + fraglength - 1
        else
            e = V[2, i]  - location.start + 1
            s = e - fraglength + 1
        end
        if strand == "-"
            s, e = fw - e + 1, fw - s + 1
        end
        w = fn(V[dataentry, i])
        inc_fun(P, s, e,  w) 
    end
    P
end


@inline function fragregion(chrom, location, FM::FragMatrixPair{T}, inc_fun=inc_meta_frag!, minfrag=0, maxfrag=1000, dataentry=4, fn=identity) where {T}
    fw = length(location)
    P = zeros(fw)
    V = get_frags(chrom, location, FM)
    
    for col in eachcol(V)
        fp = col[2] - col[1] + 1
        s = col[1] - location.start + 1
        e = col[2] - location.start + 1
        w = fn(col[dataentry])
        (minfrag ≤ fp ≤ maxfrag) && inc_fun(P, s, e, w)
    end

    P
end


#################### fragheatmap



function fragheatmap(chroms, locations, strands, FM::Vector{FragMatrixSingle{T}}, inc_fun=inc_heat_mid!, fraglength=120, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw, length(chroms))
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, length(chroms))
        fragheatmap!(TH, p, chroms, locations, strands, F, inc_fun, fraglength, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
	    elseif norm_type == :frag_region
	        H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        H .*= norm_scale/sum(F.totalecfrags[data_entry - 3] for F in FM)
    else
        H ./= length(FM)
    end
    H
end


function libsize(FM, dataentry, fn=identity, minfrag=0, maxfrag=1000000)

    σ = 0.0
    for i = 1:size(FM.FM, 2)
        fp = FM.FM[2, i] - FM.FM[1, i] + 1
        σ += ifelse(minfrag ≤ fp ≤ maxfrag, fn(FM.FM[dataentry, i]), 0.0)
    end
    σ
end


function libsize(FM, dataentry, fn=identity)
    σ = 0
    @inbounds for i = 1:size(FM.FM, 2)
        σ += fn(FM.FM[dataentry, i])
    end
    σ
end



function fragheatmap(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw, length(chroms))
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, length(chroms))
        fragheatmap!(TH, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	elseif norm_type == :frag_region
	    H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        H .*= 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM)
    else
        H ./= length(FM)
    end
    H
end





function fragheatmap!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])
    
    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
        for col in eachcol(V)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            ### flip missing?

            w = fn(V[data_entry, i])
            inc_fun(H, k, mpl, mpr,  w) 
        end
    end
    H
end

function fragheatmap!(H, p, chroms, locations, strands, FM::FragMatrixSingle{T}, inc_fun=inc_heat_mid!, fraglength=120, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])
    
    for (k, (c, l, str)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, (l.start - 2*fraglength):(l.stop + 2*fraglength), FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            strand, _ = GenomeFragments.get_strand_chrom_enc(V[3, i])
            
            if strand == STRAND_POS
                s = V[1, i] - l.start + 1
                e = s + fraglength - 1
            else
                e = V[2, i]  - l.start + 1
                s = e - fraglength + 1
            end
            if str == "-"
                s, e = fw - e + 1, fw - s + 1
            end
            w = fn(V[data_entry, i])
            inc_fun(H, k, s, e,  w) 
        end
    end
    H
end


function fragheatmap!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, P::AccParamsPair{V, W}, show_progress=true) where {T, V, W}

    fw = length(locations[1])
    
    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
        for col in eachcol(V)
            fp = V[2, i] - V[1, i] + 1
            ((fp < P.minfragsize) || (fp > P.maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            ### flip missing?

            w = P.countfun(V[data_entry, i])
            P.incfun(H, k, mpl, mpr,  w) 
        end
    end
    H
end


function fragheatmap!(H, p, chroms, locations, strands, FM::FragMatrixSingle{T}, P::AccParamsSingle{V, W}; show_progress=true) where {T, V, W}

    fw = length(locations[1])
    
    for (k, (c, l, str)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, (l.start - 2*P.fraglength):(l.stop + 2*P.fraglength), FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            strand, _ = GenomeFragments.get_strand_chrom_enc(V[3, i])
            
            if strand == STRAND_POS
                s = V[1, i] - l.start + 1
                e = s + fraglength - 1
            else
                e = V[2, i]  - l.start + 1
                s = e - fraglength + 1
            end
            if str == "-"
                s, e = fw - e + 1, fw - s + 1
            end
            w = P.countfun(V[data_entry, i])
            P.incfun(H, k, s, e,  w) 
        end
    end
    H
end

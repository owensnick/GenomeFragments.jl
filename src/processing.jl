
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




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

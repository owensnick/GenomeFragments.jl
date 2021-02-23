module GenomeFragments

## This is a module for handling high throughput sequencing reads

using GenomicFeatures, DataStructures, BioAlignments, Mmap, UnalignedVectors, XAM, IterTools


export FragMatrix, FragMatrixPair, FragMatrixSingle, get_frags, chrom_lt, load_frag_matrix, totalfrags, streambam, qualityfilt, filterpairfirst, filterpairfirstpp, 
       barcode, ishighquality, isunique,
       inc_meta_frag!, inc_meta_cut!, inc_meta_mid!,
       inc_heat_frag!, inc_heat_cut!, inc_heat_mid!,
       inc_heat_atac_cut!, inc_heat_atac_frag!, inc_meta_atac_frag!,
       inc_meta_mid_width, inc_heat_mid_width,
       inc_heat_atac_cut_width, inc_meta_atac_cut_width,
       readsintersecting, fragregion


include("fragmatrix.jl")
include("buildfm.jl")
include("processing.jl")


end # module

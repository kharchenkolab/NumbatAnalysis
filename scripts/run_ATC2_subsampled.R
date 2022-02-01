library(numbat)

out = numbat_subclone(
    count_mat_ATC2,
    ref_hca,
    df_allele_ATC2,
    gtf_hg38,
    genetic_map_hg38,
    min_cells = 20,
    t = 1e-3,
    ncores = 20,
    init_k = 3,
    max_iter = 2,
    plot = TRUE,
    out_dir = '~/paper_data/numbat_out/ATC2_subsampled'
)

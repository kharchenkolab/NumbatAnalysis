library(numbat)

out = numbat_subclone(
    count_mat_ATC2,
    ref_hca,
    df_allele_ATC2,
    gtf_transcript,
    genetic_map_hg38,
    min_cells = 20,
    t = 1e-4,
    alpha = 1e-5,
    ncores = 20,
    init_k = 3,
    max_cost = 150,
    min_LLR = 30,
    max_iter = 2,
    plot = TRUE,
    out_dir = '~/results/test'
)



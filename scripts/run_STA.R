library(Matrix)
library(glue)
home_dir = '/d0-bayes/home/tenggao'
# home_dir = '{home_dir}'
devtools::load_all(glue('{home_dir}/numbat'))

sample = 'NCI-N87'

count_mat = readRDS(glue('{home_dir}/paper_data/processed/NCI-N87_counts.rds'))
df = fread(glue('{home_dir}/external/STA/NCI-N87_allele_counts.tsv.gz')) %>% filter(cell %in% colnames(count_mat))
# ref_internal = readRDS(glue('{home_dir}/external/gastric/ref_internal.rds'))
# ref_internal = ref_internal[,'epithelial',drop=F]

bulk = get_bulk(
    count_mat,
    ref_hca,
    df,
    gtf_hg38,
    genetic_map_hg38
)

segs_loh = bulk %>% detect_loh(t = 1e-4)

segs_loh %>% fwrite(glue('{home_dir}/paper_data/numbat_out/{sample}_new/segs_loh.tsv'), sep = '\t')

out = run_numbat(
        count_mat,
        ref_hca,
        df,
        gtf_hg38,
        genetic_map_hg38,
        t = 1e-5,
        ncores = 20,
        min_cells = 50,
        ncores_nni = 20,
        multi_allelic = TRUE,
        segs_loh = segs_loh,
        out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_new')
    )
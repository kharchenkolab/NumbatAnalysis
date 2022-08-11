# library(numbat)
library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(Matrix)
library(magrittr)
library(parallel)

# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

cell_annot = fread(glue('{home_dir}/paper_data/cell_annotations/cell_annot_WASHU_march.tsv'), sep = '\t') %>%
    split(.$sample)

patient = '27522'
samples = c('27522_Primary', '27522_Remission', '27522_Relapse_1', '27522_Relapse_2')

count_mat = c()
df = c()

for (sample in samples) {
    count_mat[[sample]] = readRDS(glue('{home_dir}/paper_data/processed/{sample}_counts.rds'))
    df[[sample]] = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')
}

count_mat_combined = count_mat[samples] %>% Reduce('cbind', .)
depths = count_mat_combined %>% colSums

ref_patient = aggregate_counts(
    count_mat_combined,
    cell_annot[samples] %>% bind_rows %>% 
        filter(cell_type == 'B') %>%
        mutate(group = cell_type)
)

ref_patient = ref_patient[,'B',drop=FALSE]

tumor_cells = cell_annot[samples] %>% 
    bind_rows() %>% 
    filter(cell_type == 'Plasma') %>%
    mutate(depth = depths[cell]) %>%
    filter(depth > 1500) %>%
    pull(cell)

tumor_cells = intersect(tumor_cells, colnames(count_mat_combined))

out = run_numbat(
    count_mat_combined %>% extract(,tumor_cells),
    ref_patient,
    df[samples] %>% bind_rows %>% filter(cell %in% tumor_cells),
    gtf_hg38,
    genetic_map_hg38,
    min_cells = 50,
    t = 1e-5,
    ncores = 30,
    max_entropy = 0.5,
    multi_allelic = TRUE,
    out_dir = glue('{home_dir}/paper_data/numbat_out/{patient}')
)

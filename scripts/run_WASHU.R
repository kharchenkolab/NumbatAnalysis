library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(Matrix)
library(magrittr)
library(copykat)
library(parallel)

# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

cell_annot = fread(glue('{home_dir}/paper_data/cell_annotations/cell_annot_WASHU.tsv')) %>% 
    mutate(sample_id = str_replace(sample_id, '-', '_')) %>%
    mutate(sample_id = ifelse(sample_id == '57075_Pre_transplant', '57075_Primary', sample_id)) %>%
    mutate(cell = paste0(sample_id, '_', barcode)) %>% 
    split(.$sample_id)

con = readRDS(glue('{home_dir}/paper_data/conos_objects/conos_WASHU.rds'))

count_mat = c()
df = c()

samples = c(
    '58408_SMM', '58408_Primary', 
    '47491_SMM', '47491_Primary',
    '27522_Primary', '27522_Relapse_2',
    '37692_Primary', '59114_Relapse_1'
    )

for (sample in samples) {
    count_mat[[sample]] = t(con$samples[[sample]]$misc$rawCounts)
    cells = intersect(cell_annot[[sample]]$cell, colnames(count_mat[[sample]]))
    cell_annot[[sample]] = cell_annot[[sample]] %>% filter(cell %in% cells)
    count_mat[[sample]] = count_mat[[sample]][,cells]
    df[[sample]] = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t') %>%
        filter(cell %in% cells)
}

ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')

for (sample in samples) {

    n_cells = ncol(count_mat[[sample]])

    ncores_nni = case_when(
        n_cells < 1000 ~ 6,
        n_cells < 2000 ~ 10,
        n_cells < 4000 ~ 20,
        TRUE ~ 30
    )
        
    tryCatch(
        expr = {
            out = run_numbat(
                count_mat[[sample]],
                ref_hca[,ref_types],
                df[[sample]],
                gtf_hg38,
                genetic_map_hg38,
                min_cells = 50,
                t = 1e-5,
                ncores = 40,
                ncores_nni = ncores_nni,
                max_entropy = 0.5,
                out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_new')
            )
        },
        error = function(e){ 
            cat(glue('Error when running {sample}'))
            cat(e)
        }
    )
}

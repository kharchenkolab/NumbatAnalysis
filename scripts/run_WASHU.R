library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(Matrix)
library(parallel)

home_dir = '/d0-bayes/home/tenggao'
# home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

samples = c(
    '47491_SMM', '47491_Primary',
    '58408_SMM', '58408_Primary', 
    '27522_Primary', '27522_Relapse_2',
    '37692_Primary', '59114_Relapse_1')

ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')

for (sample in samples) {

    count_mat = readRDS(glue('{home_dir}/paper_data/processed/{sample}_counts.rds'))
    df = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')

    n_cells = ncol(count_mat)

    ncores_nni = case_when(
        n_cells < 1000 ~ 6,
        n_cells < 2000 ~ 10,
        n_cells < 4000 ~ 20,
        TRUE ~ 20
    )
        
    tryCatch(
        expr = {
            out = run_numbat(
                count_mat,
                ref_hca[,ref_types],
                df,
                gtf_hg38,
                genetic_map_hg38,
                min_cells = 50,
                t = 1e-5,
                ncores = 20,
                multi_allelic = TRUE,
                ncores_nni = ncores_nni,
                max_entropy = 0.5,
                out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_test')
            )
        },
        error = function(e){ 
            cat(glue('Error when running {sample}'))
            cat(e)
        }
    )
}

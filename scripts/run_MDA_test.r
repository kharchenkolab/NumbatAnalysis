library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(parallel)
home_dir = '/d0-bayes/home/tenggao'
# home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

# samples = c('TNBC5', 'ATC2_subsampled', 'TNBC1', 'ATC1')
samples = c('ATC5', 'TNBC2', 'TNBC3', 'TNBC4', 'ATC2', 'ATC3', 'ATC4', 'DCIS1')

## Run Numbat
for (sample in samples) {

    count_mat = readRDS(glue('{home_dir}/paper_data/processed/{sample}_counts.rds'))
    df = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')

    if (sample == 'TNBC1') {
        use_loh = TRUE 
        diploid_chroms = c(13,14,19)
    } else {
        use_loh = NULL
        diploid_chroms = NULL
    }

    n_cells = ncol(count_mat)

    ncores_nni = case_when(
        n_cells < 1000 ~ 6,
        n_cells < 2000 ~ 10,
        n_cells < 4000 ~ 20,
        TRUE ~ 30
    )

    tryCatch(
        expr = {
            out = run_numbat(
                count_mat,
                ref_hca,
                df,
                gtf_hg38,
                genetic_map_hg38,
                min_cells = 50,
                ncores = 40,
                ncores_nni = ncores_nni,
                t = 1e-5,
                multi_allelic = TRUE,
                use_loh = use_loh,
                diploid_chroms = diploid_chroms,
                out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_test')
            )
        },
        error = function(e) { 
            print(glue('Error when running {sample}'))
            print(e)
        }
    )
    
}

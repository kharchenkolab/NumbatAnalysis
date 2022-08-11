library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(parallel)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
# devtools::load_all(glue('{home_dir}/numbat'))
library(numbat)

# samples = c('TNBC5', 'ATC2_subsampled', 'TNBC1', 'ATC1', 'ATC5', 'TNBC2', 'TNBC3', 'TNBC4', 'ATC2', 'ATC3', 'ATC4', 'DCIS1')
samples = c('DCIS1')

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

    tryCatch(
        expr = {
            out = run_numbat(
                count_mat,
                ref_hca,
                df,
                gtf_hg38,
                genetic_map_hg38,
                min_cells = 50,
                ncores = 30,
                t = 1e-5,
                multi_allelic = TRUE,
                use_loh = use_loh,
                diploid_chroms = diploid_chroms,
                out_dir = glue('~/paper_data/numbat_out/{sample}')
            )
        },
        error = function(e) { 
            print(glue('Error when running {sample}'))
            print(e)
        }
    )
    
}

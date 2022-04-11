library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(purrr)
library(parallel)
library(vcfR)
library(pagoda2)
library(conos)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

# expression data
con_TNBC = readRDS(glue('{home_dir}/paper_data/conos_objects/conos_TNBC.rds'))
con_ATC = readRDS(glue('{home_dir}/paper_data/conos_objects/conos_ATC.rds'))

samples = c('TNBC5', 'TNBC1', 'ATC1', 'ATC3', 'ATC2_subsampled', 'TNBC2', 'TNBC3', 'TNBC4', 'ATC4', 'ATC5', 'ATC2')

count_mat = c()
df = c()

for (sample in samples) {
    if (sample == 'ATC2_subsampled') {
        count_mat[[sample]] = count_mat_ATC2
        df[[sample]] = df_allele_ATC2
    } else {
        if (str_detect(sample, 'ATC')) {
            con = con_ATC
        } else {
            con = con_TNBC
        }
        count_mat[[sample]] = t(con$samples[[sample]]$misc$rawCounts)
        df[[sample]] = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')
    }
}

## Run Numbat
for (sample in samples) {

    if (sample == 'TNBC1') {
        use_loh = TRUE 
        diploid_chroms = c(13,14,19)
    } else {
        use_loh = NULL
        diploid_chroms = NULL
    }

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
                ref_hca,
                df[[sample]],
                gtf_hg38,
                genetic_map_hg38,
                min_cells = 50,
                ncores = 40,
                ncores_nni = ncores_nni,
                init_k = 3,
                t = 1e-3,
                multi_allelic = TRUE,
                use_loh = use_loh,
                diploid_chroms = diploid_chroms,
                max_iter = 2,
                out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_new2')
            )
        },
        error = function(e) { 
            print(glue('Error when running {sample}'))
            print(e)
        }
    )
    
}

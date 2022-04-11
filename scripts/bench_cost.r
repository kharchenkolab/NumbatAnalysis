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

samples = c('ATC2', 'ATC1', 'TNBC5', 'TNBC1')

count_mat = c()
df = c()

for (sample in samples) {
    if (sample == 'ATC2') {
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

    for (tau in seq(0.1,0.5,0.1)) {

        tryCatch(
            expr = {
                out = run_numbat(
                    count_mat[[sample]],
                    ref_hca,
                    df[[sample]],
                    gtf_hg38,
                    genetic_map_hg38,
                    min_cells = 50,
                    ncores = 30,
                    ncores_nni = 15,
                    t = 1e-3,
                    multi_allelic = TRUE,
                    use_loh = use_loh,
                    diploid_chroms = diploid_chroms,
                    tau = tau,
                    out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_new_tau_{tau}')
                )
            },
            error = function(e) { 
                print(glue('Error when running {sample}'))
                print(e)
            }
        )

    }
    
}

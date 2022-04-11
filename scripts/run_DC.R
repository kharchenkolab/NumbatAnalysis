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
library(IRdisplay)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

# expression data
con = readRDS(glue('{home_dir}/paper_data/conos_objects/conos_BC.rds'))

samples = c('DCIS1')

count_mat = c()
df = c()

for (sample in samples) {
    count_mat[[sample]] = t(con$samples[[sample]]$misc$rawCounts)
    df[[sample]] = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')
}

## Run Numbat
for (sample in samples) {
    
    tryCatch(
        expr = {
            out = run_numbat(
                count_mat[[sample]],
                ref_hca,
                df[[sample]],
                gtf_hg38,
                genetic_map_hg38,
                min_cells = 50,
                t = 1e-3,
                ncores = 30,
                init_k = 3,
                out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_new')
            )
        },
        error = function(e) { 
            print(glue('Error when running {sample}'))
            print(e)
        }
    )
    
}

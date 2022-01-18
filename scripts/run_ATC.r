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
devtools::load_all('~/Numbat')
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'

# expression data
con = readRDS(glue('{home_dir}/external/MDA/conos_ATC.rds'))

samples = paste0('ATC', 1)

count_mat = c()
df = c()

for (sample in samples) {
    count_mat[[sample]] = as.matrix(t(con$samples[[sample]]$misc$rawCounts))
    df[[sample]] = fread(glue('{home_dir}/external/MDA/{sample}_allele_counts.tsv'), sep = '\t')
}

## Run Numbat
for (sample in samples) {
    
    tryCatch(
        expr = {
            out = numbat_subclone(
                count_mat[[sample]],
                ref_hca,
                df[[sample]],
                gtf_transcript,
                genetic_map_hg38,
                min_cells = 50,
                t = 1e-5,
                ncores = 30,
                init_k = 3,
                max_entropy = 0.6,
                out_dir = glue('{home_dir}/results/MDA/{sample}_final')
            )
        },
        error = function(e) { 
            print(glue('Error when running {sample}'))
            print(e)
        }
    )
    
}

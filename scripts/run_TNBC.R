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
library(numbat)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'

# expression data
con = readRDS(glue('{home_dir}/paper_data/conos_objects/conos_TNBC.rds'))

samples = paste0('TNBC', 1)

count_mat = c()
df = c()

for (sample in samples) {
    count_mat[[sample]] = as.matrix(t(con$samples[[sample]]$misc$rawCounts))
    df[[sample]] = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')
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
    
    tryCatch(
        expr = {
            out = numbat_subclone(
                count_mat[[sample]],
                ref_hca,
                df[[sample]],
                gtf_transcript,
                genetic_map_hg38,
                min_cells = 50,
                ncores = 30,
                init_k = 3,
                t = 1e-3,
                multi_allelic = T,
                use_loh = use_loh,
                diploid_chroms = diploid_chroms,
                max_iter = 2,
                out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}')
            )
        },
        error = function(e) { 
            print(glue('Error when running {sample}'))
            print(e)
        }
    )
    
}

library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(purrr)
library(parallel)
library(vcfR)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

ncores = 10

# samples = c('ATC2_subsampled', 'ATC1', 'TNBC5', 'TNBC1', 'NCI-N87', 
#     '47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '37692_Primary', '58408_Primary')

samples = c('NCI-N87')

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

    if (str_detect(sample, 'TNBC|ATC|NCI')) {
        ref_types = colnames(ref_hca)
    } else {
        ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')
    }

    if (str_detect(sample, 'NCI')) {
        segs_loh = fread(glue('{home_dir}/paper_data/numbat_out/NCI-N87_new/segs_loh.tsv')) %>% relevel_chrom()
    } else {
        segs_loh = NULL
    }

    n_cells = ncol(count_mat)
    ncores_nni = case_when(n_cells < 1000 ~ 6, n_cells < 2000 ~ 10, n_cells < 4000 ~ 20, TRUE ~ 30)

    for (tau in seq(0.1,0.5,0.1)) {

        tryCatch(
            expr = {
                out = run_numbat(
                    count_mat,
                    ref_hca[,ref_types],
                    df,
                    gtf_hg38,
                    genetic_map_hg38,
                    min_cells = 50,
                    ncores = ncores,
                    ncores_nni = min(ncores_nni, ncores),
                    t = 1e-5,
                    multi_allelic = TRUE,
                    use_loh = use_loh,
                    segs_loh = segs_loh,
                    diploid_chroms = diploid_chroms,
                    tau = tau,
                    out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_test_tau_{tau}')
                )
            },
            error = function(e) { 
                print(glue('Error when running {sample}'))
                print(e)
            }
        )

    }
    
}

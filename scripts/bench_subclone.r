library(dplyr)
library(data.table)
library(stringr)
library(glue)
library(Matrix)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

con_washu = readRDS(glue("{home_dir}/paper_data/conos_objects/conos_WASHU.rds"))

cell_annot = fread(glue('{home_dir}/paper_data/cell_annotations/cell_annot_WASHU_march.tsv'), sep = '\t') %>%
    split(.$sample)

samples = c('47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '37692_Primary', '58408_Primary')
ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')

for (sample in samples) {

    count_mat = readRDS(glue('{home_dir}/paper_data/processed/{sample}_counts.rds'))
    df = fread(glue('{home_dir}/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t')

    for (frac in c(0.1, 0.3, 0.5, 0.7, 0.9)) {

        for (seed in c(0)) {
            
            set.seed(seed)

            tumor_cells = cell_annot[[sample]] %>% filter(annot == 'T') %>% pull(cell)
            
            n_cells = length(tumor_cells)
            
            cells_tumor = cell_annot[[sample]] %>% 
                filter(annot == 'T') %>%
                pull(cell) %>% sample(n_cells * frac)

            cells_normal = cell_annot[[sample]] %>%
                filter(annot == 'N') %>% 
                pull(cell) %>% sample(n_cells * (1-frac))

            cells = c(cells_normal, cells_tumor)

             tryCatch(
                expr = {
                    out = run_numbat(
                        count_mat[,cells],
                        ref_hca[,ref_types],
                        df %>% filter(cell %in% cells),
                        gtf_hg38,
                        genetic_map_hg38,
                        min_cells = 0,
                        t = 1e-5,
                        ncores = 30,
                        multi_allelic = TRUE,
                        ncores_nni = 8,
                        out_dir = glue('{home_dir}/paper_data/numbat_out/{sample}_test_{frac}_{seed}')
                    )
                },
                error = function(e) { 
                    print(glue('Error when running {sample}'))
                    print(e)
                }
            )

        }

    }

}


library(dplyr)
library(data.table)
library(stringr)
library(glue)
library(Matrix)
# home_dir = '/d0-bayes/home/tenggao'
home_dir = '/home/tenggao'

con_washu = readRDS(glue("{home_dir}/paper_data/conos_objects/conos_WASHU.rds"))

cell_annot = fread(glue('{home_dir}/paper_data/cell_annotations/cell_annot_WASHU_march.tsv'), sep = '\t') %>%
    split(.$sample)

count_mat = c()
df = c()

samples = c('47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '37692_Primary', '58408_Primary')

for (sample in samples) {
    count_mat[[sample]] = t(con_washu$samples[[sample]]$misc$rawCounts)
    cells = intersect(cell_annot[[sample]]$cell, colnames(count_mat[[sample]]))
    cell_annot[[sample]] = cell_annot[[sample]] %>% filter(cell %in% cells)
    count_mat[[sample]] = count_mat[[sample]][,cells]
}

ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')

##### Copykat #####
library(copykat)

load('~/numbat/data/ref_hca_counts.rda')
ref_counts = ref_hca_counts[,rep(ref_types,2),drop=F]
colnames(ref_counts) = c(paste0(ref_types, '_', 1), paste0(ref_types, '_', 2))

for (sample in samples) {

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

            genes_common = rownames(count_mat[[sample]]) %>% intersect(rownames(ref_counts))

            out_dir = glue('{home_dir}/paper_data/copykat_out/{sample}_{frac}_{seed}')
            dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
            setwd(out_dir)

            copykat.test <- copykat(
                rawmat = as.matrix(cbind(
                    ref_counts[genes_common,],
                    count_mat[[sample]][genes_common,cells]
                )), 
                id.type="S", 
                ngene.chr=0, 
                win.size=25,
                KS.cut=0.1,
                sam.name=sample, 
                distance="euclidean", 
                norm.cell.names = colnames(ref_counts),
                n.cores=30
            )

        }
    }
}


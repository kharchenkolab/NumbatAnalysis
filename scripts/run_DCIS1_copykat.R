library(numbat)
library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(Matrix)
library(copykat)
con = readRDS('~/paper_data/conos_objects/conos_BC.rds')

ref = 'hca'
ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plamsa', 'Dendritic')
ref_counts = ref_hca_counts[,rep(ref_types,2),drop=F]

colnames(ref_counts) = c(paste0(ref_types, '_', 1), paste0(ref_types, '_', 2))

print('Running copykat')

sample = 'DCIS1'

count_mat = list()
count_mat[[sample]] = as.matrix(t(con$samples[[sample]]$misc$rawCounts))

genes_common = rownames(count_mat[[sample]]) %>% intersect(rownames(ref_counts))

out_dir = glue('~/paper_data/copykat_out/{sample}')

dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
setwd(out_dir)

copykat.test <- copykat(
    rawmat = cbind(
        ref_counts[genes_common,],
        count_mat[[sample]][genes_common,]
    ), 
    id.type="S", 
    ngene.chr=0, 
    sam.name=sample, 
    distance="euclidean", 
    norm.cell.names = colnames(ref_counts),
    n.cores=30
)

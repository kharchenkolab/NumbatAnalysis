library(dplyr)
library(data.table)
library(stringr)
library(glue)
library(pagoda2)
library(conos)

con = readRDS(glue('~/paper_data/conos_objects/conos_ATC.rds'))

cell_annot = fread(glue('~/paper_data/cell_annotations/cell_annot_MDA.tsv')) %>% 
    mutate(annot = copykat.pred) %>% 
    split(.$sample)

set.seed(0)
tumor_cells = cell_annot[['ATC2']] %>% filter(annot == 'T') %>% pull(cell)
normal_cells = cell_annot[['ATC2']] %>% filter(annot == 'N') %>% pull(cell) %>% sample(50)
cells = c(tumor_cells, normal_cells)

count_mat = as.matrix(t(con$samples[['ATC2']]$misc$rawCounts))
df = fread(glue('~/paper_data/processed/ATC2_allele_counts.tsv.gz'), sep = '\t')

count_mat_ATC2 = count_mat[,cells]
df_allele_ATC2 = df %>% filter(cell %in% cells)

# use_data(df_allele_ATC2, overwrite = T)
# use_data(count_mat_ATC2, overwrite = T)

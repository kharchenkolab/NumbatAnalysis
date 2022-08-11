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

count_mat = c()
df = c()

samples = c('47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '37692_Primary', '58408_Primary')

for (sample in samples) {
    count_mat[[sample]] = t(con_washu$samples[[sample]]$misc$rawCounts)
    cells = intersect(cell_annot[[sample]]$cell, colnames(count_mat[[sample]]))
    cell_annot[[sample]] = cell_annot[[sample]] %>% filter(cell %in% cells)
    count_mat[[sample]] = count_mat[[sample]][,cells]
}

load('~/numbat/data/ref_hca_counts.rda')

ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')

# averaging the expression profile so that inferCNV outputs consensus CNV calls for the cell mixture
# this is because inferCNV HMM outputs one set of CNV calls for each observation group (i.e. cell type)
ref_counts = rowSums(ref_hca_counts) %>% 
    as.matrix %>%
    magrittr::set_colnames('Plasma')

##### InferCNV #####

ref_types = c('Plasma')

ref_counts = ref_counts[,rep(ref_types,2),drop=F]
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

            # preparing input
            out_dir = glue('{home_dir}/paper_data/infercnv_out/{sample}_{frac}_{seed}')
            dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
            annot_file = glue('{out_dir}/{sample}_annot.tsv')

            ref_groups = paste0(str_remove(colnames(ref_counts), '_1|_2'), '_ref')

            cell_annot[[sample]] %>% 
            filter(cell %in% cells) %>%
            mutate(group = 'Plasma') %>%
            bind_rows(
                data.frame(
                    cell = colnames(ref_counts),
                    group = ref_groups
                )
            ) %>%
            select(cell, group) %>%
            fwrite(annot_file, sep = '\t', col.names = F)

            genes_common = rownames(count_mat[[sample]]) %>% intersect(rownames(ref_counts))

            # run
            infercnv_obj = infercnv::CreateInfercnvObject(
                raw_counts_matrix = cbind(
                        ref_counts[genes_common,,drop=F],
                        count_mat[[sample]][genes_common,cells]
                ),
                annotations_file = annot_file,
                delim="\t",
                gene_order_file="/home/tenggao/inferCNV/hg38_gene_order.txt",
                ref_group_names = unique(ref_groups)
            )

            infercnv_obj = infercnv::run(
                infercnv_obj,
                num_threads = 30,
                cutoff=0.1, 
                window_length = 101,
                out_dir=out_dir,
                cluster_by_groups=F, 
                denoise=TRUE,
                HMM=TRUE,
                HMM_report_by = "consensus",
                resume_mode=FALSE
            )

        }
    }
}

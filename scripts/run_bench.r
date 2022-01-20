library(numbat)
library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(Matrix)
library(magrittr)
library(copykat)
library(IRdisplay)
library(parallel)

cell_annot = fread('~/paper_data/cell_annotations/cell_annot_WASHU.tsv') %>% 
    mutate(sample_id = str_replace(sample_id, '-', '_')) %>%
    mutate(sample_id = ifelse(sample_id == '57075_Pre_transplant', '57075_Primary', sample_id)) %>%
    mutate(cell = paste0(sample_id, '_', barcode)) %>% 
    split(.$sample_id)

con = readRDS('~/paper_data/conos_objects/conos_WASHU.rds')

# samples = c('37692_Primary', '47491_SMM', '47491_Primary', '27522_Primary', '27522_Relapse_2', '58408_SMM', '58408_Primary', )
# samples = c('58408_SMM', '58408_Primary')
samples = c('59114_Relapse_1')

count_mat = c()
df = c()

for (sample in samples) {
    count_mat[[sample]] = as.matrix(t(con$samples[[sample]]$misc$rawCounts))
    cells = intersect(cell_annot[[sample]]$cell, colnames(count_mat[[sample]]))
    cell_annot[[sample]] = cell_annot[[sample]] %>% filter(cell %in% cells)
    count_mat[[sample]] = count_mat[[sample]][,cells]
    df[[sample]] = fread(glue('~/paper_data/processed/{sample}_allele_counts.tsv.gz'), sep = '\t') %>%
        filter(cell %in% cells)
}

ref = 'hca'
ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plamsa', 'Dendritic')
ref_exp = ref_hca[,str_replace(ref_types, 'Plamsa', 'Plasma')]
ref_counts = ref_hca_counts[,rep(ref_types,2),drop=F]

colnames(ref_counts) = c(paste0(ref_types, '_', 1), paste0(ref_types, '_', 2))

print('Running copykat')
for (sample in samples) {

    genes_common = rownames(count_mat[[sample]]) %>% intersect(rownames(ref_counts))

    for (smooth in c(TRUE)) {
        
        if (smooth) {
            win.size = 25
            KS.cut = 0.1
            out_dir = glue('~/copykat_results/{sample}/{ref}/smooth')
        } else {
            win.size = 1
            KS.cut = 0
            out_dir = glue('~/copykat_results/{sample}/{ref}/no_smooth')
        }

        dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
        setwd(out_dir)

        copykat.test <- copykat(
            rawmat = cbind(
                ref_counts[genes_common,],
                count_mat[[sample]][genes_common,]
            ), 
            id.type="S", 
            ngene.chr=0, 
            win.size=win.size,
            KS.cut=KS.cut,
            sam.name=sample, 
            distance="euclidean", 
            norm.cell.names = colnames(ref_counts),
            n.cores=30
        )

    }
    
}

print('Running inferCNV')
for (sample in samples) {

    for (smooth in c(TRUE)) {

        if (smooth) {
            window_length = 101
            out_dir = glue('~/inferCNV/{sample}/{ref}/smooth')
        } else {
            window_length = 1
            out_dir = glue('~/inferCNV/{sample}/{ref}/no_smooth')
        }
        
        # preparing input
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        annot_file = glue('{out_dir}/{sample}_annot.tsv')

        ref_groups = paste0(str_remove(colnames(ref_counts), '_1|_2'), '_ref')

        cell_annot[[sample]] %>% 
        mutate(group = cell_type) %>%
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
                    count_mat[[sample]][genes_common,]
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
            window_length = window_length,
            out_dir=out_dir,
            cluster_by_groups=F, 
            denoise=TRUE,
            HMM=TRUE,
            HMM_report_by = "consensus",
            resume_mode=FALSE
        )
    }
}

print('Running Numbat')
for (sample in samples) {

    bulk_tumor = fread(glue('~/paper_data/numbat_out/{sample}/bulk_tumor.tsv')) %>% 
        mutate(CHROM = factor(CHROM)) %>%
        as.data.frame()
    segs_consensus = fread(glue('~/paper_data/purple_out/{sample}/segs_consensus_dna.tsv')) %>% 
        mutate(CHROM = factor(CHROM)) %>%
        as.data.frame()

    out_dir = glue('~/paper_out/cnv_benchmark/{sample}')
    dir.create(out_dir, showWarnings = FALSE)

    exp_post = get_exp_post(
        segs_consensus,
        count_mat[[sample]],
        gtf_transcript,
        ref_exp,
        ncores = 64)$exp_post

    allele_post = get_allele_post(
        bulk_tumor,
        segs_consensus,
        df[[sample]]
    )

    joint_post = get_joint_post(
        exp_post,
        allele_post,
        segs_consensus)

    fwrite(exp_post, glue('{out_dir}/exp_post.tsv'), sep = '\t')
    fwrite(allele_post, glue('{out_dir}/allele_post.tsv'), sep = '\t')
    fwrite(joint_post, glue('{out_dir}/joint_post.tsv'), sep = '\t')
}

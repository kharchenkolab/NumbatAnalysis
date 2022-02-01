library(HoneyBADGER)
require(biomaRt) 
library(data.table)
library(dplyr)
library(glue)
library(rjags)
library(parallel)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## in order to map SNPs to genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

mart.obj = NULL
while(is.null(mart.obj)) {
    message('trying biomart')
    try(
        mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl')
    )
}

annot = fread('~/external/WASHU/sample_barcode_cell_type_annotation.txt')

cell_annot = annot %>% 
    mutate(sample_id = str_replace(sample_id, '-', '_')) %>%
    mutate(sample_id = ifelse(sample_id == '57075_Pre_transplant', '57075_Primary', sample_id)) %>%
    mutate(cell = paste0(sample_id, '_', barcode)) %>% 
    split(.$sample_id)

ref_types = c('NK', 'Macrophage', 'CD4+T', 'CD8+T', 'Myeloid', 'Monocyte', 'B', 'Plasma', 'Dendritic')

con = readRDS('~/external/WASHU/con.rds')

# samples = c('37692_Primary', '47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '58408_Primary')
samples = c('59114_Relapse_1', '58408_Primary')

count_mat = c()
df = c()

for (sample in samples) {
    count_mat[[sample]] = as.matrix(t(con$samples[[sample]]$misc$rawCounts))
    cells = intersect(cell_annot[[sample]]$cell, colnames(count_mat[[sample]]))
    cell_annot[[sample]] = cell_annot[[sample]] %>% filter(cell %in% cells)
    count_mat[[sample]] = count_mat[[sample]][,cells]
    df[[sample]] = fread(glue('~/external/WASHU/{sample}_allele_counts.tsv'), sep = '\t') %>%
        filter(cell %in% cells)
}

hb = mclapply(
    mc.cores = length(samples),
    samples,
    function(sample) {
        readRDS(glue('~/results/benchmark/hb_{sample}.rds'))
    }
) %>% setNames(samples)
    

print('Running HB')
for (sample in samples) {

    message(sample)

    segs_consensus = fread(glue('~/results/benchmark/{sample}/segs_consensus_dna.tsv')) %>% 
        mutate(CHROM = factor(CHROM)) %>%
        as.data.frame() %>%
        filter(cnv_state != 'neu')

    cells = colnames(hb[[sample]]$r.maf)

    res_all = data.frame()
    
    for (i in 1:nrow(segs_consensus)) {

        seg = segs_consensus[i,]

        message(paste0('testing seg ', seg$seg))

        region = GenomicRanges::GRanges(
            seqnames = paste0('chr', seg$CHROM),
            IRanges::IRanges(
                start = seg$seg_start,
                end = seg$seg_end
            )
        )

        res_list = mclapply(
            mc.cores = 90,
            split(cells, ceiling(seq_along(cells)/10)),
            function(cells_subset) {

                out = hb[[sample]]$calcCombCnvProb(
                    region = region,
                    r.sub = hb[[sample]]$r.maf[,cells_subset],
                    n.sc = hb[[sample]]$n.sc[,cells_subset],
                    verbose = F, 
                    quiet = T,
                    filter = T
                )

                data.frame(
                    cell = names(out$`posterior probability of amplification`),
                    p_amp = out$`posterior probability of amplification`,
                    p_del = out$`posterior probability of deletion`
                )
            }
        )

        bad = sapply(res_list, inherits, what = "try-error")

        if (any(bad)) {
            message(glue('job {paste(which(bad), collapse = ",")} failed'))
            message(res_list[bad][[1]])
        }

        res = res_list[!bad] %>%
            bind_rows() %>%
            mutate(
                CHROM = seg$CHROM,
                seg_start = seg$seg_start,
                seg_end = seg$seg_end,
                seg = seg$seg
            )

        res_all = bind_rows(res_all, res)
    }

    fwrite(res_all, glue('~/results/benchmark/{sample}_hb_sc.tsv'), sep = '\t')

}

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

ncores = 40

samples = c('37692_Primary', '47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '58408_Primary')

hb = mclapply(
    mc.cores = length(samples),
    samples,
    function(sample) {
        readRDS(glue('~/paper_data/honeybadger_out/hb_{sample}.rds'))
    }
) %>% setNames(samples)
    

print('Running HB')
for (sample in samples) {

    message(sample)

    segs_consensus = fread(glue('~/paper_data/purple_out/{sample}/segs_consensus_dna.tsv')) %>% 
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
            mc.cores = ncores,
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
            # message(glue('job {paste(which(bad), collapse = ",")} failed'))
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

    fwrite(res_all, glue('~/paper_data/honeybadger_out/{sample}_hb_sc.tsv'), sep = '\t')

}

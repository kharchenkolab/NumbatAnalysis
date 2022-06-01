library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(glue)
library(parallel)
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

run_bench = function(df, cell_annot, segs_dna, sample, n_cells = 500, ratios = seq(0, 1, 0.05), seed = 0, gamma = 20, ncores = 10, replace = FALSE) {

    ### simulate cell mixtures and infer haplotypes ### 
    message('Simulating mixture pseudobulks .. ')

    set.seed(seed)
    nodes = list()
    cells_train = cell_annot %>% sample_frac(0.7) %>% pull(cell)
    cells_test = cell_annot %>% filter(!cell %in% cells_train) %>% pull(cell)

    for (n in n_cells) {
        for (r in ratios) {
            for (s in 1:5) {

                set.seed(s)

                cells_tumor = cell_annot %>% 
                    filter(annot == 'T' & cell %in% cells_train) %>%
                    pull(cell) %>% sample(n * r, replace = replace)

                cells_normal = cell_annot %>%
                    filter(annot == 'N' & cell %in% cells_train) %>% 
                    pull(cell) %>% sample(n * (1-r), replace = replace)

                nodes[[paste0(r, '_', s, '_', n)]] = list(cells = c(cells_tumor, cells_normal), members = s, label = r, size = n)
            }
        }
    }

    bulks_mixture = mclapply(
        nodes,
        mc.cores = ncores,
        function(g) {
            get_allele_bulk(
                df_allele = df %>% filter(cell %in% g$cells),
                genetic_map = genetic_map_hg38,
                lambda = 1
            ) %>%
            mutate(
                seed = g$members,
                ratio = g$label,
                n_cells = g$n_cells
            ) %>%
            annot_consensus(
                segs_dna
            ) %>%
            group_by(CHROM, seg) %>%
            mutate(
                approx_theta_post(pAD, DP, p_s, gamma = gamma, start = 0.1)
            ) %>%
            classify_alleles() %>%
            ungroup()
    })

    ### Scoring cells ### 
    message('Scoring cells .. ')
    cell_scores = mclapply(
        bulks_mixture,
        mc.cores = ncores,
        function(bulk) {
                        
            allele_post_naive = get_allele_post(
                bulk,
                segs_dna %>% mutate(seg = factor(seg_cons)),
                df %>% filter(cell %in% cells_test),
                naive = TRUE
            )

            allele_post_phasing = get_allele_post(
                bulk,
                segs_dna %>% mutate(seg = factor(seg_cons)),
                df %>% filter(cell %in% cells_test),
                naive = FALSE
            )

            allele_post = rbind(
                allele_post_naive %>% mutate(method = 'naive'),
                allele_post_phasing %>% mutate(method = 'phasing')
            ) %>%
            mutate(
                seed = unique(bulk$seed),
                ratio = unique(bulk$ratio)
            )

            return(allele_post)

    }) %>% bind_rows()

    fwrite(
        cell_scores,
        glue('~/paper_data/phasing_benchmark/cell_scores_{sample}.tsv.gz'),
        sep = '\t', nThread = min(10, ncores)
    )

    message('Done!')

}

####### run ########

samples = c('47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '37692_Primary', '58408_Primary')

cell_annot = fread('~/paper_data/cell_annotations/cell_annot_WASHU_march.tsv') %>%
    split(.$sample)

df = c()
segs_dna = c()

for (sample in samples) {

    message(sample)

    for (ratio in c(1)) {

        df = fread(glue('~/paper_data/processed/{sample}_{ratio}_allele_counts.tsv.gz'), nThread = 4)

        segs_dna = fread(glue('~/paper_data/purple_out/{sample}/segs_consensus_dna.tsv')) %>% 
            filter(cnv_state_post %in% c('amp', 'loh', 'del')) %>%
            mutate(cnv_state = cnv_state_post) %>%
            mutate(CHROM = factor(CHROM, 1:22)) %>% 
            as.data.frame()

        run_bench(
            df, 
            cell_annot[[sample]], 
            segs_dna, 
            paste0(sample, '_', ratio),
            n_cells = 500, 
            ratios = seq(0, 1, 0.05), 
            ncores = 10,
            replace = TRUE
        )
    }
}


# ## TNBC4 ##
# cell_annot = fread('~/paper_data/cell_annotations/cell_annot_MDA.tsv', sep = ',') %>% 
#     mutate(annot = copykat.pred) %>% 
#     filter(sample == 'TNBC4')

# segs_dna = fread('~/paper_data/phasing_benchmark/segs_truth_TNBC4.tsv') %>% 
#     mutate(CHROM = factor(CHROM, 1:22)) %>%
#     as.data.frame()

# samples = c(
#     'TNBC4', 'TNBC4_1000G',
#     'TNBC4_0.25', 'TNBC4_0.5'
# )

# df = mclapply(
#     samples,
#     mc.cores = length(samples),
#     function(sample) {
#         fread(glue('~/paper_data/processed/{sample}_N_allele_counts.tsv.gz'))
#     }) %>%
#     setNames(samples)

# for (sample in samples) {

#     message(sample)

#     run_bench(
#         df[[sample]], cell_annot, segs_dna, sample,
#         n_cells = 750, 
#         ratios = c(seq(0, 0.1, 0.01), seq(0.15, 0.5, 0.05)),
#         ncores = 10,
#         replace = FALSE
#     )
# }


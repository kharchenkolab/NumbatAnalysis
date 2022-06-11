library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(glue)
library(parallel)
home_dir = '/home/tenggao'
devtools::load_all(glue('{home_dir}/numbat'))

run_bench = function(df, cell_annot, segs_dna, sample, n_cells, ratios, ls, seed = 0, gamma = 20, ncores = 10, replace = FALSE) {

    ### simulate cell mixtures ### 
    message('Simulating mixture pseudobulks .. ')

    set.seed(seed)
    nodes = list()

    for (n in n_cells) {
        for (r in ratios) {
            for (s in 1:5) {
                    
                set.seed(s)
                
                cells_tumor = cell_annot %>% 
                    filter(annot == 'T') %>%
                    pull(cell) %>% sample(n * r, replace = replace)
            
                cells_normal = cell_annot %>%
                    filter(annot == 'N') %>% 
                    pull(cell) %>% sample(n * (1-r), replace = replace)

                nodes[[paste0(n, '_', r, '_', s)]] = list(cells = c(cells_tumor, cells_normal), members = s, label = r, n_cells = n)

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
            )
    })

    # saveRDS(bulks_mixture, glue('~/paper_data/phasing_benchmark/bulks_mixture_{sample}.rds'))
    fwrite(
        bind_rows(bulks_mixture), 
        glue('~/paper_data/phasing_benchmark/bulks_mixture_{sample}.tsv.gz'),
        sep = '\t', nThread = min(10, ncores)
    )

    ### simulate CNV segments ### 
    message('Simulating CNV segments .. ')

    # function to get random segments from the profile
    sample_segs = function(bulk, l = 500, n = 10) {
        
        if (nrow(bulk) < l + n) {
            print(unname(unique(bulk$CHROM)))
            stop('not enough SNPs')
        }
        
        starts = sample(1:(nrow(bulk)-l), n)
            
        bulk_segs = lapply(
            1:length(starts),
            function(i) {
                start = starts[i]
                bulk[start:(start+l),] %>% mutate(start_random = start, end_random = start+l, id_random = i)
            }
        ) %>%
        bind_rows()
        
        return(bulk_segs)
        
    }

    sample_segs_loc = function(bulk, l = 10e6, n = 10) {

        region_size = max(bulk$POS) - min(bulk$POS)

        if (region_size < l) {
            print(unname(unique(bulk$CHROM)))
            stop('not enough SNPs')
        }

        starts = bulk$POS[sample(1:nrow(bulk[bulk$POS < max(bulk$POS) - l,]), n)]

        bulk_segs = lapply(
            1:length(starts),
            function(i) {
                start = starts[i]
                end = start + l
                
                bulk %>% filter(POS >= start & POS < end) %>%
                    mutate(start_random = start, end_random = end, id_random = i, n_snps = n())
            }
        ) %>%
        bind_rows()

        return(bulk_segs)

    }

    set.seed(seed)
    # different cell mixtures
    segs_mixture = mclapply(
        mc.cores = ncores,
        bulks_mixture,
        function(bulk) {
            # different segs                        
            lapply(
                unique(segs_dna$seg),
                function(seg) {
                    # different segment lengths
                    lapply(
                        ls,
                        function(l) {

                            tryCatch(expr = {

                                n_cells = unique(bulk$n_cells)
                                seed = unique(bulk$seed)

                                segs_mixture = bulk %>%
                                    annot_consensus(segs_dna) %>%
                                    filter(seg == UQ(seg)) %>%
                                    filter(!is.na(AD)) %>% 
                                    sample_segs_loc(n = 10, l = l)
                                                        
                                segs_normal = bulks_mixture[[paste0(n_cells, '_0_', seed)]] %>%
                                    annot_consensus(segs_dna) %>%
                                    filter(seg == UQ(seg)) %>%
                                    filter(!is.na(AD)) %>% 
                                    sample_segs_loc(n = 10, l = l) %>%
                                    mutate(ratio = unique(bulk$ratio)) # don't change

                                segs = bind_rows(
                                        segs_mixture %>% mutate(label = 1),
                                        segs_normal %>% mutate(label = 0)
                                    ) %>%
                                    mutate(MAD = ifelse(AD >= DP - AD, AD, DP-AD))

                                segs = segs %>% mutate(seg = seg, l = l)

                                return(segs)
                                
                            }, 
                            error = function(e) {
                                return(data.frame())
                            })

                        }
                    ) %>% 
                    bind_rows()
                }
            ) %>% 
            bind_rows()
        }
    ) %>% bind_rows()

    fwrite(
        segs_mixture,
        glue('~/paper_data/phasing_benchmark/segs_mixture_{sample}.tsv.gz'),
        sep = '\t', nThread = min(10, ncores)
    )

    ### Scoring segments ### 
    message('Scoring segments .. ')

    theta_neu = segs_mixture %>%
        filter(ratio == 0 & id_random == 1) %>%
        group_by(n_cells) %>%
        summarise(
            theta_mle_naive(MAD, DP),
            .groups = 'drop'
        ) %>%
        {setNames(.$theta_mle, .$n_cells)}

    message(paste0('Using baseline theta = ', signif(theta_neu,2)))

    segs_mixture = segs_mixture %>% mutate(theta_neu = theta_neu[as.character(n_cells)])

    scores = segs_mixture %>%
        mutate(seg = factor(seg)) %>%
        split(list(.$id_random, .$seg, .$seed)) %>%
        mclapply(
            mc.cores = ncores,
            function(bulk_seg) {
                                
                scores_phasing = bulk_seg %>% 
                    group_by(
                        CHROM, seg, cnv_state, 
                        id_random, ratio, n_cells, seed, l, label
                    ) %>%
                    summarise(
                        approx_theta_post(pAD, DP, p_s, gamma = gamma),
                        LLR = calc_allele_LLR(
                            pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], theta_mle, theta_0 = theta_0, gamma = unique(gamma)
                        ),
                        start_random = unique(start_random),
                        end_random = unique(end_random),
                        .groups = 'drop'
                    )
                
                scores_naive = bulk_seg %>%
                    group_by(
                        CHROM, seg, cnv_state, 
                        id_random, ratio, n_cells, seed, l, label
                    ) %>%
                    summarise(
                        theta_mle_naive(MAD, DP),
                        LLR = l_bbinom(MAD, DP, gamma*(0.5+theta_mle), gamma*(0.5-theta_mle)) - 
                            l_bbinom(
                                MAD, DP, gamma*(0.5+unique(theta_neu)),
                                gamma*(0.5-unique(theta_neu))
                            ),
                        start_random = unique(start_random),
                        end_random = unique(end_random),
                        .groups = 'drop'
                    )

                scores = bind_rows(
                    scores_phasing %>% mutate(method = 'phasing'),
                    scores_naive %>% mutate(method = 'naive')
                )
            }
        ) %>%
        bind_rows()

    fwrite(
        scores,
        glue('~/paper_data/phasing_benchmark/scores_{sample}.tsv.gz'),
        sep = '\t', nThread = min(10, ncores)
    )

    message('Done!')
}


####### run ########

## TNBC4 ##
cell_annot = fread('~/paper_data/cell_annotations/cell_annot_MDA.tsv', sep = ',') %>% 
    mutate(annot = copykat.pred) %>% 
    filter(sample == 'TNBC4')

segs_dna = fread('~/paper_data/phasing_benchmark/segs_truth_TNBC4.tsv') %>% 
    mutate(CHROM = factor(CHROM, 1:22)) %>%
    mutate(seg_length = seg_end - seg_start) %>%
    filter(seg_length > 10e6)

samples = c(
    'TNBC4', 'TNBC4_1000G',
    'TNBC4_0.25', 'TNBC4_0.5'
)

for (sample in samples) {

    df = fread(glue('~/paper_data/processed/{sample}_N_allele_counts.tsv.gz'), nThread = 4)

    message(sample)

    run_bench(
        df, 
        cell_annot,
        segs_dna,
        sample,
        n_cells = 500, 
        ratios = c(seq(0, 0.09, 0.01), seq(0.1, 0.5, 0.02)),
        ls = 10e6,
        ncores = 20,
        replace = FALSE
    )

}

# ##MM dataset ##
# samples = c('47491_Primary', '27522_Relapse_2', '59114_Relapse_1', '37692_Primary', '58408_Primary')

# cell_annot = fread('~/paper_data/cell_annotations/cell_annot_WASHU_march.tsv') %>%
#         split(.$sample)

# for (sample in samples) {

#     for (ratio in c(0.25, 0.5, 1)) {

#         df = fread(glue('~/paper_data/processed/{sample}_{ratio}_allele_counts.tsv.gz'), nThread = 4)

#         segs_dna = fread(glue('~/paper_data/purple_out/{sample}/segs_consensus_dna.tsv')) %>% 
#             filter(cnv_state_post %in% c('amp', 'loh', 'del')) %>%
#             mutate(cnv_state = cnv_state_post) %>%
#             mutate(CHROM = factor(CHROM, 1:22)) %>% 
#             as.data.frame() %>%
#             filter(seg_length > 10e6)

#         run_bench(
#             df, 
#             cell_annot[[sample]],
#             segs_dna,
#             paste0(sample, '_', ratio),
#             n_cells = 500, 
#             ratios = seq(0, 1, 0.05), 
#             ls = 10e6,
#             ncores = 20,
#             replace = TRUE
#         )

#     }
        
# }

#' @title grid_fit
#'
#' @description
#' \code{create_smudge_container} creates a smudge container (a list if points and their respective)
#' \code{run_replicate} creates
#' \code{get_centrality} calculates distance of the closest edge to the center
#' \code{test_grid_of_coverages} given minumum and maximum coverage, performs a grid search
#'
#' @export

create_smudge_container <- function(cov, cov_tab, smudge_filtering_threshold){
    smudge_container <- list()
    total_genomic_kmers <- sum(cov_tab[ , 'freq'])
    
    for (Bs in 1:8){
        cov_tab_isoB <- cov_tab[cov_tab[ , 'covB'] > cov * ifelse(Bs == 1, 0, Bs - 0.5) & cov_tab[ , 'covB'] < cov * (Bs + 0.5), ]
        # cov_tab_isoB[, 'Bs'] <- Bs
        cov_tab_isoB[, 'As'] <- round(cov_tab_isoB[, 'covA'] / cov) #these are be individual smudges cutouts given coverages
        cov_tab_isoB[cov_tab_isoB[, 'As'] == 0, 'As'] = 1
        for( As in Bs:(16 - Bs)){
            cov_tab_one_smudge <- cov_tab_isoB[cov_tab_isoB[, 'As'] == As, ]
            if (sum(cov_tab_one_smudge[, 'freq']) / total_genomic_kmers > smudge_filtering_threshold){
                label <- paste0(As, "A", Bs, "B") 
                smudge_container[[label]] <- cov_tab_one_smudge[,-which(names(cov_tab_one_smudge) %in% c('is_error', 'As'))]
            }
        }
    }
    return(smudge_container)
}

#' @export
run_replicate <- function(cov, smudge_tab, smudge_filtering_threshold){
    smudge_container <- create_smudge_container(cov, smudge_tab, smudge_filtering_threshold)
    get_centrality(cov, smudge_container)
}

#' @export
get_centrality <- function(cov, smudge_container){
    ### 1; % of kmer pairs in smudges
    per_smudge_size <- sapply(smudge_container, function(x){ sum(x[, 'freq'])})
    # prop_of_kmers <- sum(per_smudge_size) / total_genomic_kmer_pairs ##

    ### 2. measure of centrality

    As <- sapply(strsplit(names(smudge_container), 'A'), function(x){ as.numeric(x[1])})
    Bs <- sapply(strsplit(names(smudge_container), 'A'), function(x){ as.numeric(substr(x[2], 1, 1))})

    cutoffs_tab <- data.frame(name = names(smudge_container), As = As, Bs = Bs, size = per_smudge_size, centrality = NA)

    for(i in 1:length(smudge_container)){
        test_smudge <- smudge_container[[i]]
        # AAB, As = 2, Bs = 1
        A_range <- c(cov * (As[i] - 0.5), cov * (As[i] + 0.5))
        B_range <- c(cov * (Bs[i] - 0.5), cov * (Bs[i] + 0.5))
        summit <- test_smudge[which.max(test_smudge[, 'freq']), ]
        cutoffs_tab[i,'centrality'] <- sum(min(abs(summit[, 'covB'] - B_range)) + min(abs(summit[, 'covA'] - A_range))) ##
    }

    sum(cutoffs_tab[,'centrality'] * (cutoffs_tab[,'size'] / sum(cutoffs_tab[,'size'])))
}

#' @export
test_grid_of_coverages <- function(smudge_tab, smudge_filtering_threshold, min_to_explore = 10, max_to_explore = 50){
    # min_to_explore = 10
    # max_to_explore = 50 # at least for DToL genomes

    ## grid of coverages spaced by cov 2
    covs_to_test <- seq(min_to_explore + 0.05, max_to_explore + 0.05, by = 2)
    centrality_grid <- sapply(covs_to_test, run_replicate, smudge_tab, smudge_filtering_threshold)
    best_candidate <- covs_to_test[which.max(centrality_grid)]
    ## zooming in on decimals
    covs_to_test_decimals <- seq(best_candidate - 1.9, best_candidate + 1.9, by = 0.1)
    centrality_grid_decimals <- sapply(covs_to_test_decimals, run_replicate, smudge_tab, smudge_filtering_threshold)
    local_maxima <- which(diff(sign(diff(centrality_grid_decimals)))==-2)+1
    best_candidate <- covs_to_test_decimals[which.max(centrality_grid_decimals)]
    ## zooming in on second decimals
    covs_to_test_sec_dec <- seq(best_candidate - 0.09, best_candidate + 0.09, by = 0.01)
    centrality_grid_sec_dec <- sapply(covs_to_test_sec_dec, run_replicate, smudge_tab, smudge_filtering_threshold)
    return(data.frame(cov = c(covs_to_test, covs_to_test_decimals, covs_to_test_sec_dec), centrality = c(centrality_grid, centrality_grid_decimals, centrality_grid_sec_dec)))
    # covs_to_test_sec_dec[which.max(centrality_grid_sec_dec)]
}

#' @title estimate_1n_coverage_1d_subsets
#'
#' @description
#' \code{estimate_1n_coverage_1d_subsets} is used as the first proxy of 1n coverage. The estimate is subsequentially refined by the brightest smudge.
#'
#' @export

estimate_1n_coverage_1d_subsets <- function(.cov_tab){
    total_kmers <- sum(.cov_tab[, 3])
    
    minor_freq_subsets <- list(
        AB_subset = .cov_tab[, 'minor_variant_rel_cov'] > 0.48,
        AAB_subset = .cov_tab[, 'minor_variant_rel_cov'] < 0.3433333 & .cov_tab[, 'minor_variant_rel_cov'] > 0.3233333,
        AAAB_subset = .cov_tab[, 'minor_variant_rel_cov'] < 0.26 & .cov_tab[, 'minor_variant_rel_cov'] > 0.24,
        AAAAB_subset = .cov_tab[, 'minor_variant_rel_cov'] < 0.21 & .cov_tab[, 'minor_variant_rel_cov'] > 0.19,
        AAAAAB_subset = .cov_tab[, 'minor_variant_rel_cov'] < 0.1766667 & .cov_tab[, 'minor_variant_rel_cov'] > 0.1566667
    )

    peak_frame_2 <- get_1d_peaks(.cov_tab[minor_freq_subsets[[1]], ], 3)
    peak_frame_3 <- get_1d_peaks(.cov_tab[minor_freq_subsets[[2]], ], 2)
    peak_frame_4 <- get_1d_peaks(.cov_tab[minor_freq_subsets[[3]], ], 1)
    peak_frame_5 <- get_1d_peaks(.cov_tab[minor_freq_subsets[[4]], ], 1)
    peak_frame_6 <- get_1d_peaks(.cov_tab[minor_freq_subsets[[5]], ], 1)
    peak_frame_6_max_height = which.max(peak_frame_6$height)

    skip_peak_4 = FALSE
    if ( !is.na(peak_frame_2$cov[1]) & !is.na(peak_frame_4$cov[1]) ){
        if ( round(peak_frame_4$cov[1] / peak_frame_2$cov[1]) >= 2 ) { #This is characteristic of allotetraploids which do not have a high AAAB signal. In this case, peak_frame_4 may pick up on the AAAAAABB signal especially if there is high repetitiveness
            skip_peak_4 = TRUE #peak_frame_4 will not be included in the weighted mean, since it will artificially inflate the estimate
        }
    }
    skip_peak_5 = FALSE
    if ( !is.na(peak_frame_5$cov[1]) & !is.na(peak_frame_6$cov[1]) ){
        if (abs(peak_frame_6$cov[1] - peak_frame_5$cov[1]) <= 5 & peak_frame_6$height[1] > peak_frame_5$height[1]) { #This is characteristic of hexaploids in which peak_frame_5 may pick up on the AAAAAB signal
          skip_peak_5 = TRUE #peak_frame_5 will not be included in the weighted mean, since it will artificially inflate the estimate
        }
    }

    if ( !is.na(peak_frame_2$cov[1]) ){
        peak_frame_2$cov <- peak_frame_2$cov / (2 * round(peak_frame_2$cov / peak_frame_2$cov[1]))
    }
    if ( !is.na(peak_frame_3$cov[1]) ){
        peak_frame_3$cov <- peak_frame_3$cov / (3 * round(peak_frame_3$cov / peak_frame_3$cov[1]))
    }
    if ( !is.na(peak_frame_4$cov[1]) ){
        peak_frame_4$cov <- peak_frame_4$cov / (4 * round(peak_frame_4$cov / peak_frame_4$cov[1]))
    }
    if ( !is.na(peak_frame_5$cov[1]) ){
        peak_frame_5$cov <- peak_frame_5$cov / (5 * round(peak_frame_5$cov / peak_frame_5$cov[1]))
    }
    if ( !is.na(peak_frame_6$cov[1]) ){
        peak_frame_6$cov <- peak_frame_6$cov / (6 * round(peak_frame_6$cov / peak_frame_6$cov[peak_frame_6_max_height]))
    }

    if (skip_peak_4) {
    	peak_frame <- rbind(peak_frame_2, peak_frame_3, peak_frame_5, peak_frame_6)
    } else if (skip_peak_5) {
    	peak_frame <- rbind(peak_frame_2, peak_frame_3, peak_frame_4, peak_frame_6)
    } else {
    	peak_frame <- rbind(peak_frame_2, peak_frame_3, peak_frame_4, peak_frame_5, peak_frame_6)
    }

    peak_frame <- peak_frame[is.finite(rowSums(peak_frame)),]
    weighted.mean(peak_frame$cov, peak_frame$height, na.rm = T)
}

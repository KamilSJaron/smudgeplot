#' @title estimate_1n_coverage_1d_subsets
#'
#' @description
#' \code{estimate_1n_coverage_1d_subsets} is a legacy code
#'
#' @export

estimate_1n_coverage_1d_subsets <- function(.total_pair_cov, .minor_variant_rel_cov){
    total_kmers <- length(.minor_variant_rel_cov)
    minor_freq_subsets <- list(
        AB_subset = .minor_variant_rel_cov > 0.47,
        AAB_subset = .minor_variant_rel_cov < 0.35 & .minor_variant_rel_cov > 0.31,
        AAAB_subset = .minor_variant_rel_cov < 0.27 & .minor_variant_rel_cov > 0.23,
        AAAAB_subset = .minor_variant_rel_cov < 0.21 & .minor_variant_rel_cov > 0.19
    )

    peak_frame_2 <- get_1d_peaks(.total_pair_cov, minor_freq_subsets[[1]], 3)
    peak_frame_3 <- get_1d_peaks(.total_pair_cov, minor_freq_subsets[[2]], 2)
    peak_frame_4 <- get_1d_peaks(.total_pair_cov, minor_freq_subsets[[3]], 1)
    peak_frame_5 <- get_1d_peaks(.total_pair_cov, minor_freq_subsets[[4]], 1)

    peak_frame_2$cov <- peak_frame_2$cov / (2 * round(peak_frame_2$cov / peak_frame_2$cov[1]))
    peak_frame_3$cov <- peak_frame_3$cov / (3 * round(peak_frame_3$cov / peak_frame_3$cov[1]))
    peak_frame_4$cov <- peak_frame_4$cov / (4 * round(peak_frame_4$cov / peak_frame_4$cov[1]))
    peak_frame_5$cov <- peak_frame_5$cov / (5 * round(peak_frame_5$cov / peak_frame_5$cov[1]))

    peak_frame <- rbind(peak_frame_2, peak_frame_3, peak_frame_4, peak_frame_5)
    weighted.mean(peak_frame$cov, peak_frame$height, na.rm = T)
}
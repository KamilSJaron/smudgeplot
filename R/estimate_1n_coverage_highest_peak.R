#' @title estimate_1n_coverage_highest_peak
#'
#' @description
#' \code{estimate_1n_coverage_highest_peak}
#'
#' @export

estimate_1n_coverage_highest_peak <- function(.filt_peak_sizes, .minor_variant_rel_cov, .total_pair_cov){
    # 3 indications
    # AAB is presumabely highest
    # AAB is has the lowest coverages that correspondents to mutliples of two
    # the is not other peak around 1/6 and 1/2
    highest_peak <- which.max(.filt_peak_sizes$rel_size)

    trinploid_1n_est <- get_trinoploid_1n_est(.filt_peak_sizes)
    if(is.na(trinploid_1n_est)){
        # assuming that the greatest peak is the diploid one
        draft_1n <- .filt_peak_sizes$pair_cov[highest_peak] / 2
    } else {
        draft_1n <- trinploid_1n_est
    }

    highest_min_cov <- .filt_peak_sizes$minor_variant_cov_rounded[highest_peak]
    subset <- .minor_variant_rel_cov < (highest_min_cov + 0.02) & .minor_variant_rel_cov > (highest_min_cov - 0.02)
    major_peaks <- get_1d_peaks(.total_pair_cov, subset, 5)
    major_cov <- major_peaks$cov[which.max(major_peaks$height)]
    major_cov / round(major_cov / draft_1n)
}
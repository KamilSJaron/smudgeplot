#' @title get_trinoploid_1n_est
#'
#' @description
#' \code{get_trinoploid_1n_est}
#'
#' @export

get_trinoploid_1n_est <- function(.filt_peak_sizes){
    trinploids <- .filt_peak_sizes$minor_variant_cov_rounded == 0.33
    if(any(trinploids)){
        triplod_1n_guess <- min(.filt_peak_sizes$pair_cov[trinploids]) / 3
        genome_counts_by_trip_1n_est <- round(.filt_peak_sizes$pair_cov / triplod_1n_guess)
        if(sum(genome_counts_by_trip_1n_est == 3) == 1){
            # min triploid peak look like triploid
            return(triplod_1n_guess)
        }
    }
    # unable to guess
    return(NA)
}

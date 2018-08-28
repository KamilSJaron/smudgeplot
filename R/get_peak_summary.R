#' @title get_peak_summary
#'
#' @description
#' \code{get_peak_summary}
#'
#' @export

get_peak_summary <- function(.peak_points, .k){
    peak_sizes <- get_peak_sizes(.peak_points)

    filt_peak_points <- filter_peaks(.peak_points, peak_sizes)
    filt_peak_sizes <- filter_peak_sizes(peak_sizes)

    filt_peak_sizes <- merge(filt_peak_sizes, filt_peak_points[filt_peak_points$summit == T,])
    filt_peak_sizes$minor_variant_cov <- transform_minor_variant_cov(filt_peak_sizes$y, .k)
    filt_peak_sizes$pair_cov <- transform_pair_cov(filt_peak_sizes$y, .k)

    # genome_copy_nuber <- round(filt_peak_sizes$pair_cov / (min(filt_peak_sizes$pair_cov) / 2))
    filt_peak_sizes$minor_variant_cov_rounded <- sapply(filt_peak_sizes$minor_variant_cov, round_minor_variant_cov)

    return(filt_peak_sizes)

}
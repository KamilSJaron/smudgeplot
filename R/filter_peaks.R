#' @title filter_peaks
#'
#' @description
#' \code{filter_peaks}
#'
#' @export

filter_peaks <- function(.peak_points, .peak_sizes, threshold = 0.005){
    to_filter <- .peak_points$peak %in% which(.peak_sizes$rel_size < threshold)

    .peak_points$peak[to_filter] <- NA
    .peak_points$summit[to_filter] <- F
    return(.peak_points)
}
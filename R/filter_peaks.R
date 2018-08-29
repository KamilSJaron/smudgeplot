#' @title filter_peaks
#'
#' @description
#' \code{filter_peaks}
#'
#' @export

filter_peaks <- function(.peak_points, .peak_sizes){
    to_filter <- !.peak_points$peak %in% .peak_sizes$peak

    .peak_points$peak[to_filter] <- NA
    .peak_points$summit[to_filter] <- F
    return(.peak_points)
}
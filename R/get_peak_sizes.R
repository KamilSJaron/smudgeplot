#' @title get_peak_sizes
#'
#' @description
#' \code{get_peak_sizes}
#'
#' @export

get_peak_sizes <- function(.peak_points){
    peaks <- .peak_points[.peak_points$summit == T,'peak']
    peak_sizes <- sapply(peaks, function(x){sum(.peak_points[.peak_points$peak == x, 'vals'], na.rm = T)})
    data.frame(peak = peaks,
                abs_size = peak_sizes,
                rel_size = peak_sizes / sum(peak_sizes))
}
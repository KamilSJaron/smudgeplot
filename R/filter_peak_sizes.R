#' @title filter_peak_sizes
#'
#' @description
#' \code{filter_peak_sizes}
#'
#' @export

filter_peak_sizes <- function(.peak_sizes, .threshold = 0.005){
    .peak_sizes[!peak_sizes$rel_size < .threshold,]
}
#' @title plot_seq_error_line
#'
#' @description
#' \code{plot_seq_error_line} add error line to plot it's L / cov
#'
#' @export

plot_seq_error_line <- function(.total_pair_covm, .L = NA){
    if (is.na(.L)){ .L <- min(.total_pair_cov) / 2 }
    max_cov_pair <- max(.total_pair_cov)
    cov_range <- seq(2 * .L, max_cov_pair, length = 30)
    lines(.L / cov_range, cov_range, lwd = 2, lty = 2)
}
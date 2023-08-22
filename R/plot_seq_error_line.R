#' @title plot_seq_error_line
#'
#' @description
#' \code{plot_seq_error_line} add error line to plot it's L / cov
#'
#' @export

plot_seq_error_line <- function(.cov_tab, .L = NA, .col = 'black'){
    if (is.na(.L)){ .L <- min(.cov_tab[, 'covB']) }
    max_cov_pair <- max(.cov_tab[, 'total_pair_cov'])
    cov_range <- seq(2 * .L, max_cov_pair, length = 500)
    lines((.L - 1) / cov_range, cov_range, lwd = 2.5, lty = 2, col = .col)
}
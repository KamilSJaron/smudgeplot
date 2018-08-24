#' @title plot_seq_error_line
#'
#' @description
#' \code{plot_seq_error_line} add error line
#' TODO it's misplaced a bit I think, fix needed
#'
#' @export

plot_seq_error_line <- function(total_pair_cov){
    min_kmer_cov <- min(total_pair_cov) / 2
    max_cov_pair <- max(total_pair_cov)
    cov_range <- seq(2 * min_kmer_cov,max_cov_pair, length = 30)
    lines(min_kmer_cov / cov_range, cov_range, lwd = 2, lty = 2)
}
#' @title guess_genome_structure
#'
#' @description
#' \code{guess_genome_structure}
#'
#' @export

guess_genome_structure <- function(.filt_peak_sizes_line, .n){
    both_count_in_genome <- round(as.numeric(.filt_peak_sizes_line['pair_cov']) / .n)
    minor_count_in_genome <- round(as.numeric(.filt_peak_sizes_line['minor_variant_cov']) * both_count_in_genome)
    paste0(
        c(rep('A', both_count_in_genome - minor_count_in_genome),
          rep('B', minor_count_in_genome)), collapse = '')
}
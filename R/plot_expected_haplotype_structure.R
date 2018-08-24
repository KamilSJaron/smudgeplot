#' @title plot_expected_haplotype_structure
#'
#' @description
#' \code{plot_expected_haplotype_structure} adds to a plot
#' TODO do automatically till y max? or max plody?
#' TODO or to annotate only significant peaks
#'
#' @export

plot_expected_haplotype_structure <- function(.n, .cex = 1.4){
    text(1/2 - 0.027, 2 * .n, 'AB', offset = 0, cex = .cex)
    text(1/3, 3 * .n, 'AAB', offset = 0, cex = .cex)
    text(1/4, 4 * .n, 'AAAB', offset = 0, cex = .cex)
    text(1/2 - 0.04, 4 * .n, 'AABB', offset = 0, cex = .cex)
    text(2/5, 5 * .n, 'AAABB', offset = 0, cex = .cex)
    text(1/5, 5 * .n, 'AAAAB', offset = 0, cex = .cex)
    text(3/6 - 0.055, 6 * .n, 'AAABBB', offset = 0, cex = .cex)
    text(2/6, 6 * .n, 'AAAABB', offset = 0, cex = .cex)
    text(1/6, 6 * .n, 'AAAAAB', offset = 0, cex = .cex)
}
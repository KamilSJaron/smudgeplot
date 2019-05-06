#' @title plot_expected_haplotype_structure
#'
#' @description
#' \code{plot_expected_haplotype_structure} adds genome scture labels to a smudgeplot
#'
#' @export

plot_expected_haplotype_structure <- function(.n, .peak_sizes,
                                              .adjust = F, .cex = 1.3, xmax = 0.49){
    borercases <- .peak_sizes$corrected_minor_variant_cov == 0.5
    for(i in 1:nrow(.peak_sizes)){
        # xmax is in the middle of the last square in the 2d histogram,
        # which is too far from the edge, so I average it with 0.49
        # witch will pull the label bit to the edge
        text( ifelse( borercases[i] & .adjust, (xmax + 0.49) / 2, .peak_sizes$corrected_minor_variant_cov[i]),
             .peak_sizes$ploidy[i] * .n, .peak_sizes$structure[i],
             offset = 0, cex = .cex, xpd = T, pos = ifelse( borercases[i] & .adjust, 2, 1))
    }
}
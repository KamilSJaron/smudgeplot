#' @title plot_expected_haplotype_structure
#'
#' @description
#' \code{plot_expected_haplotype_structure} adds to a plot
#'
#' @export

plot_expected_haplotype_structure <- function(.n, .peak_sizes,
                                              .adjust = F, .cex = 1.3){
    borercases <- .peak_sizes$minor_variant_cov_rounded == 0.5
    if(.adjust){
        # I adjust labels on 0.5 because otherwise they would be out of plot
        .peak_sizes$minor_variant_cov_rounded[borercases] <-
              .peak_sizes$minor_variant_cov_rounded[borercases] -
              ((0.007 * .peak_sizes$ploidy[borercases]) + 0.01)
    } else {
        # this is because xlim for even squares (yeah, I don't really like this function)
        .peak_sizes$minor_variant_cov_rounded[borercases] <- 0.48
    }

    for(i in 1:nrow(.peak_sizes)){
        text(.peak_sizes$minor_variant_cov_rounded[i],
             .peak_sizes$ploidy[i] * .n, .peak_sizes$structure[i],
             offset = 0, cex = .cex, xpd = T)
    }
}
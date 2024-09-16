#' @title plot_expected_haplotype_structure
#'
#' @description
#' \code{plot_expected_haplotype_structure} adds genome scture labels to a smudgeplot
#'
#' @export

plot_expected_haplotype_structure <- function(.n, .peak_sizes,
                                              .adjust = F, .cex = 1.3, xmax = 0.49){
    borercases <- .peak_sizes$corrected_minor_variant_cov == 0.5
    structures <- reduce_structure_representation(.peak_sizes)

    for(i in 1:nrow(.peak_sizes)){
        # xmax is in the middle of the last square in the 2d histogram,
        # which is too far from the edge, so I average it with 0.49
        # witch will pull the label bit to the edge
        text( ifelse( borercases[i] & .adjust, (xmax + 0.49) / 2, .peak_sizes$corrected_minor_variant_cov[i]),
             .peak_sizes$ploidy[i] * .n, structures[i],
             offset = 0, cex = .cex, xpd = T, pos = ifelse( borercases[i] & .adjust, 2, 1))
    }
}

#' @export
plot_all_smudge_labels <- function(cov_est, ymax, xmax = 0.49, .cex = 1.3, .L = 4, err = F){
    if (err){
        for (As in 1:(floor(ymax / cov_est) - 1)){
            label <- paste0(As, "Aerr")
            text(.L / (As * cov_est), (As * cov_est) + .L, label,
                    offset = 0, cex = .cex, xpd = T, pos = ifelse(As == 1, 3, 4))
        }
    }
    for (ploidy in 2:floor(ymax / cov_est)){
        for (Bs in 1:floor(ploidy / 2)){
            As = ploidy - Bs
            label <- paste0(As, "A", Bs, "B")
            text(ifelse(As == Bs, (xmax + 0.49)/2, Bs / ploidy), ploidy * cov_est, label, 
                 offset = 0, cex = .cex, xpd = T, 
                 pos = ifelse(As == Bs, 2, 1))
        }
    }
}
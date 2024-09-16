#' @title plot_alt
#'
#' @description
#' \code{plot_alt} is an alternative plotting introduced for 0.3.0 oriel. Instead of summing all the coverages falling within boundaries of squares, we make a tiling that will cover each individual point
#'
#' @export

plot_alt <- function(cov_tab, ylim, colour_ramp){
    A_equals_B <- cov_tab[, 'covA'] == cov_tab[, 'covB']
    cov_tab[A_equals_B, 'freq'] <- cov_tab[A_equals_B, 'freq'] * 2
    cov_tab$col <- colour_ramp[1 + round(31 * cov_tab[, 'freq'] / max(cov_tab[, 'freq']))]

    # c(bottom, left, top, right)
    par(mar=c(0,0,1,1))
    plot(NULL, xlim = c(0, 0.5), ylim = ylim,
         xlab = 'Normalized minor kmer coverage: B / (A + B)',
         ylab = 'Total coverage of the kmer pair: A + B', cex.lab = 1.4, bty = 'n')
    min_cov_to_plot <- max(ylim[1],min(cov_tab[, 'total_pair_cov']))
    nothing <- sapply(min_cov_to_plot:ylim[2], plot_one_coverage, cov_tab)
    return(0)
}

plot_one_coverage <- function(cov, cov_tab){
    cov_row_to_plot <- cov_tab[cov_tab[, 'total_pair_cov'] == cov, ]
    width <- 1 / (2 * cov)
    cov_row_to_plot$left <- cov_row_to_plot[, 'minor_variant_rel_cov'] - width
    cov_row_to_plot$right <- sapply(cov_row_to_plot[, 'minor_variant_rel_cov'], function(x){ min(0.5, x + width)})
    apply(cov_row_to_plot, 1, plot_one_box, cov)
}

plot_one_box <- function(one_box_row, cov){
    left <- as.numeric(one_box_row['left'])
    right <- as.numeric(one_box_row['right'])
    rect(left, cov - 0.5, right, cov + 0.5, col = one_box_row['col'], border = NA)
}

plot_isoA_line <- function (.covA, .L, .col = "black", .ymax = 250, .lwd, .lty) {
    min_covB <- .L # min(.cov_tab[, 'covB']) # should be L really
    max_covB <- .covA
    B_covs <- seq(min_covB, max_covB, length = 500)
    isoline_x <- B_covs/ (B_covs + .covA)
    isoline_y <- B_covs + .covA
    lines(isoline_x[isoline_y < .ymax], isoline_y[isoline_y < .ymax], lwd = .lwd, lty = .lty, col = .col)
}

plot_isoB_line <- function (.covB, .ymax, .col = "black", .lwd, .lty) {
    cov_range <- seq((2 * .covB) - 2, .ymax, length = 500)
    lines((.covB)/cov_range, cov_range, lwd = .lwd, lty = .lty, col = .col)
}

plot_iso_grid <- function(.cov, .L, .ymax, .col = 'black', .lwd = 2, .lty = 2){
    for (i in 0:15){
        cov <- (i + 0.5) * .cov
        plot_isoA_line(cov, .L = .L, .ymax = .ymax, .col, .lwd = .lwd, .lty = .lty)
        if (i < 8){
            plot_isoB_line(cov, .ymax, .col, .lwd = .lwd, .lty = .lty)
        }
    }
}
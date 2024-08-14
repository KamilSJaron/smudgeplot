plot_alt <- function(cov_tab, ylim, colour_ramp, logscale = F){
    A_equals_B <- cov_tab[, 'covA'] == cov_tab[, 'covB']
    cov_tab[A_equals_B, 'freq'] <- cov_tab[A_equals_B, 'freq'] * 2
    if (logscale){
        cov_tab[, 'freq'] <- log10(cov_tab[, 'freq'])
    }
    cov_tab$col <- colour_ramp[1 + round(31 * cov_tab[, 'freq'] / max(cov_tab[, 'freq']))]

    plot(NULL, xlim = c(0, 0.5), ylim = ylim,
         xlab = 'Normalized minor kmer coverage: B / (A + B)',
         ylab = 'Total coverage of the kmer pair: A + B', cex.lab = 1.4)
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

plot_dot_smudgeplot <- function(cov_tab, colour_ramp, xlim, ylim, background_col = 'grey', cex = 0.4){
    # this is the adjustment for plotting
    cov_tab[cov_tab$covA == cov_tab$covB, 'freq'] <- cov_tab[cov_tab$covA == cov_tab$covB, 'freq'] * 2
    cov_tab$col = colour_ramp[1 + round(31 * cov_tab$freq / max(cov_tab$freq))]

    plot(NULL, xlim = xlim, ylim = ylim, xlab = 'Normalized minor kmer coverage: B / (A + B)',
         ylab = 'Total coverage of the kmer pair: A + B')
    rect(xlim[1], ylim[1], xlim[2], ylim[2], col = background_col, border = NA)
    points(cov_tab[, 'minor_variant_rel_cov'], cov_tab[, 'total_pair_cov'], col = cov_tab$col, pch = 20, cex = cex)
}

plot_peakmap <- function(cov_tab, xlim, ylim, background_col = 'grey', cex = 2){
    # this is the adjustment for plotting
    plot(NULL, xlim = xlim, ylim = ylim, xlab = 'Normalized minor kmer coverage: B / (A + B)',
         ylab = 'Total coverage of the kmer pair: A + B')
    points(cov_tab[, 'minor_variant_rel_cov'], cov_tab[, 'total_pair_cov'], col = cov_tab$peak, pch = 20, cex = cex)
    legend('bottomleft', col = 1:8, pch = 20, title = 'peak', legend = 1:8)
}

plot_seq_error_line <- function (.cov_tab, .L = NA, .col = "black") {
    if (is.na(.L)) {
        .L <- min(.cov_tab[, "covB"])
    }
    max_cov_pair <- max(.cov_tab[, "total_pair_cov"])
    cov_range <- seq((2 * .L) - 2, max_cov_pair, length = 500)
    lines((.L - 1)/cov_range, cov_range, lwd = 2.5, lty = 2, 
        col = .col)
}
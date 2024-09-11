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
    points(cov_tab[, 'minor_variant_rel_cov'], cov_tab[, 'total_pair_cov'], col = cov_tab$smudge, pch = 20, cex = cex)
    legend('bottomleft', col = 1:8, pch = 20, title = 'smudge', legend = 1:8)
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

plot_smudge_labels <- function(cov_est, ymax, xmax = 0.49, .cex = 1.3, .L = 4){
    for (As in 1:(floor(ymax / cov_est) - 1)){
        label <- paste0(As, "Aerr")
        text(.L / (As * cov_est), (As * cov_est) + .L, label, 
                 offset = 0, cex = .cex, xpd = T, pos = ifelse(As == 1, 3, 4))
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

create_smudge_container <- function(cov, cov_tab, smudge_filtering_threshold){
    smudge_container <- list()
    total_genomic_kmers <- sum(cov_tab[ , 'freq'])
    
    for (Bs in 1:8){
        cov_tab_isoB <- cov_tab[cov_tab[ , 'covB'] > cov * ifelse(Bs == 1, 0, Bs - 0.5) & cov_tab[ , 'covB'] < cov * (Bs + 0.5), ]
        # cov_tab_isoB[, 'Bs'] <- Bs
        cov_tab_isoB[, 'As'] <- round(cov_tab_isoB[, 'covA'] / cov) #these are be individual smudges cutouts given coverages
        cov_tab_isoB[cov_tab_isoB[, 'As'] == 0, 'As'] = 1
        for( As in Bs:(16 - Bs)){
            cov_tab_one_smudge <- cov_tab_isoB[cov_tab_isoB[, 'As'] == As, ]
            if (sum(cov_tab_one_smudge[, 'freq']) / total_genomic_kmers > smudge_filtering_threshold){
                label <- paste0(As, "A", Bs, "B") 
                smudge_container[[label]] <- cov_tab_one_smudge[,-which(names(cov_tab_one_smudge) %in% c('is_error', 'As'))]
            }
        }
    }
    return(smudge_container)
}
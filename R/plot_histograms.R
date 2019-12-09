#' @title plot_histograms
#'
#' @description
#' \code{plot_histograms} makes 2 plots - histograms of the both margins of smudge plot
#'
#' @export

plot_histograms <- function(.minor_variant_rel_cov, .total_pair_cov, .ymax, .smudge_summary,
                            .nbins, .fig_title = NA, .cex = 1.4, .col = NA){
    to_filter <- .total_pair_cov < .ymax - (.ymax / .nbins)
    .total_pair_cov <- .total_pair_cov[to_filter]
    .minor_variant_rel_cov <- .minor_variant_rel_cov[to_filter]
    h1 <- hist(.minor_variant_rel_cov, breaks = 100, plot = F)
    h2 <- hist(.total_pair_cov, breaks = 100, plot = F)
    top <- max(h1$counts, h2$counts)

    if( is.na(.col) ){
        .col <- rgb(0.8352, 0.2431, 0.3098)
    }

    # minor_variant_rel_cov HISTOGRAM - top
    par(mar=c(0,3.8,1,0))
    barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col = .col)
    if(!(is.na(.fig_title))){
        mtext(bquote(italic(.(.fig_title))), side=3, adj=0, line=-3, cex = .cex + 0.2)
    }

    ploidytext <- switch(.smudge_summary$genome_ploidy - 1,
                           p2 = 'diploid',
                           p3 = 'triploid',
                           p4 = 'tetraploid',
                           p5 = 'pentaploid',
                           p6 = 'hexaploid',
                           p7 = 'heptaploid',
                           p8 = 'octoploid')

    if(!(is.na(.smudge_summary$genome_ploidy))){
        mtext(paste('proposed', ploidytext), side=3, adj=0.05, line=-5, cex = .cex - 0.2)
    }

    # total pair coverage HISTOGRAM - right
    par(mar=c(3.8,0,0.5,1))
    barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col = .col, horiz = T)
    legend('bottomright', bty = 'n', paste('1n = ', round(.smudge_summary$n)), cex = .cex - 0.1)

    .peak_sizes <- .smudge_summary$peak_sizes[,c(11,3)]
    colnames(.peak_sizes) <- c('peak', 'size')
    to_remove <- sapply(.peak_sizes[,1], nchar) > 6
    .peak_sizes <- .peak_sizes[!to_remove, ]
    if( any(to_remove) ){
        .peak_sizes <- rbind(.peak_sizes, data.frame(peak = 'others', 'size' = 1 - sum(.peak_sizes[,2])) )
    }
    if(! any(is.na(.peak_sizes))){
        legend('topleft', bty = 'n', .peak_sizes[,1], cex = .cex - 0.2)
        legend('topright', bty = 'n', legend = round(.peak_sizes[,2], 2), cex = .cex - 0.2)
    }
}

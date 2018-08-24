#' @title plot_histograms
#'
#' @description
#' \code{plot_histograms} makes 2 plots - histograms of the both margins of smudge plot
#'
#' @export

plot_histograms <- function(.minor_variant_rel_cov, .total_pair_cov,
                            .fig_title = NA, .cex = 1.4){
    h1 <- hist(.minor_variant_rel_cov, breaks = 100, plot = F)
    h2 <- hist(.total_pair_cov, breaks = 100, plot = F)
    top <- max(h1$counts, h2$counts)

    # minor_variant_rel_cov HISTOGRAM - top
    par(mar=c(0,3.8,1,0))
    barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col = pal[2])
    if(!(is.na(.fig_title))){
        mtext(bquote(italic(.(.fig_title))), side=3, adj=0, line=-3, cex = .cex + 0.2)
    }

    # total pair coverage HISTOGRAM - right
    par(mar=c(3.8,0,0.5,1))
    barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col = pal[2], horiz = T)
    mtext(paste('1n = ', n), side=1, adj=0.8, line=-2, cex = .cex)
}
#' @title plot_histograms
#'
#' @description
#' \code{plot_histograms} makes 2 plots - histograms of the both margins of smudge plot
#'
#' @export

# .cov_tab <- cov_tab
# .smudge_summary <- smudge_summary
# .fig_title <- NA
# .ylim <- ylim
# .bins <- 150

plot_histograms <- function(.cov_tab, .smudge_summary,
                            .fig_title = NA, .cex = 1.4, .col = NA,
                            .ylim = NA, .bins = 100){

    if( is.na(.col) ){
        .col <- rgb(0.8352, 0.2431, 0.3098)
    }

    # removing pairs with excess coverage
    .cov_tab <- .cov_tab[.cov_tab[, 'total_pair_cov'] < .ylim[2], ]
   
    # minor_variant_rel_cov HISTOGRAM - top
    coverage_histogram(.cov_tab, bins = .bins, xlim = c(0, 0.5), 
                       column = 'minor_variant_rel_cov', horiz = FALSE, .col)

    if (!(is.na(.fig_title))) {
        mtext(bquote(italic(.(.fig_title))), side=3, adj=0, line=-3, cex = .cex + 0.2)
    }

    if (length(.smudge_summary$genome_ploidy) > 0){
        if ( .smudge_summary$genome_ploidy < 9){
            ploidytext <- switch(.smudge_summary$genome_ploidy - 1,
                                p2 = 'diploid',
                                p3 = 'triploid',
                                p4 = 'tetraploid',
                                p5 = 'pentaploid',
                                p6 = 'hexaploid',
                                p7 = 'heptaploid',
                                p8 = 'octoploid')
        } else {
            ploidytext <- paste0(.smudge_summary$genome_ploidy, '-ploid')
        }

        if(!(is.na(.smudge_summary$genome_ploidy))){
            mtext(paste('proposed', ploidytext), side=3, adj=0.05, line=-5, cex = .cex - 0.2)
        }
    } 


    # total pair coverage HISTOGRAM - right
    coverage_histogram(.cov_tab, bins = .bins, xlim = .ylim, # this is because horizonal is TRUE 
                       column = 'total_pair_cov', horiz = TRUE, .col)

    if (.smudge_summary$n > 0){ #0 means it was not infered
        legend('bottomright', bty = 'n', paste('1n = ', round(.smudge_summary$n)), cex = .cex - 0.1)

        .peak_sizes <- data.frame(peak = reduce_structure_representation(.smudge_summary$peak_sizes),
                                size = .smudge_summary$peak_sizes[,3])

        if(! any(is.na(.peak_sizes))){
            legend('topleft', bty = 'n', .peak_sizes[,1], cex = .cex - 0.2)
            legend('topright', bty = 'n', legend = round(.peak_sizes[,2], 2), cex = .cex - 0.2)
        }
    }                       
}

#' @title plot_smudgeplot
#'
#' @description
#' \code{plot_smudgeplot} is the core smudgeplot function to plot the core of the smudgeplot
#' the cental 2D histogram/ landscape plot / heatmap / howeveryouwantyoucallthis plot
#'
#' @param .k - the matrix of densities or list of x,y axies and z, the matrix. Usually output of function
#'
#' @param .n - the 1n coverage, the coverage of the haploid genome (can be estimated by TODO function)
#'
#' @param .colour_ramp - the colour pallete for the heat colours
#'
#' @param .cex - the size of the axis labels [default 1.4]
#'
#' @author Kamil Jaron \email{kamiljaron at gmail.com}
#'
#' @export

plot_smudgeplot <- function(.k, .n, .colour_ramp, .cex = 1.4){
    # margins 'c(bottom, left, top, right)'
    par(mar=c(4.8,4.8,1,1))
    # 2D HISTOGRAM
    image(.k, col = .colour_ramp,
        xlab = 'Normalized minor kmer coverage: B / (A + B)',
        ylab = 'Total coverage of the kmer pair: A + B', cex.lab = 1.4,
        axes=F
    )

    if (.n == 0){
        axis(2)
    } else {
        axis(2, at = 2:8 * .n, labels = paste(2:8, 'n'))
    }
    axis(1, at = c(1/5, 1/4, 1/3, max(.k$x)),
            labels = c('1/5', '1/4', '1/3', '1/2'))
}
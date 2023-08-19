#' @title coverage_histogram
#'
#' @description
#' \code{coverage_histogram} plots a histogram from a table that contains both values (what is normally x) and their respective frequencies (a column called 'freq');
#'
#' @export

# coverage_histogram(.cov_tab, bins = .bins, xlim = c(0, 0.5), 
#                     column = 'minor_variant_rel_cov', horiz = FALSE, .col)
# bins <- .bins
# xlim <- c(0, 0.5)
# column <- 'minor_variant_rel_cov'
# horiz <- FALSE

coverage_histogram <- function(.cov_tab, bins, xlim, column, horiz, .col){

    # define breaks of the histogram
    prettified_bins <- length(pretty(c(xlim[1],xlim[2]), n = bins))
    x_breaks <- seq(0, xlim[2], length = prettified_bins)
    # identify buckets of each kmer pair
    # Note: these are right open intervals so the last ratio bucket contain 0.5
    kmerpairs2buckets <- list(findInterval(.cov_tab[, column], x_breaks, left.open = TRUE))
    # and for each bucket aggregate their frequencies
    agregated_freq <- aggregate(.cov_tab[, 'freq'], by = kmerpairs2buckets, sum)

    # from, to are the borders of the bucket, pos is for plotting, and freq is the agregated value 
    ratio_hist <- data.frame(from = x_breaks[1:(prettified_bins-1)], to = x_breaks[2:prettified_bins], freq = 0)
    ratio_hist[, 'pos'] <- (ratio_hist[, 'from'] + ratio_hist[, 'to']) / 2
    ratio_hist[agregated_freq[, 1], 'freq'] <- agregated_freq[, 2] 

    width = (ratio_hist[1, 'to'] - ratio_hist[1, 'from']) / 2

    # horizontal plot
    if (horiz){
        # c(bottom, left, top, right)
        par(mar=c(3.8,0,0,1))
        plot(NULL, type="n", xlab="", ylab="", axes=FALSE,
             ylim=xlim, xlim=c(0, max(ratio_hist[, 'freq'])))
        for ( i in 1:prettified_bins) {
            rect(0, ratio_hist[i, 'pos'] - width, ratio_hist[i, 'freq'], ratio_hist[i, 'pos'] + width, col = .col, border = FALSE)
        }
    } else {
        # vertical plot
        par(mar=c(0,3.8,1,0))
        plot(NULL, type="n", xlab="", ylab="", axes=FALSE,
             ylim=c(0, max(ratio_hist[, 'freq'])), xlim=xlim)
        for ( i in 1:prettified_bins) {
            rect(ratio_hist[i, 'pos'] - width, 0, ratio_hist[i, 'pos'] + width, ratio_hist[i, 'freq'], col = .col, border = FALSE)
        }
    }
}


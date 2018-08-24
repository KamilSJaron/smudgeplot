#' @title plot_legend
#'
#' @description
#' \code{plot_legend} generate topright corener in the smudgeplot
#'
#' @export

plot_legend <- function(.k, .total_pair_cov, .colour_ramp, .sqrt_scale = T){
    par(mar=c(0,0,2,1))
    plot.new()
    title('kmers pairs')
    for(i in 1:32){
        rect(0,(i - 0.01) / 33, 0.5, (i + 0.99) / 33, col = .colour_ramp[i])
    }
    if( .sqrt_scale == T ){
        # TODO correct this scale
        kmer_max <- (length(.total_pair_cov) * max((.k$z)^2)) / sum((.k$z)^2)
        for(i in 0:6){
            text(0.75, i / 6, round((sqrt(kmer_max) * i)^2 / 6000) * 1000, offset = 0)
        }
    } else {
        kmer_max <- (length(.total_pair_cov) * max(.k$z)) / sum(.k$z)
        for(i in 0:6){
            text(0.75, i / 6, round((kmer_max * i) / 6000) * 1000, offset = 0)
        }
    }
}
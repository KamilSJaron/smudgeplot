#' @title plot_legend
#'
#' @description
#' \code{plot_legend} generate topright corener in the smudgeplot
#'
#' @export

plot_legend <- function(.k, .colour_ramp, .log_scale = T){
    par(mar=c(0,0,2,1))
    plot.new()
    print_title <- ifelse(.log_scale, 'log kmers pairs', 'kmers pairs')
    title(print_title)
    for(i in 1:32){
        rect(0,(i - 0.01) / 33, 0.5, (i + 0.99) / 33, col = .colour_ramp[i])
    }
    kmer_max <- max(smudge_container$dens)
    if( .log_scale == T ){
        for(i in 0:6){
            text(0.75, i / 6, rounding(10^(log10(kmer_max) * i / 6)), offset = 0)
        }
    } else {
        for(i in 0:6){
            text(0.75, i / 6, rounding(kmer_max * i / 6), offset = 0)
        }
    }
}
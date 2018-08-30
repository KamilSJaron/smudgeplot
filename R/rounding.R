#' @title rounding
#'
#' @description
#' \code{rounding} rounds to thousands or hundreds if the number is smaller than 1000
#'
#' @export

rounding <- function(number){
    if(number > 1000){
        round(number / 1000) * 1000
    } else {
        round(number / 100) * 100
    }
}

plot_legend <- function(.k, .colour_ramp, .sqrt_scale = T){
    par(mar=c(0,0,2,1))
    plot.new()
    title('kmers pairs')
    for(i in 1:32){
        rect(0,(i - 0.01) / 33, 0.5, (i + 0.99) / 33, col = .colour_ramp[i])
    }
    kmer_max <- max(smudge_container$dens)
    if( .sqrt_scale == T ){
        for(i in 0:6){
            text(0.75, i / 6, rounding((sqrt(kmer_max) * i / 6)^2), offset = 0)
        }
    } else {
        for(i in 0:6){
            text(0.75, i / 6, rounding(kmer_max * i / 6), offset = 0)
        }
    }
}
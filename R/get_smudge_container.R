#' @title get_smudge_container
#'
#' @description
#' \code{get_smudge_container} is inspired by https://stackoverflow.com/a/18103689/2962344
#'
#' @export

get_smudge_container <- function(.minor_variant_rel_cov, .total_pair_cov,
                                 .nbins = 20, .xlim = c(0, 0.5), .ylim = NA){
    if( any(is.na(.ylim)) ){
        .ylim = c(0,max(.total_pair_cov))
    }
    smudge_container <- list()
    smudge_container$x <- seq(.xlim[1], ((.nbins - 1) / .nbins) * .xlim[2], length = .nbins)
    smudge_container$y <- c(seq(.ylim[1], ((.nbins - 1) / .nbins) * .ylim[2], length = .nbins), .ylim[2])

    freq <-  as.data.frame(table(findInterval(.minor_variant_rel_cov, smudge_container$x, left.open = T),
                                 findInterval(.total_pair_cov, smudge_container$y, left.open = T)),
                           stringsAsFactors = F)
    freq[,1] <- as.numeric(freq[,1])
    freq[,2] <- as.numeric(freq[,2])
    freq <- freq[!freq$Var2 == .nbins+1,]

    smudge_container$dens <- matrix(0, .nbins, .nbins)
    smudge_container$dens[cbind(freq[,1], freq[,2])] <- freq[,3]

    smudge_container$y <- smudge_container$y[1:.nbins]
    smudge_container$z <- sqrt(smudge_container$dens)

    return(smudge_container)
}
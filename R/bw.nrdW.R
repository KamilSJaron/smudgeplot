#' @title bw.nrdW
#'
#' @description
#' \code{bw.nrdW} is a function that mimics bw.nrd0 - Silverman's ‘rule of thumb’ for estimating the kernel width, but using weighted values;
#' this can be useful if one wants to mimics kerned smoothing to a histogram using counts as weights
#'
#' @export

bw.nrdW <- function (x, w)
{
    if (length(x) < 2L)
        stop("need at least 2 data points")
    Ex <- weighted.mean(x, w)
    hi <- sqrt(1 / (sum(w) - 1) * sum(w * (x - Ex)^2)) #wieighted SD
    if (!(lo <- min(hi, wtd.iqr(x, w)/1.34)))
        (lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1)
    0.9 * lo * sum(w)^(-0.2)
}

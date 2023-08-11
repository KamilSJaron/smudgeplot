#' @title wtd.quantile & wtd.iqr
#'
#' @description
#' \code{wtd.quantile}
#' calculates a quantile using "wighted data", specifically designed for "histogram-like" data
#' \code{wtd.iqr}
#' 0.75 - 0.25 weighted quantile
#'
#' @export

wtd.quantile <- function(x, q=0.25, weight=NULL) {
  o <- order(x)
  n <- sum(weight)
  order <- 1 + (n - 1) * q
  low  <- pmax(floor(order), 1)
  high <- pmin(ceiling(order), n)
  low_contribution <- high - order
  allq <- approx(x=cumsum(weight[o])/sum(weight), y=x[o], xout = c(low, high)/n, method = "constant",
      f = 1, rule = 2)$y
  low_contribution * allq[1] + (1 - low_contribution) * allq[2]
}

wtd.iqr <- function(x, w=NULL) {
  wtd.quantile(x, q=0.75, weight=w) - wtd.quantile(x, q=0.25, weight=w)
}

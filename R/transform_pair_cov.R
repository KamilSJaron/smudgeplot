#' @title transform_pair_cov
#'
#' @description
#' \code{transform_pair_cov}
#' to transform 1 ... 30 coordinates to y coordinates of the 2d hisrtogram
#'
#' @export

transform_pair_cov <- function(.y, .k){
    .k$y[.y]
}
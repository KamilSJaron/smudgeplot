#' @title transform_pair_cov
#'
#' @description
#' \code{transform_pair_cov}
#' to transform coordinates of smudgeplot to total_pair_cov
#'
#' @export

transform_pair_cov <- function(.y, .smudge_container){
    .smudge_container$y[.y] + (diff(.smudge_container$y)[1] / 2)
}
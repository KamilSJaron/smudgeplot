#' @title transform_minor_variant_cov
#'
#' @description
#' \code{transform_minor_variant_cov}
#' to transform coordinates of smudgeplot to minor_variant_cov
#'
#' @export

transform_minor_variant_cov <- function(.x, .k){
    .k$x[.x]
}
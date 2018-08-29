#' @title transform_minor_variant_cov
#'
#' @description
#' \code{transform_minor_variant_cov}
#' to transform coordinates of smudgeplot to minor_variant_rel_cov
#'
#' @export

transform_minor_variant_cov <- function(.x, .smudge_container){
    # .x corresponds to <.x, .x+1> interval
    # I add 1/2 of the size of the interval to return mean interval values
    .smudge_container$x[.x] + (diff(.smudge_container$x)[1] / 2)
}
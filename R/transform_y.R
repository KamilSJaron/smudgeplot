#' @title transform_y
#'
#' @description
#' \code{transform_y}
#' to transform 1 ... 30 coordinates to y coordinates of the 2d hisrtogram
#'
#' @export

transform_y <- function(.y, .min_kmerpair_cov, .max_cov_pair){
    (((.y - 1) / 29) * (.max_cov_pair - .min_kmerpair_cov)) + .min_kmerpair_cov
}
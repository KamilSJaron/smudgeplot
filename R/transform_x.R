#' @title transform_x
#'
#' @description
#' \code{transform_x}
#' to transform 1 ... 30 coordinates to x coordinates of the 2d hisrtogram
#'
#' @export

transform_x <- function(x, for_plot = T){
    orig <- (((x - 1) / 29) / 2)
    if ( for_plot ){
        ((orig / 0.5) * 0.46) + 0.02
    } else {
        orig
    }
}
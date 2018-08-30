#' @title smuge_summary
#'
#' @description
#' \code{smuge_summary} a function that is used to report all the summary output messages
#' like this I will be able to easily redirect the waring message to file as well
#'
#' @export

smuge_summary <- function(...){
    write(paste(..., sep = '\t'), stdout())
}

#' @title smuge_warn
#'
#' @description
#' \code{smuge_warn} a function that is used to report all the wanrning messages
#' like this I will be able to easily redirect the waring message to file as well
#'
#' @export

smuge_warn <- function(...){
    write(paste(..., sep = '\t'), stderr())
}

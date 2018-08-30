#' @title smudge_warn
#'
#' @description
#' \code{smudge_warn} a function that is used to report all the wanrning messages
#' like this I will be able to easily redirect the waring message to file as well
#'
#' @export

smudge_warn <- function(.filename, ...){
    write(paste(..., sep = '\t'), stderr())
    write(paste(..., sep = '\t'), paste0(.filename, "_warnings.txt"), append = T)
}
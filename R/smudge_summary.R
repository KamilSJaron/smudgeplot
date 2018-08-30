#' @title smudge_summary
#'
#' @description
#' \code{smudge_summary} a function that is used to report all the summary output messages
#' like this I will be able to easily redirect the waring message to file as well
#'
#' @export

smudge_summary <- function(.filename, ...){
    write(paste(..., sep = '\t'), paste0(.filename, "_verbose_summary.txt"), append = T)
}

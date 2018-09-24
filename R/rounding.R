#' @title rounding
#'
#' @description
#' \code{rounding} rounds to thousands or hundreds if the number is smaller than 1000
#'
#' @export

rounding <- function(number){
    if(number > 1000){
        round(number / 1000) * 1000
    } else if (number > 100){
        round(number / 100) * 100
    } else {
        round(number / 10) * 10
    }
}
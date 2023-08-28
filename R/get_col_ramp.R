#' @title get_col_ramp
#'
#' @description
#' \code{get_col_ramp} returns a colour ramp used by smudgeplots; 
#' args need to contain args$col_ramp and args$invert_cols features
#' the first needs to be a ramp-generating function (e.g. "viridis", "grey.colors", "magma" or "mako")
#' args$invert_cols is just TRUE or FALSE (if to revert the colours in the palete)
#'
#' @export

get_col_ramp <- function(.args, delay = 0){
    colour_ramp <- eval(parse(text = paste0(.args$col_ramp,"(", 32 - delay, ")")))
    if (.args$invert_cols){
        colour_ramp <- rev(colour_ramp)
    }
    colour_ramp <- c(rep(colour_ramp[1], delay), colour_ramp)
    return(colour_ramp)
}


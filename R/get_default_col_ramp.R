#' @title get_default_col_ramp
#'
#' @description
#' \code{get_default_col_ramp} returns default colour ramp used by smudgeplots
#'
#' @export

get_default_col_ramp <- function(){
    pal <- brewer.pal(11,'Spectral')
    rf <- colorRampPalette(rev(pal[3:11]))
    rf(32)
}


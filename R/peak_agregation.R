#' @title peak_agregation
#'
#' @description
#' \code{peak_agregation} is the function that agregates tiles from the smudge_container
#' into smudges - individual distributions
#'
#' @param .smudge_container - a smudge_container object; a list with `x`,`y` axies of the 2d histogram and `dens` matrix of 2h histogram values.
#'        Usually output of function \code{get_smudge_container}
#'
#' @export

peak_agregation <- function(.smudge_container, .L = NA){
    if ( is.na(.L) ){ .L <- min(.smudge_container$y) / 2 } # guess L if not specified
    .L <- max(c(.L, 10))
    nbins <- length(.smudge_container$x)
    peak_points <- data.frame(vals = as.vector(.smudge_container$dens),
                              x = rep(1:nbins, nbins),
                              y = rep(1:nbins, each = nbins))
    peak_points$peak <- NA
    peak_points$summit <- NA
    peak_points <- peak_points[order(peak_points$vals, decreasing = T),]
    #TODO no go zone x,y <- .smudge_container
    nogo_x <- findInterval(.L / .smudge_container$y, .smudge_container$x, left.open = T)
    # mark the first
    peak_points$peak[1] <- 1
    peak_points$summit[1] <- T

    for (point in 2:nrow(peak_points)) {
        parsed_points <- peak_points[1:point-1,]
        x_n <- abs(parsed_points[,'x'] - peak_points[point,'x']) <= 1
        y_n <- abs(parsed_points[,'y'] - peak_points[point,'y']) <= 1
        if (any(x_n & y_n)) {
            tiles_around <- parsed_points[x_n & y_n,]
            peak_points$peak[point] <- tiles_around$peak[which.max(tiles_around$vals)]
            peak_points$summit[point] <- F
        } else {
            # filt_peak_sizes$y -> pair_coverage
            # parse only point above eror line
            if ( peak_points[point,'x'] > nogo_x[peak_points[point,'y']] ) {
                peak_points$peak[point] <- max(parsed_points$peak, na.rm = T) + 1
                peak_points$summit[point] <- T
            } else {
                peak_points$summit[point] <- F
            }
        }
    }
    return(peak_points)
}

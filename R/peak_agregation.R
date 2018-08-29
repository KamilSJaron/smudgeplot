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

peak_agregation <- function(.smudge_container){
    nbins <- length(.smudge_container$x)
    peak_points <- data.frame(vals = as.vector(.smudge_container$dens),
                              x = rep(1:nbins, nbins),
                              y = rep(1:nbins, each = nbins))
    peak_points$peak <- NA
    peak_points$summit <- NA
    peak_points <- peak_points[order(peak_points$vals, decreasing = T),]
    # mark the first
    peak_points$peak[1] <- 1
    peak_points$summit[1] <- T

    for(point in 2:nrow(peak_points)){
        parsed_points <- peak_points[1:point-1,]
        x_n <- abs(parsed_points[,'x'] - peak_points[point,'x']) <= 1
        y_n <- abs(parsed_points[,'y'] - peak_points[point,'y']) <= 1
        if(any(x_n & y_n)){
            tiles_around <- parsed_points[x_n & y_n,]
            peak_points$peak[point] <- tiles_around$peak[which.max(tiles_around$vals)]
            peak_points$summit[point] <- F
        } else {
            peak_points$peak[point] <- max(parsed_points$peak) + 1
            peak_points$summit[point] <- T
        }
    }
    return(peak_points)
}

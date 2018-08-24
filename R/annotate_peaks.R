#' @title annotate_peaks
#'
#' @description
#' \code{annotate_peaks}
#'
#' @export

# annotate peaks
annotate_peaks <- function(peak_points, min_kmerpair_cov, max_kmerpair_cov){
    for (summit in which(peak_points$summit == T)){
        peak <- peak_points$peak[summit]
        to_plot <- peak_points[peak_points$peak == peak, c('x','y')]
        to_plot <- to_plot[!is.na(to_plot$x),]
        points(transform_x(to_plot$y), transform_y(to_plot$x, min_kmerpair_cov, max_kmerpair_cov), pch = peak)
    }
}
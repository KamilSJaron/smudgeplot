#' @title annotate_summits
#'
#' @description
#' \code{annotate_summits}
#'
#' @export

annotate_summits <- function(peak_points, peak_sizes, min_kmerpair_cov, max_kmerpair_cov, col = 'red'){
    for (summit in which(peak_points$summit == T)){
         peak <- peak_points$peak[summit]
         peak_size <- peak_sizes$rel_size[peak_sizes$peak == peak]
         # points(transform_x(peak_points$y[summit]),
         #        transform_y(peak_points$x[summit], min_kmerpair_cov, max_kmerpair_cov))
         text(transform_x(peak_points$y[summit]),
              transform_y(peak_points$x[summit], min_kmerpair_cov, max_kmerpair_cov),
              round(peak_size, 3),
              offset = 0, cex = 1.2, col = col, xpd=NA)
    }
}
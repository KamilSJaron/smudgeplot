#' @title get_1d_peaks
#'
#' @description
#' \code{get_1d_peaks} is a legacy code
#'
#' @export

get_1d_peaks <- function(.total_pair_cov, .subset, .num_of_peaks = 3, .adjust = 10, .peak_frame = data.frame(), .depth = 1){
    if (sum(.subset) < 2){
        return(data.frame(cov = NA, height = NA))
    }
    if (nrow(.peak_frame) < .num_of_peaks & .depth < 11){
        d <- density(.total_pair_cov[.subset], adjust = .adjust) # returns the density data
        # plot(d)
        selected_points <- which(diff(sign(diff(d$y))) == -2)
        .peak_frame <- data.frame(cov = d$x[selected_points], height = d$y[selected_points] * sum(.subset))
        .peak_frame <- get_1d_peaks(.total_pair_cov, .subset, .num_of_peaks, .adjust - 1, .peak_frame, .depth + 1)
    }
    .peak_frame
}
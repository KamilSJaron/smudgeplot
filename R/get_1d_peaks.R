#' @title get_1d_peaks
#'
#' @description
#' \code{get_1d_peaks}
#'
#' @export

get_1d_peaks <- function(subset, num_of_peaks = 3, adjust = 10, peak_frame = data.frame(), depth = 1){
    if (nrow(peak_frame) < num_of_peaks & depth < 11){
        d <- density(total_pair_cov[subset], adjust = adjust) # returns the density data
        # plot(d)
        selected_points <- which(diff(sign(diff(d$y))) == -2)
        peak_frame <- data.frame(cov = d$x[selected_points], height = d$y[selected_points] * sum(subset))
        peak_frame <- get_1d_peaks(subset, num_of_peaks, adjust - 1, peak_frame, depth + 1)
    }
    peak_frame
}
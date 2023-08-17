#' @title get_1d_peaks
#'
#' @description
#' \code{get_1d_peaks} is a legacy code
#'
#' @export

get_1d_peaks <- function(.cov_tab, .num_of_peaks = 3, .adjust = 10, .peak_frame = data.frame(), .depth = 1){
    if (nrow(.cov_tab) < 2){
        return(data.frame(cov = NA, height = NA))
    }
    if (nrow(.peak_frame) < .num_of_peaks & .depth < 11){
        kmer_pairs <- sum(.cov_tab[, 'freq'])
        # used weighted version of silverman's rule of thumb
        bw <- bw.nrdW(.cov_tab[, 'total_pair_cov'], .cov_tab[, 'freq'])
        # density wants to have the weights relative, I will divide it and them multiply the fitted heights (alternatively, I could suppres warnings)
        rel_weights <- .cov_tab[, 'freq'] / kmer_pairs
        d <- density(.cov_tab[, 'total_pair_cov'], weights = rel_weights, adjust = .adjust, bw = bw) # returns the density data
        selected_points <- which(diff(sign(diff(d$y))) == -2) + 1
        .peak_frame <- data.frame(cov = d$x[selected_points], height = d$y[selected_points] * kmer_pairs)
        .peak_frame <- get_1d_peaks(.cov_tab, .num_of_peaks, .adjust - 1, .peak_frame, .depth + 1)
    }
    .peak_frame
}
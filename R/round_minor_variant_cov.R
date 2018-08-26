#' @title round_minor_variant_cov
#'
#' @description
#' \code{round_minor_variant_cov}
#'
#' @export

round_minor_variant_cov <- function(.min_cov){
    # 0.5   0.4   0.33    0.25   0.2
    # 1/2   2/5   1/3     1/4    1/5
    expected_freq <- c(0.5, 0.4, 0.33, 0.25, 0.2)
    expected_freq[which.min(abs(expected_freq - .min_cov))]
}
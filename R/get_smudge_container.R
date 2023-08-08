#' @title get_smudge_container
#'
#' @description
#' \code{get_smudge_container} is inspired by https://stackoverflow.com/a/18103689/2962344
#'
#' @export

get_smudge_container <- function(.cov_tab, .nbins = 20, 
                                 .xlim = c(0, 0.5), .ylim = NA){

  total_pair_cov <- .cov_tab[, 1] + .cov_tab[, 2]
  minor_variant_rel_cov <- .cov_tab[, 1] / total_pair_cov
  
  if (any(is.na(.ylim))) {
    .ylim <- c(0, max(.total_pair_cov))
  }

  smudge_container <- list()
  smudge_container$x <- seq(.xlim[1], ((.nbins - 1) / .nbins) * .xlim[2], length = .nbins)
  smudge_container$y <- c(seq(.ylim[1], ((.nbins - 1) / .nbins) * .ylim[2], length = .nbins), .ylim[2])

  xcoords <- findInterval(minor_variant_rel_cov, smudge_container$x, left.open = TRUE)
  ycoords <- findInterval(total_pair_cov, smudge_container$y, left.open = TRUE) # all >.ylim[2] will be in nbin+1 th bin

  if (ncol(.cov_tab) == 3) { # ploidy-plot style "cov1 cov2 density" table
    group_by <- list(paste(xcoords, ycoords))
    agregated_counts <- aggregate(.cov_tab[, 3], by = group_by, sum)
    agregated_coordinates <- matrix(as.numeric(unlist(strsplit(agregated_counts[, 1], " "))), ncol = 2, byrow = TRUE)
    coordinates_in_range <- !agregated_coordinates[, 2] == .nbins+1 # the rows that belong to nbins+1 bin need to be removed; this is identify which to keep
    agregated_counts <- agregated_counts[coordinates_in_range, 2] # remove nbin+1 in counts, and get rid keep only the counts (not cooredinates)
    agregated_coordinates <- agregated_coordinates[coordinates_in_range, ] # and subset the separated coordinates too
  } else { # original smudgepot style "cov1 cov2" table;
    agregated_counts <-  as.data.frame(table(xcoords, ycoords), stringsAsFactors = FALSE)
    agregated_counts[, 1] <- as.numeric(agregated_counts[, 1])
    agregated_counts[, 2] <- as.numeric(agregated_counts[, 2])
    agregated_counts <- agregated_counts[!agregated_counts[, 2] == .nbins + 1, ]
    agregated_coordinates <- cbind(agregated_counts[, 1], agregated_counts[, 2])
    agregated_counts <- agregated_counts[, 3]
  }
  
  smudge_container$dens <- matrix(0, .nbins, .nbins)
  smudge_container$dens[agregated_coordinates] <- agregated_counts

  smudge_container$y <- smudge_container$y[1:nbins]
  smudge_container$z <- log10(smudge_container$dens)
  smudge_container$z[is.infinite(smudge_container$z)] <- 0

  return(smudge_container)
}
####### SAMPLE

library(smudgeplot)

cov_tab_sample <- read.table('data/testing_Scer.tsv', col.names = c('covB', 'covA', 'freq'))
cov_tab_sample[, 'total_pair_cov'] <- cov_tab_sample[, 1] + cov_tab_sample[, 2]
cov_tab_sample[, 'minor_variant_rel_cov'] <- cov_tab_sample[, 1] / cov_tab_sample[, 'total_pair_cov']

xlim <- c(0, 0.5)
ylim <- c(0, 150)
nbins <- 40

smudge_container_sample <- get_smudge_container(cov_tab_sample, nbins, xlim, ylim)

args <- list()
args$col_ramp <- "viridis"
args$invert_cols <- FALSE
colour_ramp <- get_col_ramp(args) # get the default colour ramp (Spectral, 11)

plot_smudgeplot(smudge_container_sample, 15.5, colour_ramp)

######## SCER

cov_tab <- read.table("data/Scer/PloidyPlot3_text.smu", col.names = c('covB', 'covA', 'freq')) #nolint
cov_tab[, 'total_pair_cov'] <- cov_tab[, 1] + cov_tab[, 2]
cov_tab[, 'minor_variant_rel_cov'] <- cov_tab[, 1] / cov_tab[, 'total_pair_cov']

smudge_container <- get_smudge_container(cov_tab, nbins, xlim, ylim)
plot_smudgeplot(smudge_container, 15.5, colour_ramp)
rect(0.4375, smudge_container$y[6], 0.4500, smudge_container$y[7], col = 'white')
rect(0, 15.5*2, 0.500, 15.5*2, col = 'grey')
rect(0.49, 0, 0.49, 475, col = 'grey')
rect(0.4935, 0, 0.49, 475, col = 'grey')

########

L <- 12 
y_min <- 22.50  
y_max <- 26.25


bracket_size <- diff(smudge_container$y)[1]

allowed_ranges <- lapply(smudge_container$y, function(ymin){if(ymin + bracket_size <= 2*L){ return(c()) }; seq(max(ceiling(ymin), 2*L), ymin + bracket_size, by = 1)})

# this will generate a table of all possible k-mer combinations of covA and covB
impossible_ratios <- function(allowed_range, L, x_brackets){
  impossible <- rep(TRUE, length(smudge_container$x))
  if (is.null(allowed_range)) {
    return (impossible)
  }
  minicov_tab <- merge(L:floor(max(allowed_range) / 2), (ceiling(max(allowed_range) / 2):(max(allowed_range) - L)))
  colnames(minicov_tab) <- c("covB", "covA")
  minicov_tab <- minicov_tab[minicov_tab[, 'covB'] <= minicov_tab[, 'covA'], ] 
  minicov_tab[, 'cov_ratio'] <- minicov_tab[, 'covB'] / (minicov_tab[, 'covB'] + minicov_tab[, 'covA'])
  minicov_tab[, 'bracket'] <- findInterval(minicov_tab[, 'cov_ratio'], x_brackets, left.open = TRUE)
  xcoords <- unique(minicov_tab[, 'bracket'])  
  impossible[xcoords] <- FALSE
  return (impossible)
}

rect(0.4375, smudge_container$y[6], 0.4500, smudge_container$y[7], col = 'white')

rect(0, 0, 0.5, 10, col = 'white')
rect(0, 0, 0.05, 250, col = 'white')

allowed_range <- allowed_ranges[[8]]
impossible_ratios(allowed_ranges[[7]], L, smudge_container$x)
impossible_ratios(allowed_ranges[[8]], L, smudge_container$x)
matrix_of_impossibles <- sapply(allowed_ranges, impossible_ratios, L, smudge_container$x)

smudge_container <- get_smudge_container(cov_tab, nbins, xlim, ylim)
smudge_container$z[matrix_of_impossibles] <- NA
plot_smudgeplot(smudge_container, 15.5, colour_ramp)


L <- 12
cov_tab_mini <- merge(c(L:20), c(L:20))
colnames(cov_tab_mini) <- c("covA", "covB")
cov_tab_mini <- cov_tab_mini[cov_tab_mini$covB <= cov_tab_mini$covA, ] 
# now cov_tab_mini is a realistic cov table

cov_tab_mini$cov_sum <- cov_tab_mini$covB + cov_tab_mini$covA
cov_tab_mini$cov_ratio <- cov_tab_mini$covB / cov_tab_mini$cov_sum

cov_tab_mini_first_covsum_bracket <- cov_tab_mini[cov_tab_mini$cov_sum %in% allowed_range, ]
table(cov_tab_mini_first_covsum_bracket$cov_ratio)


library(smudgeplot)

args <- ArgumentParser()$parse_args()
args$output <- "data/Scer/sploidyplot_test"
args$nbins <- 40
args$kmer_size <- 21
args$homozygous <- FALSE
args$L <- c()
args$col_ramp <- 'viridis'
args$invert_cols <- TRUE


  L <- 12
  smudge_container$x
  smudge_container$y
  is_possible <- function(x_min, x_max, y_min, y_max, L){
    y_min <- 22.50  
    y_max <- 26.25
    seq(ceiling(y_min), y_max, by = 1)
    x_min <- 0.4500
    x_max <- 0.4625
    if (x_max < L / y_max)
  }


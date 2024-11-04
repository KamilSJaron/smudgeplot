library(smudgeplot)

args <- ArgumentParser()$parse_args()
args$output <- "data/Scer/sploidyplot_test"
args$nbins <- 40
args$kmer_size <- 21
args$homozygous <- FALSE
args$L <- c()
args$col_ramp <- 'viridis'
args$invert_cols <- TRUE

cov_tab <- read.table("data/Scer/PloidyPlot3_text.smu", col.names = c('covB', 'covA', 'freq'), skip = 2) #nolint
cov_tab[, 'total_pair_cov'] <- cov_tab[, 1] + cov_tab[, 2]
cov_tab[, 'minor_variant_rel_cov'] <- cov_tab[, 1] / cov_tab[, 'total_pair_cov']

cov_tab_1n_est <- round(estimate_1n_coverage_1d_subsets(cov_tab), 1)

xlim <- c(0, 0.5)
# max(total_pair_cov); 10*draft_n
ylim <- c(0, 150)
nbins <- 40

smudge_container <- get_smudge_container(cov_tab, nbins, xlim, ylim)
smudge_container$z <- smudge_container$dens

plot_popart <- function(cov_tab, ylim, colour_ramp){
    A_equals_B <- cov_tab[, 'covA'] == cov_tab[, 'covB']
    cov_tab[A_equals_B, 'freq'] <- cov_tab[A_equals_B, 'freq'] * 2
    cov_tab$col <- colour_ramp[1 + round(31 * cov_tab[, 'freq'] / max(cov_tab[, 'freq']))]

    plot(NULL, xlim = c(0.1, 0.5), ylim = ylim, xaxt = "n", yaxt = "n", xlab = '', ylab = '', bty = 'n')
    min_cov_to_plot <- max(ylim[1],min(cov_tab[, 'total_pair_cov']))
    nothing <- sapply(min_cov_to_plot:ylim[2], plot_one_coverage, cov_tab)
    return(0)
}

par(mfrow = c(2, 5))

args$col_ramp <- "viridis"
args$invert_cols <- FALSE
colour_ramp <- get_col_ramp(args) # get the default colour ramp (Spectral, 11)
# plot_smudgeplot(smudge_container, 15.5, colour_ramp)
plot_popart(cov_tab, c(20, 120), colour_ramp)

args$invert_cols <- TRUE
colour_ramp <- get_col_ramp(args) # get the default colour ramp (Spectral, 11)
# plot_smudgeplot(smudge_container, 15.5, colour_ramp)
plot_popart(cov_tab, c(20, 120), colour_ramp)


for (ramp in c("grey.colors", "magma", "plasma", "mako", "inferno", "rocket", "heat.colors", "cm.colors")){
    args$col_ramp <- ramp
    colour_ramp <- get_col_ramp(args) # get the default colour ramp (Spectral, 11)
    plot_popart(cov_tab, c(20, 120), colour_ramp)
}


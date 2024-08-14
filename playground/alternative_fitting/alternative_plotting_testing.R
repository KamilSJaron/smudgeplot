library(smudgeplot)
library(argparse)
source('playground/alternative_fitting/alternative_plotting_functions.R')

parser <- ArgumentParser()
parser$add_argument("-i", "-infile", help="Input file")
parser$add_argument("-o", "-outfile", help="Output file")
args <- parser$parse_args()

args$col_ramp <- 'viridis'
args$invert_cols <- F

# cov_tab_daAchMill1 <- read.table('data/dicots/smu_files/daAchMill1.k31_ploidy.smu.txt', col.names = c('covB', 'covA', 'freq'))
# cov_tab <- read.table(args$file, col.names = c('covB', 'covA', 'freq'))
# args <- list()
# args$i <- 'data/ddSalArbu1/ddSalArbu1.k31_ploidy_converted.smu_with_peaks.txt'
# args$o <- 'data/ddSalArbu1/smudge_with_peaks'

cov_tab <- read.table(args$i, col.names = c('covB', 'covA', 'freq','peak'))

xlim = c(0, 0.5)
ylim = c(0, 150)


cov_tab[, 'total_pair_cov'] <- cov_tab[, 1] + cov_tab[, 2]
cov_tab[, 'minor_variant_rel_cov'] <- cov_tab[, 1] / cov_tab[, 'total_pair_cov']

colour_ramp <- viridis(32)
colour_ramp_log <- get_col_ramp(args, 16)
# cols = colour_ramp[1 + round(31 * cov_tab$freq / max(cov_tab$freq))]

# solving this "density" problem; cov1 cov1; have twice lower probability than cov1 cov2; we need to multiply these points, but it needs to be somehow corrected for the fit / summaries

plot_dot_smudgeplot(cov_tab, colour_ramp, xlim, ylim)

pdf(paste0(args$o, '_background.pdf'))
    plot_alt(cov_tab, ylim, colour_ramp, F)
dev.off()

pdf(paste0(args$o, '_log_background.pdf'))
    plot_alt(cov_tab, ylim, colour_ramp, T)
dev.off()

pdf(paste0(args$o, '_peaks.pdf'))
    plot_peakmap(cov_tab, xlim = xlim, ylim = ylim)
dev.off()

# plots the line where there will be nothing
# plot_seq_error_line(cov_tab, .col = 'red')

# head(cov_tab[order(cov_tab[,'total_pair_cov']), ], 12)
# colour_ramp
# 3 / 8:13

####

# cov_tab_Fiin_ideal <- read.table('data/Fiin/idealised/kmerpairs_idealised_with_transformations.tsv', header = T)
# head(cov_tab_Fiin_ideal)

# xlim = c(0, 0.5)
# ylim = c(0, max(cov_tab_Fiin_ideal[, 'total_pair_cov']))

# pdf('data/Fiin/idealised/straw_plot_test3.pdf')
#     plot_dot_smudgeplot(cov_tab_Fiin_ideal, colour_ramp, xlim, ylim)
# dev.off()

# pdf('data/Fiin/idealised/straw_plot_test2.pdf')
#     plot_unsquared_smudgeplot(cov_tab_Fiin_ideal, colour_ramp, xlim, ylim)
# dev.off()

# # testing the packaged version

# library(smudgeplot)
# args <- list()
# args$col_ramp <- 'viridis'
# args$invert_cols <- F
# colour_ramp <- get_col_ramp(args)

# cov_tab_Fiin_ideal <- read.table('data/Fiin/idealised/kmerpairs_idealised_with_transformations.tsv', header = T)

# pdf('data/Fiin/idealised/straw_plot_test.pdf')
#     plot_alt(cov_tab_Fiin_ideal, c(50, 700), colour_ramp)
# dev.off()

# source('playground/alternative_fitting/alternative_plotting_functions.R')

# pdf('data/Fiin/idealised/straw_plot_test2.pdf')
#     plot_unsquared_smudgeplot(cov_tab_Fiin_ideal, colour_ramp, c(0, 0.5), c(50, 700))
# dev.off()


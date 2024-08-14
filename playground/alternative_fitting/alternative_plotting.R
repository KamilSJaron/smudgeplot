# HM Revenue and custom Tax Return form;
# Needs to be 2023 Jan to 2024 Jan



library(smudgeplot)
source('playground/alternative_plotting_functions.R')

cov_tab_daAchMill1 <- read.table('data/dicots/smu_files/daAchMill1.k31_ploidy.smu.txt', col.names = c('covB', 'covA', 'freq'))
ylim = c(0, 250)

xlim = c(0, 0.5)



cov_tab_daAchMill1[, 'total_pair_cov'] <- cov_tab_daAchMill1[, 1] + cov_tab_daAchMill1[, 2]
cov_tab_daAchMill1[, 'minor_variant_rel_cov'] <- cov_tab_daAchMill1[, 1] / cov_tab_daAchMill1[, 'total_pair_cov']

args <- list()
args$col_ramp <- 'viridis'
args$invert_cols <- F
colour_ramp <- get_col_ramp(args)
cols = colour_ramp[1 + round(31 * cov_tab_daAchMill1$freq / max(cov_tab_daAchMill1$freq))]

# solving this "density" problem; cov1 cov1; have twice lower probability than cov1 cov2; we need to multiply these points, but it needs to be somehow corrected for the fit / summaries

plot_dot_smudgeplot(cov_tab_daAchMill1, colour_ramp, xlim, ylim)

plot_unsquared_smudgeplot(cov_tab_daAchMill1, colour_ramp, xlim, ylim)

# plots the line where there will be nothing
plot_seq_error_line(cov_tab_daAchMill1, .col = 'red')

head(cov_tab_daAchMill1[order(cov_tab_daAchMill1[,'total_pair_cov']), ], 12)
colour_ramp
3 / 8:13

####

cov_tab_Fiin_ideal <- read.table('data/Fiin/idealised/kmerpairs_idealised_with_transformations.tsv', header = T)
head(cov_tab_Fiin_ideal)

xlim = c(0, 0.5)
ylim = c(0, max(cov_tab_Fiin_ideal[, 'total_pair_cov']))

plot_dot_smudgeplot(cov_tab_Fiin_ideal, colour_ramp, xlim, ylim)

plot_unsquared_smudgeplot(cov_tab_Fiin_ideal, colour_ramp, xlim, ylim)

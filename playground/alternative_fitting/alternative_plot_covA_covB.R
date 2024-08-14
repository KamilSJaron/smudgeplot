library(ggplot2)

plot_unsquared_smudgeplot <- function(cov_tab, colour_ramp, xlim, ylim){
    # this is the adjustment for plotting
    # cov_tab[cov_tab$covA == cov_tab$covB, 'freq'] <- cov_tab[cov_tab$covA == cov_tab$covB, 'freq'] * 2
    cov_tab$col = colour_ramp[1 + round(31 * cov_tab$freq / max(cov_tab$freq))]

    plot(NULL, xlim = xlim, ylim = ylim,
         xlab = 'covA',
         ylab = 'covB', cex.lab = 1.4)

    ## This might bite me in the a.., instead of taking L as an argument, I guess it from the data
    # L = floor(min(cov_tab_daAchMill1[, 'total_pair_cov']) / 2)
    ggplot(cov_tab, aes(x=covA, y=covB, weight = weight)) +
    geom_bin2d() +
    theme_bw()

}

real_clean <- read.table('data/Fiin/kmerpairs_k51_text.smu', col.names = c('covA', 'covB', 'freq'))
real_clean$weight <- real_clean$freq / sum(real_clean$freq)

# plot(real_clean[, 'covA'], real_clean[, 'covB'])

xlim <- range(real_clean[, 'covA'])
ylim <- range(real_clean[, 'covB'])

library(smudgeplot)
args <- list()
args$col_ramp <- 'viridis'
args$invert_cols <- F
colour_ramp <- get_col_ramp(args)
real_clean$col <- colour_ramp[1 + round(31 * real_clean$freq / max(real_clean$freq))]

plot(NULL, xlim = xlim, ylim = ylim,
     xlab = 'covA',
     ylab = 'covB', cex.lab = 1.4)

ggplot(real_clean, aes(x=covA, y=covB, weight = weight)) +
  geom_bin2d() +
  theme_bw()

head(real_clean)


# plotSquare <- function(row){
#     rect(as.numeric(row['covA']) - 0.5, as.numeric(row['covB']) - 0.5, as.numeric(row['covA']) + 0.5, as.numeric(row['covB']) + 0.5, col = row['col'], border = NA)
# }
# apply(real_clean, 1, plotSquare)

files <- c('data/Avag1/coverages_2.tsv',
           'data/Lcla1/Lcla1_pairs_coverages_2.tsv',
           'data/Mflo2/coverages_2.tsv',
           'data/Rvar1/Rvar1_pairs_coverages_2.tsv',
           'data/Ps791/Ps791_pairs_coverages_2.tsv',
           'data/Aric1/Aric1_pairs_coverages_2.tsv',
           "data/Rmag1/Rmag1_pairs_coverages_2.tsv")

###
library(smudgeplot)
args <- list()
args$input <- 'data/Mflo2/Mflo2_coverages_2.tsv'
args$output <- "figures/Mflo2_v0.1.0"
args$nbins <- 40
args$kmer_size <- 21
args$homozygous <- F

###

i <- 7
n <- NA
cov <- read.table(files[i])

# run bits of smudgeplot.R to get k, and peak summary

filter <- total_pair_cov < 350
total_pair_cov_filt <- total_pair_cov[filter]
minor_variant_rel_cov_filt <- minor_variant_rel_cov[filter]

ymax <- max(total_pair_cov_filt)
ymin <- min(total_pair_cov_filt)

# the lims trick will make sure that the last column of squares will have the same width as the other squares
smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov, .nbins = 40)

x <- seq(xlim[1], ((nbins - 1) / nbins) * xlim[2], length = nbins)
y <- c(seq(ylim[1] - 0.1, ((nbins - 1) / nbins) * ylim[2], length = nbins), ylim[2])

.peak_points <- peak_points
.smudge_container <- smudge_container
.total_pair_cov <- total_pair_cov
.treshold <- 0.05
fig_title <- 'test'

image(smudge_container, col = colour_ramp)
# contour(x.bin, y.bin, freq2D, add=TRUE, col=rgb(1,1,1,.7))

#######
# PLOT
#######

library(plotly)
packageVersion('plotly')

p <- plot_ly(x = k_toplot$x, y = k_toplot$y, z = k_toplot$z) %>% add_surface()
htmlwidgets::saveWidget(p, "Ps791_smudge_surface.html")
# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = api_create(p, filename="Ps791_smudge_surface-2")
chart_link

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(k_toplot, n, colour_ramp)
plot_expected_haplotype_structure(n, peak_sizes, T)
# annotate_peaks(peak_points, ymin, ymax)
# annotate_summits(peak_points, peak_sizes, ymin, ymax, 'black')
# TODO fix plot_seq_error_line(total_pair_cov)
# 2,3 hist
# TODO rescale histogram axis by the scale of the smudgeplot
plot_histograms(minor_variant_rel_cov, total_pair_cov)
# 4 legend
plot_legend(k_toplot, total_pair_cov, colour_ramp)

# findInterval(c(0.1, 0.2, 0.33, 0.5), seq(0, 0.5, length = 41))

##########################################################
## TEST
## idea here was to propagate from the highest point and expand the peak till it's monotonic
# starting_point <- which(dens_m == max(dens_m), arr.ind = TRUE)
# starting_val <- dens_m[starting_point]
# peak_points <- data.frame(x = starting_point[,2], y = starting_point[,1], value = starting_val)
#
# points_to_explore <- get_neibours(starting_val)
# val_to_comp <- starting_val
#
# for(one_point in 1:nrow(points_to_explore)){
#     one_point <- points_to_explore[one_point,]
#     point_val <- dens_m[t(one_point)]
#     if(point_val < val_to_comp){
#         peak_points <- rbind(peak_points,
#             data.frame(x = one_point[2], y = one_point[1], value = point_val))
#     }
# }
#
# get_neibours <- function(point){
#     neibours_vec <- matrix(rep(starting_point,8) + c(-1,-1,0,-1,1,-1,-1,0,+1,0,-1,1,0,1,1,1),
#                            ncol = 2, byrow = T)
#     neibours_vec[rowSums(neibours_vec <= 30 & neibours_vec >= 1) == 2,]
# }
#

##########################
###ALTERNATIVE PLTTING ###
##########################
# library('spatialfil')
# msnFit(high_cov_filt, minor_variant_rel_cov)

## alternative plotting
# library(hexbin) # honeycomb plot
# h <- hexbin(df)
# # h@count <- sqrt(h@count)
# plot(h, colramp=rf)
# gplot.hexbin(h, colorcut=10, colramp=rf)


### TEST  plot lines at expected coverages
#
# for(i in 2:6){
#       lines(c(0, 0.6), rep(i * n, 2), lwd = 1.4)
#       text(0.1, i * n, paste0(i,'x'), pos = 3)
# }

# FUTURE - wrapper
# smudgeplot < - function(.k, .minor_variant_rel_cov, .total_pair_cov, .n,
#                             .sqrt_scale = T, .cex = 1.4, .fig_title = NA){
#     if( .sqrt_scale == T ){
#         # to display densities on squared root scale (bit like log scale but less agressive)
#         .k$z <- sqrt(.k$z)
#     }
#
#     pal <- brewer.pal(11,'Spectral')
#     rf <- colorRampPalette(rev(pal[3:11]))
#     colour_ramp <- rf(32)
#
#     layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
#
#     # 2D HISTOGRAM
#     plot_smudgeplot(...)
#
#     # 1D historgram - minor_variant_rel_cov on top
#     plot_histogram(...)
#
#     # 1D historgram - total pair coverage - right
#     plot_histogram(...)
#
#     # LEGEND (topright corener)
#     plot_legend(...)
#
# }

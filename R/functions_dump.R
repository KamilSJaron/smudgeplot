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

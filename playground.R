files <- c('data/Avag1/coverages_2.tsv',
           'data/Lcla1/Lcla1_pairs_coverages_2.tsv',
           'data/Mflo2/coverages_2.tsv',
           'data/Rvar1/Rvar1_pairs_coverages_2.tsv',
           'data/Ps791/Ps791_pairs_coverages_2.tsv')

cov <- lapply(files, read.table)
i <- 5

######
# INPUT PROCESSING
######
# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov[[i]]$V1 / (cov[[i]]$V1 + cov[[i]]$V2)
# total covarate of the kmer pair
total_pair_cov <- cov[[i]]$V1 + cov[[i]]$V2

# source all the functions

high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

######
# SUMMARY CALCULATION
######

# estimation method 1 - using weighted mean of all detected peaks and their expected ploidy
n <- round(estimate_1n_coverage(),1)

# estimation method 2 - using global maxima in the 2d plot
k <- kde2d(minor_variant_rel_cov, total_pair_cov, n=30,
           lims = c(0.02, 0.48, min(total_pair_cov), min(10*n, max(total_pair_cov))))

peak_points <- peak_agregation(k)
peak_sizes <- get_peak_sizes(peak_points)

# filt_peak_points <- filter_peaks(peak_points, peak_sizes)
# filt_peak_sizes <- filter_peak_sizes(peak_sizes)

#######
# PLOT
#######

k_toplot <- k
k_toplot$z <- sqrt(k_toplot$z)

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(k_toplot, n, colour_ramp)
plot_expected_haplotype_structure(n)
annotate_peaks(peak_points, total_pair_cov)
annotate_summits(peak_points, peak_sizes, total_pair_cov)
# TODO fix plot_seq_error_line(total_pair_cov)
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov)
# 4 legend
plot_legend(k_toplot, total_pair_cov, colour_ramp)


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

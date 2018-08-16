files <- c('data/Avag1/coverages_2.tsv',
           'data/Lcla1/Lcla1_pairs_coverages_2.tsv',
           'data/Mflo2/coverages_2.tsv',
           'data/Rvar1/Rvar1_pairs_coverages_2.tsv')

cov <- lapply(files, read.table)
i <- 3
# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov[[i]]$V1 / (cov[[i]]$V1 + cov[[i]]$V2)
# total covarate of the kmer pair
total_pair_cov <- cov[[i]]$V1 + cov[[i]]$V2

high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

# estimation method 1 - using weighted mean of all detected peaks and their expected ploidy
n <- round(estimate_1n_coverage(),1)

# estimation method 2 - using global maxima in the 2d plot
k <- kde2d(minor_variant_rel_cov, total_pair_cov, n=30,
           lims = c(0.02, 0.48, min(total_pair_cov), max(total_pair_cov)))

peak_points <- data.frame(vals = as.vector(k$z), x = rep(1:30, each = 30), y = rep(1:30, 30))
peak_points$peak <- NA
peak_points$summit <- NA
peak_points <- peak_points[order(peak_points$vals, decreasing = T),]
# mark the first
peak_points$peak[1] <- 1
peak_points$summit[1] <- T

for(point in 2:nrow(peak_points)){
    parsed_points <- peak_points[1:point-1,]
    x_n <- abs(parsed_points[,'x'] - peak_points[point,'x']) <= 1
    y_n <- abs(parsed_points[,'y'] - peak_points[point,'y']) <= 1
    if(any(x_n & y_n)){
        tiles_around <- parsed_points[x_n & y_n,]
        peak_points$peak[point] <- tiles_around$peak[which.max(tiles_around$vals)]
        peak_points$summit[point] <- F
    } else {
        peak_points$peak[point] <- max(parsed_points$peak) + 1
        peak_points$summit[point] <- T
    }
}

peak_sizes <- c()
for (peak in unique(peak_points$peak)){
    peak_sizes <- c(peak_sizes, sum(peak_points[peak_points$peak == peak, 'vals'], na.rm = T))
}
peak_sizes <- peak_sizes / sum(peak_sizes)

to_filter <- peak_points$peak %in% which(peak_sizes < 0.005)
peak_sizes <- peak_sizes[!peak_sizes < 0.005]

peak_points$peak[to_filter] <- NA
peak_points$summit[to_filter] <- F
# unique(peak_points$peak)

peak_points <- peak_points[with(peak_points, order(x,y)),]
# still works if x goes 1 1 ... 1 1 2 2 ... 2 2 3 3... 3 3 4 4 4 ... 30 30 30
# and y is              1 2 3 4 5 ... 30 1 2 3 4 5 ... 30 ... 1 2 3 4  ... 30
image(matrix(peak_points$vals^0.25, ncol = 30))

for (summit in which(peak_points$summit == T)){
     peak <- peak_points$peak[summit]
     text(peak_points$y[summit] / 30 ,
          peak_points$x[summit] / 30,
          round(peak_sizes[peak], 3),
          pos = 2)

    to_plot <- peak_points[peak_points$peak == peak, c('x','y')]
    to_plot <- to_plot[!is.na(to_plot$x),]
    points(to_plot$y / 30, to_plot$x / 30, pch = peak)
}

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

#########################
####### FUNCTIONS #######
#########################


#### APPROACH ONE ####
# estimate 1n coverage by smothing stripes of points at minor_variant_rel_cov 0.5, 0.33, 0.25, 0.2
# and assuming that the firt 0.5 is the diploid peak (which is the problem I guess)
# get local maxima
get_1d_peaks <- function(subset, num_of_peaks = 3, adjust = 10, peak_frame = data.frame(), depth = 1){
    if (nrow(peak_frame) < num_of_peaks & depth < 11){
        d <- density(total_pair_cov[subset], adjust = adjust) # returns the density data
        # plot(d)
        selected_points <- which(diff(sign(diff(d$y))) == -2)
        peak_frame <- data.frame(cov = d$x[selected_points], height = d$y[selected_points] * sum(subset))
        peak_frame <- get_1d_peaks(subset, num_of_peaks, adjust - 1, peak_frame, depth + 1)
    }
    peak_frame
}

estimate_1n_coverage <- function(){
    total_kmers <- length(minor_variant_rel_cov)
    minor_freq_subsets <- list(
        AB_subset = minor_variant_rel_cov > 0.47,
        AAB_subset = minor_variant_rel_cov < 0.35 & minor_variant_rel_cov > 0.31,
        AAAB_subset = minor_variant_rel_cov < 0.27 & minor_variant_rel_cov > 0.23,
        AAAAB_subset = minor_variant_rel_cov < 0.21 & minor_variant_rel_cov > 0.19
    )

    peak_frame_2 <- get_1d_peaks(minor_freq_subsets[[1]], 3)
    peak_frame_3 <- get_1d_peaks(minor_freq_subsets[[2]], 2)
    peak_frame_4 <- get_1d_peaks(minor_freq_subsets[[3]], 1)
    peak_frame_5 <- get_1d_peaks(minor_freq_subsets[[4]], 1)

    peak_frame_2$cov <- peak_frame_2$cov / (2 * round(peak_frame_2$cov / peak_frame_2$cov[1]))
    peak_frame_3$cov <- peak_frame_3$cov / (3 * round(peak_frame_3$cov / peak_frame_3$cov[1]))
    peak_frame_4$cov <- peak_frame_4$cov / (4 * round(peak_frame_4$cov / peak_frame_4$cov[1]))
    peak_frame_5$cov <- peak_frame_5$cov / (5 * round(peak_frame_5$cov / peak_frame_5$cov[1]))

    peak_frame <- rbind(peak_frame_2, peak_frame_3, peak_frame_4, peak_frame_5)
    weighted.mean(peak_frame$cov, peak_frame$height)
}

#### APPROACH TWO #####
# propagation of clusters from summits of the 2d histogram

### LATEST peak_sizes function

get_peak_sizes <- function(.peak_points){
    peaks <- .peak_points[peak_points$summit == T,'peak']
    summits <- .peak_points[peak_points$summit == T,'vals']
    peak_sizes <- sapply(peaks, function(x){sum(.peak_points[.peak_points$peak == x, 'vals'])})
    data.frame(peak = peaks,
                abs_size = peak_sizes,
                rel_size = peak_sizes / sum(peak_sizes),
                summit_height = summits)
}

### LATEST peak agregation function

peak_agregation <- function(.k){
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
    return(peak_points)
}

filter_peaks <- function(.peak_points, .peak_sizes, threshold = 0.005){
    to_filter <- .peak_points$peak %in% which(.peak_sizes$rel_size < threshold)

    .peak_points$peak[to_filter] <- NA
    .peak_points$summit[to_filter] <- F
    return(.peak_points)
}

filter_peak_sizes <- function(.peak_sizes, threshold = 0.005){
    .peak_sizes[!peak_sizes$rel_size < threshold,]
}

get_peak_annotation <- function(peak_points, peak_sizes){

}


# to transform 1 ... 30 coordinates to x coordinates of the 2d hisrtogram
transform_x <- function(x, for_plot = T){
    orig <- (((x - 1) / 29) / 2)
    if ( for_plot ){
        ((orig / 0.5) * 0.46) + 0.02
    } else {
        orig
    }
}

# to transform 1 ... 30 coordinates to y coordinates of the 2d hisrtogram
transform_y <- function(y, min_kmerpair_cov, max_cov_pair){
    (((y - 1) / 29) * (max_cov_pair - min_kmerpair_cov)) + min_kmerpair_cov
}

#
# VISUALIZATION
#  - smudgeplot: weapper
#
#  - plot_smudgeplot: core function for 2d histogram
#  - annotate_peaks:

# FUTURE - wrapper
# smudgeplot <- function(.k, .minor_variant_rel_cov, .total_pair_cov, .n,
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

# SMUDGEPLOT - 2D histogram - landscape plot - howeveryouwantyoucallthis plot
# makes a plot
plot_smudgeplot <- function(.k, .n, .colour_ramp, .cex = 1.4){
    # margins 'c(bottom, left, top, right)'
    par(mar=c(4.8,4.8,1,1))
    # 2D HISTOGRAM
    image(.k, col = colour_ramp,
        xlab = 'Normalized minor kmer coverage: B / (A + B)',
        ylab = 'Total coverage of the kmer pair: A + B', cex.lab = 1.4,
        axes=F
    )

    axis(2, at = 2:8 * .n, labels = paste(2:8, 'n'))
    axis(1, at = c(1/5, 1/4, 1/3, 2/5, 0.487),
            labels = c('1:4', '1:3', '1:2', '2:3', '1:1'))
}

# EXPECTED COMPOSITIONS - bettern than lines
# adds to a plot
# TODO do automatically till y max? or max plody?
plot_expected_haplotype_structure <- function(.n, .cex = 1.4){
    text(1/2 - 0.027, 2 * .n, 'AB', offset = 0, cex = .cex)
    text(1/3, 3 * .n, 'AAB', offset = 0, cex = .cex)
    text(1/4, 4 * .n, 'AAAB', offset = 0, cex = .cex)
    text(1/2 - 0.04, 4 * .n, 'AABB', offset = 0, cex = .cex)
    text(2/5, 5 * .n, 'AAABB', offset = 0, cex = .cex)
    text(1/5, 5 * .n, 'AAAAB', offset = 0, cex = .cex)
    text(3/6 - 0.055, 6 * .n, 'AAABBB', offset = 0, cex = .cex)
    text(2/6, 6 * .n, 'AAAABB', offset = 0, cex = .cex)
    text(1/6, 6 * .n, 'AAAAAB', offset = 0, cex = .cex)
}

# plot HISTOGRAMS on mrgins
# makes 2 plots
plot_histograms <- function(.minor_variant_rel_cov, .total_pair_cov,
                            .fig_title = NA, .cex = 1.4){
    h1 <- hist(.minor_variant_rel_cov, breaks = 100, plot = F)
    h2 <- hist(.total_pair_cov, breaks = 100, plot = F)
    top <- max(h1$counts, h2$counts)

    # minor_variant_rel_cov HISTOGRAM - top
    par(mar=c(0,3.8,1,0))
    barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col = pal[2])
    if(!(is.na(.fig_title))){
        mtext(bquote(italic(.(.fig_title))), side=3, adj=0, line=-3, cex = .cex + 0.2)
    }

    # total pair coverage HISTOGRAM - right
    par(mar=c(3.8,0,0.5,1))
    barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col = pal[2], horiz = T)
    mtext(paste('1n = ', n), side=1, adj=0.8, line=-2, cex = .cex)
}

# LEGEND (topright corener in the main plot)
# makes a plot
plot_legend <- function(.k, .total_pair_cov, .colour_ramp, .sqrt_scale = T){
    par(mar=c(0,0,2,1))
    plot.new()
    title('kmers pairs')
    for(i in 1:32){
        rect(0,(i - 0.01) / 33, 0.5, (i + 0.99) / 33, col = .colour_ramp[i])
    }
    if( .sqrt_scale == T ){
        # TODO correct this scale
        kmer_max <- (length(.total_pair_cov) * max((.k$z)^2)) / sum((.k$z)^2)
        for(i in 0:6){
            text(0.75, i / 6, round((sqrt(kmer_max) * i)^2 / 6000) * 1000, offset = 0)
        }
    } else {
        kmer_max <- (length(.total_pair_cov) * max(.k$z)) / sum(.k$z)
        for(i in 0:6){
            text(0.75, i / 6, round((kmer_max * i) / 6000) * 1000, offset = 0)
        }
    }
}

# add error line
# TODO it's misplaced a bit I think, fix needed
plot_seq_error_line <- function(total_pair_cov){
    min_kmer_cov <- min(total_pair_cov) / 2
    max_cov_pair <- max(total_pair_cov)
    cov_range <- seq(2 * min_kmer_cov,max_cov_pair, length = 30)
    lines(min_kmer_cov / cov_range, cov_range, lwd = 2, lty = 2)
}

# annotate peaks
annotate_peaks <- function(peak_points, min_kmerpair_cov, max_kmerpair_cov){
    for (summit in which(peak_points$summit == T)){
        peak <- peak_points$peak[summit]
        to_plot <- peak_points[peak_points$peak == peak, c('x','y')]
        to_plot <- to_plot[!is.na(to_plot$x),]
        points(transform_x(to_plot$y), transform_y(to_plot$x, min_kmerpair_cov, max_kmerpair_cov), pch = peak)
    }
}

# annotate summits
annotate_summits <- function(peak_points, peak_sizes, min_kmerpair_cov, max_kmerpair_cov, col = 'red'){
    for (summit in which(peak_points$summit == T)){
         peak <- peak_points$peak[summit]
         peak_size <- peak_sizes$rel_size[peak_sizes$peak == peak]
         # points(transform_x(peak_points$y[summit]),
         #        transform_y(peak_points$x[summit], min_kmerpair_cov, max_kmerpair_cov))
         text(transform_x(peak_points$y[summit]),
              transform_y(peak_points$x[summit], min_kmerpair_cov, max_kmerpair_cov),
              round(peak_size, 3),
              offset = 0, cex = 1.2, col = col, xpd=NA)
    }
}



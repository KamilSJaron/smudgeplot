#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("viridis"))

# suppressPackageStartupMessages(library("smudgeplot"))

#################
### funcitons ###
#################
# retirying the smudgeplot R package
get_col_ramp <- function(.args, delay = 0){
    colour_ramp <- eval(parse(text = paste0(.args$col_ramp,"(", 32 - delay, ")")))
    if (.args$invert_cols){
        colour_ramp <- rev(colour_ramp)
    }
    colour_ramp <- c(rep(colour_ramp[1], delay), colour_ramp)
    return(colour_ramp)
}

wtd.quantile <- function(x, q=0.25, weight=NULL) {
  o <- order(x)
  n <- sum(weight)
  order <- 1 + (n - 1) * q
  low  <- pmax(floor(order), 1)
  high <- pmin(ceiling(order), n)
  low_contribution <- high - order
  allq <- approx(x=cumsum(weight[o])/sum(weight), y=x[o], xout = c(low, high)/n, method = "constant",
      f = 1, rule = 2)$y
  low_contribution * allq[1] + (1 - low_contribution) * allq[2]
}

wtd.iqr <- function(x, w=NULL) {
  wtd.quantile(x, q=0.75, weight=w) - wtd.quantile(x, q=0.25, weight=w)
}

plot_alt <- function(cov_tab, ylim, colour_ramp, log = F){
    A_equals_B <- cov_tab[, 'covA'] == cov_tab[, 'covB']
    cov_tab[A_equals_B, 'freq'] <- cov_tab[A_equals_B, 'freq'] * 2
    if (log){
        cov_tab[, 'freq'] <- log10(cov_tab[, 'freq'])
    }
    cov_tab$col <- colour_ramp[1 + round(31 * cov_tab[, 'freq'] / max(cov_tab[, 'freq']))]

    # c(bottom, left, top, right)
    par(mar=c(4.8,4.8,1,1))
    plot(NULL, xlim = c(0, 0.5), ylim = ylim,
         xlab = 'Normalized minor kmer coverage: B / (A + B)',
         ylab = 'Total coverage of the kmer pair: A + B', cex.lab = 1.4, bty = 'n')
    min_cov_to_plot <- max(ylim[1],min(cov_tab[, 'total_pair_cov']))
    nothing <- sapply(min_cov_to_plot:ylim[2], plot_one_coverage, cov_tab)
    return(0)
}

plot_one_coverage <- function(cov, cov_tab){
    cov_row_to_plot <- cov_tab[cov_tab[, 'total_pair_cov'] == cov, ]
    width <- 1 / (2 * cov)
    cov_row_to_plot$left <- cov_row_to_plot[, 'minor_variant_rel_cov'] - width
    cov_row_to_plot$right <- sapply(cov_row_to_plot[, 'minor_variant_rel_cov'], function(x){ min(0.5, x + width)})
    apply(cov_row_to_plot, 1, plot_one_box, cov)
}

plot_one_box <- function(one_box_row, cov){
    left <- as.numeric(one_box_row['left'])
    right <- as.numeric(one_box_row['right'])
    rect(left, cov - 0.5, right, cov + 0.5, col = one_box_row['col'], border = NA)
}

plot_isoA_line <- function (.covA, .L, .col = "black", .ymax = 250, .lwd, .lty) {
    min_covB <- .L # min(.cov_tab[, 'covB']) # should be L really
    max_covB <- .covA
    B_covs <- seq(min_covB, max_covB, length = 500)
    isoline_x <- B_covs/ (B_covs + .covA)
    isoline_y <- B_covs + .covA
    lines(isoline_x[isoline_y < .ymax], isoline_y[isoline_y < .ymax], lwd = .lwd, lty = .lty, col = .col)
}

plot_isoB_line <- function (.covB, .ymax, .col = "black", .lwd, .lty) {
    cov_range <- seq((2 * .covB) - 2, .ymax, length = 500)
    lines((.covB)/cov_range, cov_range, lwd = .lwd, lty = .lty, col = .col)
}

plot_iso_grid <- function(.cov, .L, .ymax, .col = 'black', .lwd = 2, .lty = 2){
    for (i in 0:15){
        cov <- (i + 0.5) * .cov
        plot_isoA_line(cov, .L = .L, .ymax = .ymax, .col, .lwd = .lwd, .lty = .lty)
        if (i < 8){
            plot_isoB_line(cov, .ymax, .col, .lwd = .lwd, .lty = .lty)
        }
    }
}

plot_expected_haplotype_structure <- function(.n, .peak_sizes,
                                              .adjust = F, .cex = 1.3, xmax = 0.49){
    .peak_sizes <- .peak_sizes[.peak_sizes[, 'size'] > 0.05, ]
    .peak_sizes[, 'ploidy'] <- nchar(.peak_sizes[, 'structure'])

    decomposed_struct <- strsplit(.peak_sizes[, 'structure'], '')
    .peak_sizes[, 'corrected_minor_variant_cov'] <- sapply(decomposed_struct, function(x){ sum(x == 'B') } ) / .peak_sizes[, 'ploidy']
    .peak_sizes[, 'label'] <- reduce_structure_representation(.peak_sizes[, 'structure'])

    borercases <- .peak_sizes$corrected_minor_variant_cov == 0.5

    for(i in 1:nrow(.peak_sizes)){
        # xmax is in the middle of the last square in the 2d histogram,
        # which is too far from the edge, so I average it with 0.49
        # witch will pull the label bit to the edge
        text( ifelse( borercases[i] & .adjust, (xmax + 0.49) / 2, .peak_sizes$corrected_minor_variant_cov[i]),
             .peak_sizes$ploidy[i] * .n, .peak_sizes[i, 'label'],
             offset = 0, cex = .cex, xpd = T, pos = ifelse( borercases[i] & .adjust, 2, 1))
    }
}

reduce_structure_representation <- function(smudge_labels){
    structures_to_adjust <- (sapply(smudge_labels, nchar) > 4)

    if ( any(structures_to_adjust) ) {
        decomposed_struct <- strsplit(smudge_labels[structures_to_adjust], '')
        As <- sapply(decomposed_struct, function(x){ sum(x == 'A') } )
        Bs <- sapply(decomposed_struct, length) - As
        smudge_labels[structures_to_adjust] <- paste0(As, 'A', Bs, 'B')
    }
    return(smudge_labels)
}

plot_legend <- function(kmer_max, .colour_ramp, .log_scale = T){
    par(mar=c(0,0,2,1))
    plot.new()
    print_title <- ifelse(.log_scale, 'log kmers pairs', 'kmers pairs')
    title(print_title)
    for(i in 1:32){
        rect(0,(i - 0.01) / 33, 0.5, (i + 0.99) / 33, col = .colour_ramp[i])
    }
    # kmer_max <- max(smudge_container$dens)
    if( .log_scale == T ){
        for(i in 0:6){
            text(0.75, i / 6, rounding(10^(log10(kmer_max) * i / 6)), offset = 0)
        }
    } else {
        for(i in 0:6){
            text(0.75, i / 6, rounding(kmer_max * i / 6), offset = 0)
        }
    }
}

rounding <- function(number){
    if(number > 1000){
        round(number / 1000) * 1000
    } else if (number > 100){
        round(number / 100) * 100
    } else {
        round(number / 10) * 10
    }
}

#############
## SETTING ##
#############

parser <- ArgumentParser()
parser$add_argument("--homozygous", action="store_true", default = F,
                    help="Assume no heterozygosity in the genome - plotting a paralog structure; [default FALSE]")
parser$add_argument("-i", "--input", default = "*_smu.txt",
                    help="name of the input tsv file with covarages [default \"*_smu.txt\"]")
parser$add_argument("-s", "--smudges", default = "not_specified",
                    help="name of the input tsv file with annotated smudges and their respective sizes")
parser$add_argument("-o", "--output", default = "smudgeplot",
                    help="name pattern used for the output files (OUTPUT_smudgeplot.png, OUTPUT_summary.txt, OUTPUT_warrnings.txt) [default \"smudgeplot\"]")
parser$add_argument("-t", "--title",
                    help="name printed at the top of the smudgeplot [default none]")
parser$add_argument("-q", "--quantile_filt", type = "double",
                    help="Remove kmer pairs with coverage over the specified quantile; [default none]")
parser$add_argument("-n", "--n_cov", type = "double",
                    help="the haploid coverage of the sequencing data [default inference from data]")
parser$add_argument("-c", "-cov_filter", type = "integer",
                    help="Filter pairs with one of them having coverage bellow specified threshold [default 0]")
parser$add_argument("-ylim", type = "integer", 
                    help="The upper limit for the coverage sum (the y axis)")
parser$add_argument("-col_ramp", default = "viridis",
                    help="A colour ramp available in your R session [viridis]")
parser$add_argument("--invert_cols", action="store_true", default = F,
                    help="Set this flag to invert colorus of Smudgeplot (dark for high, light for low densities)")
 
args <- parser$parse_args()

colour_ramp_log <- get_col_ramp(args, 16) # create palette for the log plots
colour_ramp <- get_col_ramp(args) # create palette for the linear plots

if ( !file.exists(args$input) ) {
    stop("The input file not found. Please use --help to get help", call.=FALSE)
}

cov_tab <- read.table(args$input, header = T) # col.names = c('covB', 'covA', 'freq', 'is_error'), 
smudge_tab <- read.table(args$smudges, col.names = c('structure', 'size'))

# total covarate of the kmer pair
cov_tab[, 'total_pair_cov'] <- cov_tab[, 'covA'] + cov_tab[, 'covB']
# calcualte relative coverage of the minor allele
cov_tab[, 'minor_variant_rel_cov'] <- cov_tab[, 'covB'] / cov_tab[, 'total_pair_cov']

##### coverage filtering

if ( !is.null(args$c) ){
    threshold <- args$c
    low_cov_filt <- cov_tab[, 'covA'] < threshold | cov_tab[, 'covB'] < threshold 
    # smudge_warn(args$output, "Removing", sum(cov_tab[low_cov_filt, 'freq']), 
    #             "kmer pairs for which one of the pair had coverage below",
    #             threshold, paste0("(Specified by argument -c ", args$c, ")"))
    cov_tab <- cov_tab[!low_cov_filt, ]
    # smudge_warn(args$output, "Processing", sum(cov_tab[, 'freq']), "kmer pairs")    
}

##### quantile filtering
if ( !is.null(args$q) ){
    # quantile filtering (remove top q%, it's not really informative)    
    threshold <- wtd.quantile(cov_tab[, 'total_pair_cov'], args$q, cov_tab[, 'freq'])
    high_cov_filt <- cov_tab[, 'total_pair_cov'] < threshold
    # smudge_warn(args$output, "Removing", sum(cov_tab[!high_cov_filt, 'freq']), 
    #             "kmer pairs with coverage higher than",
    #             threshold, paste0("(", args$q, " quantile)"))
    cov_tab <- cov_tab[high_cov_filt, ]
}

cov <- args$n_cov
if (cov == wtd.quantile(cov_tab[, 'total_pair_cov'], 0.95, cov_tab[, 'freq'])){
    ylim <- c(min(cov_tab[, 'total_pair_cov']), max(cov_tab[, 'total_pair_cov']))
} else {
    ylim <- c(min(cov_tab[, 'total_pair_cov']) - 1, # or 0?
          min(max(100, 10*cov), max(cov_tab[, 'total_pair_cov'])))
}

xlim <- c(0, 0.5)
error_fraction <- sum(cov_tab[, 'is_error'] * cov_tab[, 'freq']) / sum(cov_tab[, 'freq']) * 100
error_string <- paste("err =", round(error_fraction, 1), "%")
cov_string <- paste0("1n = ", cov)

if (!is.null(args$ylim)){ # if ylim is specified, set the boundary by the argument instead
    ylim[2] <- args$ylim
}

fig_title <- ifelse(length(args$title) == 0, NA, args$title[1])
# histogram_bins = max(30, args$nbins)

##########
# LINEAR #
##########
pdf(paste0(args$output,'_smudgeplot.pdf'))

# layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
layout(matrix(c(4,2,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_alt(cov_tab, ylim, colour_ramp_log)
if (cov > 0){
    plot_expected_haplotype_structure(cov, smudge_tab, T, xmax = 0.49)
}


# 4 legend
plot_legend(max(cov_tab[, 'freq']), colour_ramp, F)

### add annotation
# print smudge sizes
plot.new()
if (cov > 0){
    legend('topleft', bty = 'n', reduce_structure_representation(smudge_tab[,'structure']), cex = 1.1)
    legend('top', bty = 'n', legend = round(smudge_tab[,2], 2), cex = 1.1)
    legend('bottomleft', bty = 'n', legend = c(cov_string, error_string), cex = 1.1)
} else {
    legend('bottomleft', bty = 'n', legend = error_string, cex = 1.1)
}

plot.new()
mtext(bquote(italic(.(fig_title))), side=3, adj=0.1, line=-3, cex = 1.6)


dev.off()

############
# log plot #
############

pdf(paste0(args$output,'_smudgeplot_log10.pdf'))

layout(matrix(c(4,2,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# cov_tab[, 'freq'] <- log10(cov_tab[, 'freq'])
# 1 smudge plot
plot_alt(cov_tab, ylim, colour_ramp_log, log = T)

if (cov > 0){
    plot_expected_haplotype_structure(cov, smudge_tab, T, xmax = 0.49)
}

# 4 legend
plot_legend(max(cov_tab[, 'freq']), colour_ramp_log, T)

# print smudge sizes
plot.new()
if (cov > 0){
    legend('topleft', bty = 'n', reduce_structure_representation(smudge_tab[,'structure']), cex = 1.1)
    legend('top', bty = 'n', legend = round(smudge_tab[,2], 2), cex = 1.1)
    legend('bottomleft', bty = 'n', legend = c(cov_string, error_string), cex = 1.1)
} else {
    legend('bottomleft', bty = 'n', legend = error_string, cex = 1.1)
}


plot.new()
mtext(bquote(italic(.(fig_title))), side=3, adj=0.1, line=-3, cex = 1.6)

dev.off()
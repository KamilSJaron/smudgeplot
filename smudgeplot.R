#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
# args[1] - input tsv file
# args[2] - output [default: smudgeplot.png]

if(length(args) > 0) {
    if ( args[1] == "--help"){
        stop("Usage: smudgeplot.R [input.tsv] [output.png] [plot_title] [haplod_cov]", call.=FALSE)
    }
}

infile <- ifelse(length(args) < 1, 'coverages_2.tsv', args[1])
outfile <- ifelse(length(args) < 2, 'smudgeplot.png', args[2])
fig_title <- ifelse(length(args) < 3, NA, args[3])
n <- ifelse(length(args) < 4, NA, as.numeric(args[4]))

library(methods)
library(smudgeplot)

# prepare colours
pal <- brewer.pal(11,'Spectral')
rf <- colorRampPalette(rev(pal[3:11]))
colour_ramp <- rf(32)

cov <- read.table(infile)

# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
# total covarate of the kmer pair
total_pair_cov <- cov$V1 + cov$V2

# quantile filtering (remove top 1%, it's not really informative)
# THIS SHOULD JUST RESCALE AXIS, NOT BE KICKED OUT
high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

if( is.na(n) ){
    n <- round(estimate_1n_coverage(),1)
}

# the lims trick will make sure that the last column of squares will have the same width as the other squares
k <- kde2d(minor_variant_rel_cov, total_pair_cov, n=30,
           lims = c(0.02, 0.48, min(total_pair_cov), max(total_pair_cov)))
k_toplot <- k
k_toplot$z <- sqrt(k_toplot$z)

png(outfile)

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(k_toplot, n, colour_ramp)
plot_expected_haplotype_structure(n)
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov)
# 4 legend
plot_legend(k_toplot, total_pair_cov, colour_ramp)

dev.off()
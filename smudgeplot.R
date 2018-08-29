#!/usr/bin/env Rscript

#############
## SETTING ##
#############

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

######################
## INPUT PROCESSING ##
######################

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
    draft_n <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)
} else {
    draft_n <- n
}

ymax <- min(10*draft_n, max(total_pair_cov))
ymin <- min(total_pair_cov) - 1

smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov,
                                         .nbins = 40, .ylim = c(ymin, ymax))

#############
## SUMMARY ##
#############

peak_points <- peak_agregation(smudge_container)
peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)
n <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov)
peak_sizes$structure <- apply(peak_sizes, 1,
                              function(x){ guess_genome_structure(x, n)})
peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)

genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]

##########
## PLOT ##
##########

png(outfile)

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, n, colour_ramp)
plot_expected_haplotype_structure(n, peak_sizes, T)
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov, ymax, fig_title, genome_ploidy)
# 4 legend
plot_legend(smudge_container, total_pair_cov, colour_ramp)

dev.off()
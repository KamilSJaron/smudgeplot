#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("argparse"))
library("smudgeplot")

#############
## SETTING ##
#############

parser <- ArgumentParser()
parser$add_argument("-v", "--version", action="store_true", default = FALSE,
                    help="print the version and exit")
parser$add_argument("-i", "--input", default = "coverages_2.tsv",
                    help="name of the input tsv file with covarages.")
parser$add_argument("-o", "--output", default = "smudeplot",
                    help="name pattern used for the output files.")
parser$add_argument("-t", "--title",
                    help="name printed at the top of the smudgeplot.")
parser$add_argument("-n", "--n_cov", type = "integer",
                    help="the haploid coverage of the sequencing data (inference from data by default - recommended)")
args <- parser$parse_args()

if ( args$version ) {
    stop("Smudgeplot v0.1.0", call.=FALSE)
}

# prepare colours
pal <- brewer.pal(11,'Spectral')
rf <- colorRampPalette(rev(pal[3:11]))
colour_ramp <- rf(32)

######################
## INPUT PROCESSING ##
######################

n <- args$n_cov
cov <- read.table(args$input)

# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
# total covarate of the kmer pair
total_pair_cov <- cov$V1 + cov$V2

# quantile filtering (remove top 1%, it's not really informative)
# THIS SHOULD JUST RESCALE AXIS, NOT BE KICKED OUT
high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

if( length(n) == 0 ){
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

# TODO add warning
to_filter <- peak_sizes$ploidy == 1
smuge_warn("WARNING: XX peaks detected with unexpected coverage.\n")
smuge_warn(paste0(capture.output(peak_sizes[to_filter,]), collapse = "\n"))

peak_sizes <- peak_sizes[!to_filter,]

genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]

# TODO write summary

##########
## PLOT ##
##########

png(paste0(args$output,'_smudgeplot.png'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, n, colour_ramp)
plot_expected_haplotype_structure(n, peak_sizes, T)
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov, ymax, args$title, genome_ploidy)
# 4 legend
plot_legend(smudge_container, total_pair_cov, colour_ramp)

dev.off()
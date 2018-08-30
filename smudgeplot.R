#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("smudgeplot"))

#############
## SETTING ##
#############

parser <- ArgumentParser()
parser$add_argument("-v", "--version", action="store_true", default = FALSE,
                    help="print the version and exit")
parser$add_argument("-i", "--input", default = "coverages_2.tsv",
                    help="name of the input tsv file with covarages [default \"coverages_2.tsv\"]")
parser$add_argument("-o", "--output", default = "smudgeplot",
                    help="name pattern used for the output files (OUTPUT_smudgeplot.png, OUTPUT_summary.txt, OUTPUT_warrnings.txt) [default \"smudgeplot\"]")
parser$add_argument("-t", "--title",
                    help="name printed at the top of the smudgeplot [default none]")
parser$add_argument("-n", "--n_cov", type = "integer",
                    help="the haploid coverage of the sequencing data [default inference from data]")
parser$add_argument("-L", "--low_cutoff", type = "integer",
                    help="the lower boundary used when dumping kmers from jellyfish [default min(total_pair_cov) / 2]")
parser$add_argument("-nbins", type = "integer", default = 40,
                    help="the number of nbins used for smudgeplot matrix (nbins x nbins) [default 40]")
parser$add_argument("-k", "--kmer_size", type = "integer", default = 21,
                    help="The kmer size used to calculate kmer spectra [default 21]")

args <- parser$parse_args()
version_message <- "Smudgeplot v0.1.0"

if ( args$version ) {
    stop(version_message, call.=FALSE)
} else {
    smuge_warn("running", version_message)
}

# prepare colours TODO integrate in package
pal <- brewer.pal(11,'Spectral')
rf <- colorRampPalette(rev(pal[3:11]))
colour_ramp <- rf(32)

smuge_warn("\n######################")
smuge_warn("## INPUT PROCESSING ##")
smuge_warn("######################")

cov <- read.table(args$input)

# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
# total covarate of the kmer pair
total_pair_cov <- cov$V1 + cov$V2

# quantile filtering (remove top 1%, it's not really informative)
high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
smuge_warn("Removing", sum(!high_cov_filt), "kmer pairs with coverage higher than",
           quantile(total_pair_cov, 0.99), "(0.99 quantile)")
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

L <- ifelse( length(args$L) == 0, min(total_pair_cov) / 2, args$L)
n_subset_est <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)

draft_n <- ifelse(length(args$n_cov) == 0, n_subset_est, args$n_cov)

ymax <- min(10*draft_n, max(total_pair_cov))
ymin <- min(total_pair_cov) - 1

smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov,
                                         .nbins = args$nbins, .ylim = c(ymin, ymax))

smuge_warn("\n#############")
smuge_warn("## SUMMARY ##")
smuge_warn("#############")

peak_points <- peak_agregation(smudge_container)
peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)
n_peak_est <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov)

n <- ifelse(length(args$n_cov) == 0, n_peak_est, args$n_cov)

peak_sizes$structure <- apply(peak_sizes, 1,
                              function(x){ guess_genome_structure(x, n)})
peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)

to_filter <- peak_sizes$ploidy == 1
if( any(to_filter) ){
    smuge_warn(paste(sum(to_filter), "peaks of kmer pairs detected with coverage < (1n_coverage * 2) =", round(n * 2, 1)))
    tab_to_print <- peak_sizes[to_filter,c(2,3,8,9)]
    tab_to_print <- round(tab_to_print, 2)
    colnames(tab_to_print) <- c('kmers_in_peak[#]', 'kmers_in_peak[proportion]', 'summit B / (A + B)', 'summit A + B')
    smuge_warn(paste0(capture.output(tab_to_print), collapse = "\n"))
    peak_sizes <- peak_sizes[!to_filter,]
}

peak_sizes$rel_size <- peak_sizes$rel_size / sum(peak_sizes$rel_size)
genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]

smuge_warn("* User defined 1n coverage:", args$n_cov)
smuge_warn("* Subset 1n coverage estimate:", n_subset_est)
smuge_warn("* Highest peak 1n coverage estimate:", round(n_peak_est, 1))
smuge_warn("1n coverage used in smudgeplot (one of the three above):", round(n, 1))

# TODO write summary

for_sure_heterozygous <- sapply(peak_sizes$structure, function(x){sum(unlist(strsplit(x, split = '')) == 'B') == 1})
minimal_number_of_heterozygous_loci <- ceiling(sum(peak_sizes$abs_size[for_sure_heterozygous]) / args$kmer_size)

smuge_summary("* Minimal number of heterozygous loci:", minimal_number_of_heterozygous_loci)
smuge_summary("Note: This number is NOT an estimate of the total number heterozygous loci, it's merly setting the lower boundary if the inference of heterozygosity peaks is correct.")

smuge_summary("* Proportion of heterozygosity carried by pairs in different genome copies (table)")
smuge_summary(paste0(capture.output(
    data.frame(
        genome_copies = 2:max(peak_sizes$ploidy),
        propotion_of_heterozygosity = round(sapply(2:max(peak_sizes$ploidy),
                                 function(x){ sum(peak_sizes$rel_size[peak_sizes$ploidy == x]) }
                                ), 2))), collapse = "\n"))

carried_by_paralogs <- sum(peak_sizes$rel_size[peak_sizes$ploidy > genome_ploidy])
smuge_summary("* Proportion of heterozygosity carried by paralogs:", carried_by_paralogs)

smuge_summary("* Summary of all detected peaks (table)")
tab_to_print <- peak_sizes[,c(11,2,3,8,9)]
tab_to_print[,c(2:5)] <- round(tab_to_print[,c(2:5)], 2)
colnames(tab_to_print) <- c('peak','kmers [#]', 'kmers [proportion]', 'summit B / (A + B)', 'summit A + B')
smuge_summary(paste0(capture.output(tab_to_print), collapse = "\n"))

smuge_warn("\n##########")
smuge_warn("## PLOT ##")
smuge_warn("##########")

fig_title <- ifelse(length(args$title) == 0, NA, args$title[1])

png(paste0(args$output,'_smudgeplot.png'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, n, colour_ramp)
plot_expected_haplotype_structure(n, peak_sizes, T)
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov,
                ymax, fig_title, genome_ploidy, peak_sizes[,c(11,3)])
# 4 legend
plot_legend(smudge_container, total_pair_cov, colour_ramp)

dev.off()
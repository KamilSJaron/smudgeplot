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
parser$add_argument("--homozygous", action="store_true", default = F,
                    help="Assume no heterozygosity in the genome - plotting a paralog structure; [default FALSE]")
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
    smudge_warn(args$output, "running", version_message)
}

# for easier manipulation I store estimated things in a list
smudge_summary <- list()

smudge_warn(args$output, "\n######################")
smudge_warn(args$output, "## INPUT PROCESSING ##")
smudge_warn(args$output, "######################")

cov <- read.table(args$input)

# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
# total covarate of the kmer pair
total_pair_cov <- cov$V1 + cov$V2

# quantile filtering (remove top 1%, it's not really informative)
high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
smudge_warn(args$output, "Removing", sum(!high_cov_filt), "kmer pairs with coverage higher than",
           quantile(total_pair_cov, 0.99), "(0.99 quantile)")
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

L <- ifelse( length(args$L) == 0, min(total_pair_cov) / 2, args$L)
smudge_summary$n_subset_est <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)

draft_n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_subset_est, args$n_cov)

ymax <- min(10*draft_n, max(total_pair_cov))
ymin <- min(total_pair_cov) - 1

smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov,
                                         .nbins = args$nbins, .ylim = c(ymin, ymax))

smudge_warn(args$output, "\n#############")
smudge_warn(args$output, "## SUMMARY ##")
smudge_warn(args$output, "#############")

peak_points <- peak_agregation(smudge_container)
peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)
smudge_summary$n_peak_est <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov)

smudge_summary$n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_peak_est, args$n_cov)

# if the organism is completely homozygous, all the detected kmer pairs are corresponding to paralogs
# therefore ther inference will confuse AB peak to AABB peak etc.
# that is recolvable just telling to guess half of the coverage instead
if( args$homozygous ){
    smudge_summary$n <- smudge_summary$n / 2
    smudge_summary$n_subset_est <- smudge_summary$n_subset_est / 2
}

if( L > (smudge_summary$n / 2) ){
    smudge_warn(args$output, "!! Careful, your coverage filter on the lower end (L = ", L,
                ") is higher than half of the 1n coverage estimate ( 1n / 2 = ", round(smudge_summary$n / 2, 2))
    smudge_warn(args$output, "If the real 1n coverage is half of your estimate you would not picked it up due to the filtering.")
    smudge_warn(args$output, "Consider reruning the analysis with lover L as well (sothing like (1n / 2) - 5 should do the job)")
    smudge_warn(args$output, "Another good way for verificaiton would be to compare it to GenomeScope estimate of haploid coverage")
}

peak_sizes$structure <- apply(peak_sizes, 1,
                              function(x){ guess_genome_structure(x, smudge_summary$n)})
peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)

to_filter <- peak_sizes$ploidy == 1
if( any(to_filter) ){
    smudge_warn(args$output, paste(sum(to_filter), "peaks of kmer pairs detected with coverage < (1n_coverage * 2) =", round(n * 2, 1)))
    tab_to_print <- peak_sizes[to_filter,c(2,3,8,9)]
    tab_to_print <- round(tab_to_print, 2)
    colnames(tab_to_print) <- c('kmers_in_peak[#]', 'kmers_in_peak[proportion]', 'summit B / (A + B)', 'summit A + B')
    smudge_warn(args$output, paste0(capture.output(tab_to_print), collapse = "\n"))
    peak_sizes <- peak_sizes[!to_filter,]
}

peak_sizes$rel_size <- peak_sizes$rel_size / sum(peak_sizes$rel_size)
smudge_summary$peak_sizes <- peak_sizes
smudge_summary$genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]

# this will be probably diploid,
# but theoretically one can imagine a species that si completely homozygous tetraploid
if( args$homozygous ){
    smudge_summary$genome_ploidy <- smudge_summary$genome_ploidy / 2
}

generate_summary(args, smudge_summary)

smudge_warn(args$output, "\n##########")
smudge_warn(args$output, "## PLOT ##")
smudge_warn(args$output, "##########")

fig_title <- ifelse(length(args$title) == 0, NA, args$title[1])
colour_ramp <- get_default_col_ramp() # get the default colour ramp (Spectral, 11)

png(paste0(args$output,'_smudgeplot.png'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, smudge_summary$n, colour_ramp)
plot_expected_haplotype_structure(smudge_summary$n, peak_sizes, T)
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov,
                ymax, fig_title, smudge_summary$genome_ploidy,
                peak_sizes[,c(11,3)], smudge_summary$n)
# 4 legend
plot_legend(smudge_container, total_pair_cov, colour_ramp)

dev.off()
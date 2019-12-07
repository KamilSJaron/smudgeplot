#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("smudgeplot"))

#############
## SETTING ##
#############

parser <- ArgumentParser()
parser$add_argument("--homozygous", action="store_true", default = F,
                    help="Assume no heterozygosity in the genome - plotting a paralog structure; [default FALSE]")
parser$add_argument("-i", "--input", default = "coverages_2.tsv",
                    help="name of the input tsv file with covarages [default \"coverages_2.tsv\"]")
parser$add_argument("-o", "--output", default = "smudgeplot",
                    help="name pattern used for the output files (OUTPUT_smudgeplot.png, OUTPUT_summary.txt, OUTPUT_warrnings.txt) [default \"smudgeplot\"]")
parser$add_argument("-t", "--title",
                    help="name printed at the top of the smudgeplot [default none]")
parser$add_argument("-q", "--quantile_filt", type = "double",
                    help="Remove kmer pairs with coverage over the specified quantile; [default none]")
parser$add_argument("-n", "--n_cov", type = "double",
                    help="the haploid coverage of the sequencing data [default inference from data]")
parser$add_argument("-L", "--low_cutoff", type = "integer",
                    help="the lower boundary used when dumping kmers [default min(total_pair_cov) / 2]")
parser$add_argument("-nbins", type = "integer",
                    help="the number of nbins used for smudgeplot matrix (nbins x nbins) [default autodetection]")
parser$add_argument("-k", "--kmer_size", type = "integer", default = 21,
                    help="The kmer size used to calculate kmer spectra [default 21]")

args <- parser$parse_args()

iterative_nbins <- F
if( is.null(args$nbins) ){
    args$nbins <- 40
    iterative_nbins <- T
}

# for easier manipulation I store estimated things in a list
smudge_summary <- list()

smudge_warn(args$output, "\n######################")
smudge_warn(args$output, "## INPUT PROCESSING ##")
smudge_warn(args$output, "######################")

if ( !file.exists(args$input) ) {
    stop("The input file not found. Please use --help to get help", call.=FALSE)
}

cov <- read.table(args$input)

# calcualte relative coverage of the minor allele
minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
# total covarate of the kmer pair
total_pair_cov <- cov$V1 + cov$V2

if ( !is.null(args$q) ){
    # quantile filtering (remove top 1%, it's not really informative)
    high_cov_filt <- quantile(total_pair_cov, args$q) > total_pair_cov
    smudge_warn(args$output, "Removing", sum(!high_cov_filt), "kmer pairs with coverage higher than",
               quantile(total_pair_cov, args$q), paste0("(", args$q, " quantile)"))
    minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
    total_pair_cov <- total_pair_cov[high_cov_filt]
}

L <- ifelse( length(args$L) == 0, min(total_pair_cov) / 2, args$L)
smudge_summary$n_subset_est <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)

draft_n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_subset_est, args$n_cov)

ymax <- min(10*draft_n, max(total_pair_cov))
ymin <- min(total_pair_cov) - 1

smudge_warn(args$output, "\n#############")
smudge_warn(args$output, "## SUMMARY ##")
smudge_warn(args$output, "#############")

dulpicit_structures <- T
repeat {
    smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov,
                                             .nbins = args$nbins, .ylim = c(ymin, ymax))

    peak_points <- peak_agregation(smudge_container)
    peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)

    the_smallest_n <- min(get_trinoploid_1n_est(peak_sizes), draft_n)
    smudge_summary$n_peak_est <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov, the_smallest_n)

    if (length(args$n_cov) == 0) {
        if( abs(log2(smudge_summary$n_subset_est / smudge_summary$n_peak_est)) > 1 & !args$homozygous){
            smudge_summary$n <- smudge_summary$n_subset_est
        } else {
            smudge_summary$n <- smudge_summary$n_peak_est
        }
    } else {
        smudge_summary$n <- args$n_cov 
    }

    # if the organism is completely homozygous, all the detected kmer pairs are corresponding to paralogs
    # therefore ther inference will confuse AB peak to AABB peak etc.
    # that is recolvable just telling to guess half of the coverage instead
    if( args$homozygous ){
        smudge_summary$n <- smudge_summary$n / 2
        smudge_summary$n_subset_est <- smudge_summary$n_subset_est / 2
    }

    peak_sizes$structure <- apply(peak_sizes, 1,
                                  function(x){ guess_genome_structure(x, smudge_summary$n)})

    dulpicit_structures <- any(table(peak_sizes$structure) > 1)
    # if there are more smudges on the same location & if user have not specified nbins
    if(dulpicit_structures & iterative_nbins){
        if(args$nbins > 20){
            args$nbins <- args$nbins - 5
        } else {
            args$nbins <- args$nbins - 2
        }
        smudge_warn(args$output, "detecting two smudges at the same positions, not enough data for this number of bins lowering number of bins to ", args$nbins)
    } else {
        break
    }
}

if( abs(log2(smudge_summary$n_subset_est / smudge_summary$n_peak_est)) > 1 & !args$homozygous){
    smudge_warn(args$output, "!! Careful, the two types of estimates of 1n coverage differ a lot (",
                smudge_summary$n_subset_est, "and", smudge_summary$n_peak_est, ")")
    smudge_warn(args$output, "meaning that at least of one of the smudgeplot methods to estimate the haploid coverage got it wrong")
    smudge_warn(args$output, "Using subset estimate instead of highest peak estimate (less precise but also less often completely wrong)")
    smudge_warn(args$output, "Does the smudgeplot look sane? Is at least one of the 1n estimates close a GenomeScope estimate?")
    smudge_warn(args$output, "You can help us imrove this software by sharing this strange smudgeplot on https://github.com/KamilSJaron/smudgeplot/issues.")
}

if( L > (smudge_summary$n / 2) & !args$homozygous ){
    smudge_warn(args$output, "!! Careful, your coverage filter on the lower end (L = ", L,
                ") is higher than half of the 1n coverage estimate ( 1n / 2 = ", round(smudge_summary$n / 2, 2))
    smudge_warn(args$output, "If the real 1n coverage is half of your estimate you would not picked it up due to the filtering.")
    smudge_warn(args$output, "If you have sufficient coverage, consider reruning the analysis with lower L (something like (1n / 2) - 5)")
    smudge_warn(args$output, "One good way for verificaiton would be to compare it to GenomeScope estimate of haploid coverage")
}

peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)

to_filter <- peak_sizes$ploidy == 1
if( any(to_filter) ){
    smudge_warn(args$output, paste(sum(to_filter), "peaks of kmer pairs detected with coverage < (1n_coverage * 2) =", round(smudge_summary$n * 2, 1)))
    tab_to_print <- peak_sizes[to_filter,c(2,3,8,9)]
    tab_to_print <- round(tab_to_print, 2)
    colnames(tab_to_print) <- c('kmers_in_peak[#]', 'kmers_in_peak[proportion]', 'summit B / (A + B)', 'summit A + B')
    smudge_warn(args$output, paste0(capture.output(tab_to_print), collapse = "\n"))
    peak_sizes <- peak_sizes[!to_filter,]
}

peak_sizes$rel_size <- peak_sizes$rel_size / sum(peak_sizes$rel_size)
peak_sizes <- peak_sizes[order(peak_sizes$rel_size, decreasing = T),]
smudge_summary$peak_sizes <- peak_sizes

# genome ploidy is the ploidy with highest number of corresponding kmer pairs regardless of topology
considered_ploidies <- unique(peak_sizes$ploidy)
ploidy_with_most_smudges <- which.max(sapply(considered_ploidies, function(x){ sum(peak_sizes[peak_sizes$ploidy == x,'rel_size']) }) )
smudge_summary$genome_ploidy <- considered_ploidies[ploidy_with_most_smudges]
# smudge_summary$genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]

# this will be probably diploid,
# but theoretically one can imagine a species that si completely homozygous tetraploid
if( args$homozygous ){
    smudge_summary$genome_ploidy <- smudge_summary$genome_ploidy / 2
    if(!smudge_summary$genome_ploidy %in% c(2,4,6,8)){
        smudge_warn(args$output, "Guessing really strange ploidy", smudge_summary$genome_ploidy, "perhaps there is not enough coverage for a good inference.")
        smudge_warn(args$output, "You can trust the plot, but guessed ploidy or peak detection might be completely off.")
    }
}

generate_summary(args, smudge_summary)

smudge_warn(args$output, "\n##########")
smudge_warn(args$output, "## PLOT ##")
smudge_warn(args$output, "##########")

fig_title <- ifelse(length(args$title) == 0, NA, args$title[1])
colour_ramp <- get_default_col_ramp() # get the default colour ramp (Spectral, 11)

png(paste0(args$output,'_smudgeplot_log10.png'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, smudge_summary$n, colour_ramp)
plot_expected_haplotype_structure(smudge_summary$n, peak_sizes, T, xmax = max(smudge_container$x))
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov,
                ymax, smudge_summary, args$nbins, fig_title)
# 4 legend
plot_legend(smudge_container, colour_ramp)

dev.off()

# replace the log transformed values by non-transformed
smudge_container$z <- smudge_container$dens

png(paste0(args$output,'_smudgeplot.png'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, smudge_summary$n, colour_ramp)
plot_expected_haplotype_structure(smudge_summary$n, peak_sizes, T, xmax = max(smudge_container$x))
# 2,3 hist
plot_histograms(minor_variant_rel_cov, total_pair_cov,
                ymax, smudge_summary, args$nbins, fig_title)
# 4 legend
plot_legend(smudge_container, colour_ramp, F)

dev.off()
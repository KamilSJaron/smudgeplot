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
parser$add_argument("-i", "--input", default = "*_smu.txt",
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
parser$add_argument("-c", "-cov_filter", type = "integer",
                    help="Filter pairs with one of them having coverage bellow specified threshold [default 0; disables parameter L]")
parser$add_argument("-ylim", type = "integer", 
                    help="The upper limit for the coverage sum (the y axis)")
parser$add_argument("-nbins", type = "integer",
                    help="the number of nbins used for smudgeplot matrix (nbins x nbins) [default autodetection]")
parser$add_argument("-col_ramp", default = "viridis",
                    help="A colour ramp available in your R session [viridis]")
parser$add_argument("--invert_cols", action="store_true", default = F,
                    help="Set this flag to invert colorus of Smudgeplot (dark for high, light for low densities)")
parser$add_argument("--plot_err_line", action="store_true", default = F,
                    help="Set this flag to add a line of theh higher expected occurance of errors paired with genomic k-mers")
parser$add_argument("--just_plot", action="store_true", default = F,
                    help="Turns off the inference of coverage and annotation of smudges; simply generates smudgeplot. (default False)")
parser$add_argument("--alt_plot", action="store_true", default = F,
                    help="Uses a new way to plot smudgeplots using tiling strategy, which is likely to be the default for the Oriel 0.3.0 release (default False)")


args <- parser$parse_args()

colour_ramp_log <- get_col_ramp(args, 16) # create palette for the log plots
colour_ramp <- get_col_ramp(args) # create palette for the linear plots

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

cov_tab <- read.table(args$input, col.names = c('covB', 'covA', 'freq'))

# total covarate of the kmer pair
cov_tab[, 'total_pair_cov'] <- cov_tab[, 'covA'] + cov_tab[, 'covB']
# calcualte relative coverage of the minor allele
cov_tab[, 'minor_variant_rel_cov'] <- cov_tab[, 'covB'] / cov_tab[, 'total_pair_cov']

if ( !is.null(args$q) ){
    # quantile filtering (remove top q%, it's not really informative)
    threshold <- quantile(cov_tab[, 'total_pair_cov'], args$q)
    high_cov_filt <- cov_tab[, 'total_pair_cov'] < threshold
    smudge_warn(args$output, "Removing", sum(cov_tab[!high_cov_filt, 'freq']), 
                "kmer pairs with coverage higher than",
                threshold, paste0("(", args$q, " quantile)"))
    cov_tab <- cov_tab[high_cov_filt, ]
}

##### alt plot
if( args$alt_plot ){
    ylim <- c(0, max(cov_tab[, 'total_pair_cov']))
    if (!is.null(args$ylim)){ # if ylim is specified, set the bounday by the argument instead
        ylim[2] <- args$ylim
    }
    plot_name <- paste0(args$output,'_alt_smudgeplot.pdf')
    pdf(plot_name)
        plot_alt(cov_tab, ylim, colour_ramp)
    dev.off()
    stop(paste("Alt plot has been generated: ", plot_name, " we are done here!"), call.=FALSE)
}

#####

if ( !is.null(args$c) ){
    threshold <- args$c
    low_cov_filt <- cov_tab[, 'covA'] < threshold | cov_tab[, 'covB'] < threshold 
    smudge_warn(args$output, "Removing", sum(cov_tab[low_cov_filt, 'freq']), 
                "kmer pairs for which one of the pair had coverage below",
                threshold, paste0("(Specified by argument -c ", args$c, ")"))
    cov_tab <- cov_tab[!low_cov_filt, ]
    smudge_warn(args$output, "Processing", sum(cov_tab[, 'freq']), "kmer pairs")
    L <- min(cov_tab[, 'covB'])
} else {
    L <- ifelse(length(args$L) == 0, min(cov_tab[, 'covB']), args$L)
}
smudge_summary$n_subset_est <- round(estimate_1n_coverage_1d_subsets(cov_tab), 1)

if (!args$just_plot){
    draft_n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_subset_est, args$n_cov)
    ylim <- c(min(cov_tab[, 'total_pair_cov']) - 1, # or 0?
          min(10*draft_n, max(cov_tab[, 'total_pair_cov'])))
} else {
    ylim <- c(0,max(cov_tab[, 'total_pair_cov'])) # if there is no inference, the default coverage is max cov
}

xlim <- c(0, 0.5)

if (!is.null(args$ylim)){ # if ylim is specified, set the bounday by the argument instead
    ylim[2] <- args$ylim
}


smudge_warn(args$output, "\n#############")
smudge_warn(args$output, "## SUMMARY ##")
smudge_warn(args$output, "#############")

dulpicit_structures <- T
repeat {
    smudge_container <- get_smudge_container(cov_tab, .nbins = args$nbins, .xlim = xlim, .ylim = ylim)

    peak_points <- peak_agregation(smudge_container)
    if (args$just_plot){
        break
    }
    peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)

    the_smallest_n <- min(get_trinoploid_1n_est(peak_sizes), draft_n)
    smudge_summary$n_peak_est <- estimate_1n_coverage_highest_peak(peak_sizes, cov_tab, the_smallest_n)
    if (length(args$n_cov) == 0) {
        if( (abs(log2(smudge_summary$n_subset_est / smudge_summary$n_peak_est)) > .95 | smudge_summary$n_peak_est < smudge_summary$n_subset_est) & !args$homozygous){
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

if (!args$just_plot){
    # checks / warnings for potential inference problems
    if( abs(log2(smudge_summary$n_subset_est / smudge_summary$n_peak_est)) > 1 & !args$homozygous){
        smudge_warn(args$output, "!! Warning, the two types of estimates of 1n coverage differ a lot (",
                    smudge_summary$n_subset_est, "and", smudge_summary$n_peak_est, ")")
        smudge_warn(args$output, "i.e. at least of one of the smudgeplot methods to estimate the haploid coverage got it wrong")
        smudge_warn(args$output, "Using subset estimate instead of highest peak estimate (less precise but also less often completely wrong)")
        smudge_warn(args$output, "Does the smudgeplot look sane? Is at least one of the 1n estimates close a GenomeScope estimate?")
        smudge_warn(args$output, "You can help us imrove this software by sharing this strange smudgeplot on https://github.com/KamilSJaron/smudgeplot/issues.")
    }
    if( L > (smudge_summary$n / 2) & !args$homozygous ){
        smudge_warn(args$output, "!! Warning, your coverage filter on the lower end (L = ", L,
                    ") is higher than half of the 1n coverage estimate ( 1n / 2 = ", round(smudge_summary$n / 2, 2))
        smudge_warn(args$output, "If the real 1n coverage is half of your estimate you would not picked it up due to the filtering.")
        smudge_warn(args$output, "If you have sufficient coverage, consider reruning the analysis with lower L (something like (1n / 2) - 5)")
        smudge_warn(args$output, "One good way for verificaiton would be to compare it to GenomeScope estimate of haploid coverage")
    }

    peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
    peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)

    to_filter <- peak_sizes$ploidy <= 1
    if( any(to_filter) ){
        smudge_warn(args$output, paste(sum(to_filter), "peaks of kmer pairs detected with coverage < (1n_coverage * 2) =", round(smudge_summary$n * 2, 1)))
        tab_to_print <- peak_sizes[to_filter,c(2,3,8,9)]
        tab_to_print <- round(tab_to_print, 2)
        colnames(tab_to_print) <- c('kmers_in_peak[#]', 'kmers_in_peak[proportion]', 'summit B / (A + B)', 'summit A + B')
        smudge_warn(args$output, paste0(capture.output(tab_to_print), collapse = "\n"))
        peak_sizes <- peak_sizes[!to_filter,]
        smudge_warn(args$output, "This might be due to kmers with sequencing errors present in the kmer dump.")
        smudge_warn(args$output, "Increasing L might help remove erroneous kmers.")
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
} else {
    smudge_summary$n <- 0 # this is interpreted as "no infered"
}

smudge_warn(args$output, "\n##########")
smudge_warn(args$output, "## PLOT ##")
smudge_warn(args$output, "##########")

fig_title <- ifelse(length(args$title) == 0, NA, args$title[1])

pdf(paste0(args$output,'_smudgeplot_log10.pdf'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, smudge_summary$n, colour_ramp_log)
if (!args$just_plot){
plot_expected_haplotype_structure(smudge_summary$n, peak_sizes, T, xmax = max(smudge_container$x))
}
if (args$plot_err_line){
    plot_seq_error_line(cov_tab)
}

histogram_bins = max(30, args$nbins)
# 2,3 hist
plot_histograms(cov_tab, smudge_summary, fig_title, .ylim = ylim, .bins = histogram_bins) # I am testing here setting the number of bars to the same number as the number of squares
# 4 legend
plot_legend(smudge_container, colour_ramp_log)

dev.off()

# replace the log transformed values by non-transformed
smudge_container$z <- smudge_container$dens

pdf(paste0(args$output,'_smudgeplot.pdf'))

layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# 1 smudge plot
plot_smudgeplot(smudge_container, smudge_summary$n, colour_ramp)
if (!args$just_plot){
    plot_expected_haplotype_structure(smudge_summary$n, peak_sizes, T, xmax = max(smudge_container$x))
}
# plot error line L - 1 / cov ~ cov
if (args$plot_err_line){
    plot_seq_error_line(cov_tab)
}

# 2,3 hist
plot_histograms(cov_tab, smudge_summary, fig_title,
                .ylim = ylim, .bins = histogram_bins)

# 4 legend
plot_legend(smudge_container, colour_ramp, F)

dev.off()

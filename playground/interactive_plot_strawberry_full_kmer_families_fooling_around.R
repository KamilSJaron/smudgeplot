library("methods")
library("argparse")
library("smudgeplot")
library("hexbin")

# preprocessing
# to get simply number of memebers / family (exploration)
# cat data/strawberry_iinumae/kmer_counts_L109_U.tsv  | cut -f 1 > data/strawberry_iinumae/kmer_counts_L109_U_family_members.tsv
# awk '{row_sum = 0; row_max = 0; row_min = 10000; for (i=2; i <= NF; i++){ row_sum += $i; if ($i > row_max){row_max = $i} if ($i < row_min){row_min = $i} } print row_sum "\t" row_min "\t" row_max }' data/strawberry_iinumae/kmer_counts_L109_U.tsv > data/strawberry_iinumae/kmer_counts_L109_U_sums_min_max.tsv
# (exploration)
#
#

args <- ArgumentParser()$parse_args()
args$homozygous <- F
args$input <- './data/strawberry_iinumae/kmer_counts_L109_U.tsv'
args$output = './data/strawberry_iinumae/straw'
args$title = 'F. iinumae'
args$nbins <- 40
args$L <- NULL
args$n_cov <- NULL
args$k <- 21

family_members <- table(read.table('data/strawberry_iinumae/kmer_counts_L109_U_family_members.tsv')$V1)
family_members / sum(family_members)

#          2          3          4          5          6          7          8 
# 0.60463486 0.18252481 0.08881352 0.05260407 0.03399928 0.02262155 0.01480190 

tranf_cov <- read.table('data/strawberry_iinumae/kmer_counts_L109_U_sums_min_max.tsv', col.names = c('sum', 'min', 'max'))

iterative_nbins <- T

smudge_summary <- list()

pairs <- c(tranf_cov$min + tranf_cov$max == tranf_cov$sum)
cov <- tranf_cov[pairs,c(2,3)]
minor_variant_rel_cov <- cov$min / (cov$min + cov$max)
total_pair_cov <- cov$min + cov$max

plot_smudge <- function(minor_variant_rel_cov, total_pair_cov, add_to_title = ''){
	L <- ifelse( length(args$L) == 0, min(total_pair_cov) / 2, args$L)
	smudge_summary$n_subset_est <- 144
	draft_n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_subset_est, args$n_cov)

	ymax <- min(10*draft_n, max(total_pair_cov))
	ymin <- min(total_pair_cov) - 1

	dulpicit_structures <- T
	repeat {
	    smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov, .nbins = args$nbins, .ylim = c(ymin, ymax))
	    peak_points <- peak_agregation(smudge_container)
	    peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)

	    the_smallest_n <- min(get_trinoploid_1n_est(peak_sizes), draft_n)
	    smudge_summary$n_peak_est <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov, the_smallest_n)

	    smudge_summary$n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_peak_est, args$n_cov)

	    peak_sizes$structure <- apply(peak_sizes, 1, function(x){ guess_genome_structure(x, smudge_summary$n)})

	    dulpicit_structures <- any(table(peak_sizes$structure) > 1)
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

	peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
	peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)


	peak_sizes$rel_size <- peak_sizes$rel_size / sum(peak_sizes$rel_size)
	peak_sizes <- peak_sizes[order(peak_sizes$rel_size, decreasing = T),]
	smudge_summary$peak_sizes <- peak_sizes

	considered_ploidies <- unique(peak_sizes$ploidy)
	ploidy_with_most_smudges <- which.max(sapply(considered_ploidies, function(x){ sum(peak_sizes[peak_sizes$ploidy == x,'rel_size']) }) )
	smudge_summary$genome_ploidy <- considered_ploidies[ploidy_with_most_smudges]

	generate_summary(args, smudge_summary)
	fig_title <- paste(args$title, add_to_title)

	layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
	# 1 smudge plot
	plot_smudgeplot(smudge_container, smudge_summary$n, colour_ramp)
	plot_expected_haplotype_structure(smudge_summary$n, peak_sizes, T, xmax = max(smudge_container$x))
	# 2,3 hist
	plot_histograms(minor_variant_rel_cov, total_pair_cov,
	                ymax, smudge_summary, args$nbins, fig_title)
	# 4 legend
	plot_legend(smudge_container, colour_ramp)
}

plot_smudge(minor_variant_rel_cov, total_pair_cov, 'just pairs')
## OK, the C++ software kind of works,
# but the smudgeplot of the families of two DOES NOT look exactly the same as when produced by python kmer pairs script
# it's hard to say from here where is the difference comes from, but... it remains true that the smudgeplots are different


# Now, let look at the non diploid ones
cov <- tranf_cov[!pairs,]
minor_variant_rel_cov <- cov$min / cov$sum
total_pair_cov <- cov$sum

pdf('figures/multi_family_covsum_hist.pdf')
	plot_hist <- hist(total_pair_cov, breaks = 2000, xlim = c(0,5000))
	for(i in c(4, 6,8,10, 16,32) * 144){
		lines(c(i,i), c(0, 100000), col = 'red')
		text(i, max(plot_hist$counts), labels = i / 144, pos = 4)
	}
dev.off()

pdf('figures/multi_family_smudgeplot.pdf')
	plot_smudge(minor_variant_rel_cov, total_pair_cov, 'only_multi_families')
dev.off()

minor_variant_rel_cov <- cov$max / cov$sum
total_pair_cov <- cov$sum

low_cov_data <- c(total_pair_cov < 1300)
minor_variant_rel_cov <- minor_variant_rel_cov[low_cov_data]
total_pair_cov <- total_pair_cov[low_cov_data]

# Create hexbin object and plot
h <- hexbin(data.frame(max_rel_cov = minor_variant_rel_cov,total_pair_cov))

pdf('figures/multi_family_max_rel_cov.pdf')
	plot(h, colramp=viridis)
dev.off()
#plot_smudge(minor_variant_rel_cov, total_pair_cov, 'only_multi_families')

tranf_cov$memebers <- read.table('data/strawberry_iinumae/kmer_counts_L109_U_family_members.tsv')$V1
cov <- tranf_cov[tranf_cov$memebers == 3,]

cov$mid <- cov$sum - (cov$min + cov$max)

cov$min <- round(cov$min / 144)
cov$mid <- round(cov$mid / 144)
cov$max <- round(cov$max / 144)

tri_family_tab <- table(paste(cov$max, cov$mid, cov$min))

tri_represented <- tri_family_tab[tri_family_tab > 5000]
tri_family_overview <- data.frame(name = names(tri_represented), count = as.numeric(tri_represented), stringsAsFactors = F)

tri_family_overview$As <- sapply(strsplit(tri_family_overview$name, ' '), function(x)( as.numeric(x[1]) ) )
tri_family_overview$Bs <- sapply(strsplit(tri_family_overview$name, ' '), function(x)( as.numeric(x[2]) ) )
tri_family_overview$Cs <- sapply(strsplit(tri_family_overview$name, ' '), function(x)( as.numeric(x[3]) ) )

tri_family_overview$ploidy <- rowSums(tri_family_overview[,c('As', 'Bs', 'Cs')])
tri_family_overview$smudge <- paste0(sapply(tri_family_overview$As, function(x) { paste0(rep('A', x), collapse = '') }),
									 sapply(tri_family_overview$Bs, function(x) { paste0(rep('B', x), collapse = '') }),
									 sapply(tri_family_overview$Cs, function(x) { paste0(rep('C', x), collapse = '') }))

tri_family_overview <- tri_family_overview[order(tri_family_overview$ploidy),]
tri_barplot <- tri_family_overview$count
tri_spaces <- c(0.5, ifelse(tri_family_overview$ploidy[-1] == tri_family_overview$ploidy[-length(tri_barplot)], 0.2, 0.5))
names(tri_barplot) <- tri_family_overview$ploidy
annot <- tri_family_overview$count > 5e4

pdf('figures/triallelic_smudge_sizes.pdf')
	bar_pos <- barplot(tri_barplot, space = tri_spaces)
	text(bar_pos[annot], tri_family_overview$count[annot], tri_family_overview$smudge[annot])
dev.off()
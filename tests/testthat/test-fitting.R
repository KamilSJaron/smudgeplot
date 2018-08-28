context("integration")

cov <- read.table('fmlo_smaple_cov_file.tsv')

pal <- brewer.pal(11,'Spectral')
rf <- colorRampPalette(rev(pal[3:11]))
colour_ramp <- rf(32)

minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
total_pair_cov <- cov$V1 + cov$V2

high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

draft_n <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)
ymax <- min(10 * draft_n, max(total_pair_cov))
ymin <- min(total_pair_cov) - 1

nbins <- 15
smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov, .nbins = nbins)

test_that("2d fitting somehow works", {
        expect_gt(sum(smudge_container$dens), 0)
    }
)

# peak_points <- peak_agregation(smudge_container)
# peak_sizes <- get_peak_summary(peak_points, smudge_container)
# n <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov)
# peak_sizes$structure <- apply(peak_sizes, 1,
#                               function(x){ guess_genome_structure(x, n)})
# peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
# peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)
#
# genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]
#
# test_that("M. flo is expected to be estimated as triploid", {
#         expect_equal(genome_ploidy, 3)
#     }
# )

# layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
# plot_smudgeplot(smudge_container, n, colour_ramp)
# plot_expected_haplotype_structure(n, peak_sizes, T)
# plot_histograms(minor_variant_rel_cov, total_pair_cov, ymax, "TEST", genome_ploidy)
# plot_legend(smudge_container, total_pair_cov, colour_ramp)

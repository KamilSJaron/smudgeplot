context("generation of peak summary")

cov <- read.table('fmlo_smaple_cov_file.tsv')

minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
total_pair_cov <- cov$V1 + cov$V2

high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

draft_n <- 205
ymax <- min(10 * draft_n, max(total_pair_cov))
ymin <- min(total_pair_cov) - 1

nbins <- 30
smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov, .nbins = nbins)

peak_points <- peak_agregation(smudge_container)
# peak_points[peak_points$summit == T,]

test_that("peak agregation creates more than one smudge", {
        expect_true(sum(peak_points$summit) > 1)
        expect_true(any(!is.na(peak_points$peak)))
    }
)

peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.05)
the_smallest_n <- get_trinoploid_1n_est(peak_sizes)

n <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov, the_smallest_n)

test_that("test that polished estimate of coverage for mflo makes sense", {
        expect_true(abs(n - 207) < 10)
    }
)

peak_sizes$structure <- apply(peak_sizes, 1,
                              function(x){ guess_genome_structure(x, n)})

peak_sizes$corrected_minor_variant_cov <- sapply(peak_sizes$structure, function(x){round(mean(unlist(strsplit(x, split = '')) == 'B'), 2)})
peak_sizes$ploidy <- sapply(peak_sizes$structure, nchar)

genome_ploidy <- peak_sizes$ploidy[which.max(peak_sizes$rel_size)]

test_that("mflo was estimated triploid", {
        expect_true(genome_ploidy == 3)
    }
)

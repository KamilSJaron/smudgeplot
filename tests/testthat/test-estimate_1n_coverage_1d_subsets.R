context("1d subsets fitting")

cov <- read.table('fmlo_smaple_cov_file.tsv')

minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
total_pair_cov <- cov$V1 + cov$V2

high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

draft_n <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)
expected_n <- 207 # from genomescope

test_that("approximate coverage estimation", {
        expect_true(abs(draft_n - expected_n) < 10)
    }
)

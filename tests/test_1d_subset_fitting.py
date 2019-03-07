import argparse
import unittest
from smudgeplot import smudgedata

import os

class TestFitting(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.infile = open('tests/data/mflo_smaple_cov_file.tsv')
        args.L = 0
        args.nbins = 0
        args.q = 0.95
        self.smudge = smudgedata.smudgedata(args)
        self.smudge.loadData()

    def testLoadData(self):
        self.assertEqual(len(self.smudge.sum_cov), 1000)

    def testInitialNEstimate(self):
        self.smudge.initialNEstimate()
        self.assertAlmostEqual(self.smudge.n_init, 206.76, 1)
#
# minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
# total_pair_cov <- cov$V1 + cov$V2
#
# high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
# minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
# total_pair_cov <- total_pair_cov[high_cov_filt]
#
# draft_n <- round(estimate_1n_coverage_1d_subsets(total_pair_cov, minor_variant_rel_cov), 1)
# expected_n <- 207 # from genomescope
#
# test_that("approximate coverage estimation", {
#         expect_true(abs(draft_n - expected_n) < 10)
#     }
# )
#

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
        args.n = 0
        self.smudge = smudgedata.smudgedata(args)
        self.smudge.loadData()

    def testLoadData(self):
        self.assertEqual(len(self.smudge.sum_cov), 1000)

    def testInitialNEstimate(self):
        self.smudge.initialNEstimate()
        self.assertAlmostEqual(self.smudge.n_init, 206.76, 1)

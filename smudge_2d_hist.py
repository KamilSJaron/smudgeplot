#!/usr/bin/env python3

import argparse
import sys

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='name of the input tsv file with covarages [default \"coverages_2.tsv\"]"')
  parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='smudgeplot')
  # parser.add_argument('-k', help='The length of the kmer.', default=21)
  # parser.add_argument('-v', '--version', action="store_true", default = FALSE, help="print the version and exit")
  args = parser.parse_args()
  kmer_pairs_file = args.infile

  minor_variant_rel_cov = []
  total_pair_cov = []

  for line in kmer_pairs_file:
    c1, c2 = line.split()
    c1 = int(c1)
    c2 = int(c2)
    sum_cov = c1 + c2
    minor_variant_rel_cov.append(c1 / sum_cov)
    total_pair_cov.append(sum_cov)

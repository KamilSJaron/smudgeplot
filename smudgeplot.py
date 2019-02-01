#!/usr/bin/env python3

import argparse
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
from smudgedata import smudgedata

def main():
  parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='name of the input tsv file with covarages [default \"coverages_2.tsv\"]"')
  parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='smudgeplot')
  # parser.add_argument('-k', help='The length of the kmer.', default=21)
  # parser.add_argument('-v', '--version', action="store_true", default = FALSE, help="print the version and exit")
  # # if not specified by user
  # nbins = 40
  # #    L <- ifelse( length(args$L) == 0, min(total_pair_cov) / 2, args$L)
  # # this is an integer division but that should be alright
  # L = int(min(total_pair_cov) / 2)

  logging.basicConfig(level=logging.INFO)
  # parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
  # args = parser.parse_args()
  # args.infile = open('./data/samples/Mflo2_sample_50000_coverages_2.tsv')
  # smudge = smudgedata(args)
  smudge = smudgedata(parser.parse_args())

  logging.info('loading')
  smudge.loadData()
  logging.info('done')

  # quantile filtering
  # -> user specified of qunatiles to be filtered out by the sum of coverages
  # -> turned off by default
  # smudge.quantile_filt()

  # LOG some stats
  # logging.info('Running smudgeplot')
  # logging.info('Processing ' + str(len(total_pair_cov)) + ' kmer pairs.')
  # logging.info('Parameters:')
  # logging.info('L=' + str(L))
  # logging.info('k=21')

  # smudge.est_n_init()
  # ymin = min(total_pair_cov) - 1
  # ymax = min([10*smudge.n_init, max(total_pair_cov)])

  # LOG
  # logging.info('Smudgeplot coverage range ' + str(ymin) + ' - ' + str(ymax))

  logging.info('calculating hist')
  smudge.calculateHist(40,200,1000)
  logging.info('done')
  # LOG some stats
  # while True:
  #   logging.info('running with nbins: ' + str(nbins))
  #   smudges = smudge_matrix(minor_variant_rel_cov, total_pair_cov, nbins, ymin, ymax)
  #   if there are more smudges on the same location & if user have not specified nbins
  #   if smudges.hasDuplicitSmudges() and nbins > 10:
  #     if nbins > 20 :
  #       nbins = nbins - 5
  #     else :
  #       nbins = nbins - 2
  #     logging.warn(args$output, "detecting two smudges at the same positions, not enough data for this number of bins lowering number of bins to ", args$nbins)
  #   else :
  #     break

  logging.info('plotting *_smudgeplot_pythonic.png')
  smudge.plot()
  logging.info('done')

  logging.info('saving *_smudgematrix.tsv')
  smudge.saveMatrix()
  logging.info('done')

    # peak agregation : peak_agregation(smudge_container)
    # peak summary : get_peak_summary(peak_points, smudge_container, 0.02)
    # the smallest n : min(get_trinoploid_1n_est(peak_sizes), draft_n)
    # n peak est : estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov, the_smallest_n)
    # final n est : ifelse(length(args$n_cov) == 0, smudge_summary$n_peak_est, args$n_cov)
    # # if the organism is completely homozygous, all the detected kmer pairs are corresponding to paralogs
    # # therefore ther inference will confuse AB peak to AABB peak etc.
    # # that is recolvable just telling to guess half of the coverage instead
    # # this could be handeled internally, smudge contains all the arguments
    # if args$homozygous:
    #   smudge_summary$n <- smudge_summary$n / 2
    #   smudge_summary$n_subset_est <- smudge_summary$n_subset_est / 2
    # guess ploidy : apply(peak_sizes, 1, function(x){ guess_genome_structure(x, smudge_summary$n)})

if __name__=='__main__':
  main()

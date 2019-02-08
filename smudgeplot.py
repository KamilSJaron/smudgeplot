#!/usr/bin/env python3

import argparse
import sys
import logging
from smudgedata import smudgedata

# Development:
#  - logging is done on the standard error stream using logging package (note different levels info, warning...)
#  - majority opetations are performed though class smudgedata
#  - methods are camel case; variables snake case

################
###  SCRIPT  ###
################

def main():
  parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='name of the input tsv file with covarages [default \"coverages_2.tsv\"]"')
  parser.add_argument('-o', help='The pattern used to name the output [default \"smudgeplot\"].', default='smudgeplot')
  parser.add_argument('-L', help='The lower boundary used when dumping kmers [default min(total_pair_cov) / 2]', type=int, default=0)
  parser.add_argument('-nbins', help='The number of nbins used for smudgeplot matrix (nbins x nbins) [default autodetection]', type=int, default=0)
  # parser.add_argument('-k', help='The length of the kmer.', default=21)
  # parser.add_argument('-v', '--version', action="store_true", default = FALSE, help="print the version and exit")
  # # if not specified by user
  # nbins = 40

  logging.basicConfig(level=logging.INFO)
  smudge = smudgedata(parser.parse_args())

  logging.info('loading kmers')
  smudge.loadData()
  logging.info('done')

  logging.info('L = ' + str(smudge.args.L))
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

  smudge.initialNEstimate()
  ymin = min(smudge.sum_cov) - 1
  ymax = int(min([10*smudge.n_init, max(smudge.sum_cov)]))

  # LOG
  logging.info('smudgeplot coverage range ' + str(ymin) + ' - ' + str(ymax))
  logging.info('initial 1n coverage estimate {:0.2f}'.format(smudge.n_init))

  # LOG some stats
  while True:
    logging.info('calculating hist with nbins: ' + str(smudge.nbins))
    smudge.calculateHist(ymin, ymax)
    smudge.agregateSmudges()
    smudge.countSmudgeSize(treshold = 0.02)
    smudge.brightestSmudgeNEstimate()
    # if the organism is completely homozygous, all the detected kmer pairs are corresponding to paralogs
    # therefore ther inference will confuse AB peak to AABB peak etc.
    # that is recolvable just telling to guess half of the coverage instead
    # if( args$homozygous ){
    #     smudge_summary$n <- smudge_summary$n / 2
    #     smudge_summary$n_subset_est <- smudge_summary$n_subset_est / 2
    # }
    # smudge.guessGenomeStructure()
    # peak_sizes$structure <- apply(peak_sizes, 1,
    #                               function(x){ guess_genome_structure(x, smudge_summary$n)})

    # smudge.nbins != smudge.args.nbins -> if user have set up number of bins, don't repeat, just plot what user have told us to
    if smudge.hasDuplicitSmudges() and smudge.nbins > 10 and smudge.nbins != smudge.args.nbins:
      smudge.lowerNbins()
      logging.warning("detecting two smudges at the same positions, not enough data for this number of bins lowering number of bins to " + str(smudge.nbins))
    else :
      break
  logging.info('done')

  logging.info("saving " + smudge.args.o + "_smudgeplot_pythonic.png")
  smudge.plot()
  logging.info('done')

  logging.info("saving " + smudge.args.o + "_smudgematrix.tsv")
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

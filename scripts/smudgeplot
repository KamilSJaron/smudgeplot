#!/usr/bin/env python3

import argparse
import sys
import logging
from smudgedata import smudgedata

# number name(-development); meaning with more commits over the mentioned version
version = '0.1.3 beta3-development'

# Development:
#  - logging is done on the standard error stream using logging package (note different levels info, warning...)
#  - majority opetations are performed though class smudgedata
#  - methods are camel case; variables snake case

################
###  SCRIPT  ###
################

def main():
  parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
  # parser.add_argument('-v', '--version', action="store_true", default = FALSE, help="print the version and exit")
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='name of the input tsv file with covarages [default \"coverages_2.tsv\"]"')
  parser.add_argument('-o', help='The pattern used to name the output [default \"smudgeplot\"].', default='smudgeplot')
  parser.add_argument('-q', help='Remove kmer pairs with coverage over the specified quantile; [default none]', default=1)
  parser.add_argument('-L', help='The lower boundary used when dumping kmers [default min(total_pair_cov) / 2]', type=int, default=0)
  parser.add_argument('-nbins', help='The number of nbins used for smudgeplot matrix (nbins x nbins) [default autodetection]', type=int, default=0)
  # parser.add_argument('-k', help='The length of the kmer.', default=21)
  parser.add_argument('--homozygous', action="store_true", default = False, help="Assume no heterozygosity in the genome - plotting a paralog structure; [default False]")

  logging.basicConfig(level=logging.INFO)
  smudge = smudgedata(parser.parse_args())

  logging.info('Running smudgeplot v' + version)
  # LOG arguments by user

  logging.info('loading kmers')
  smudge.loadData()
  smudge.quantile_filt()
  logging.info('done')

  # LOG some stats

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
    if smudge.args.homozygous:
      smudge.brightest_smudge_n = smudge.brightest_smudge_n / 2
    smudge.guessGenomeStructure()

    # smudge.nbins != smudge.args.nbins -> if user have set up number of bins, don't repeat, just plot what user have told us to
    if smudge.hasDuplicitSmudges() and smudge.nbins > 10 and smudge.nbins != smudge.args.nbins:
      smudge.lowerNbins()
      logging.warning("detecting two smudges at the same positions, not enough data for this number of bins lowering number of bins to " + str(smudge.nbins))
    else :
      break
  logging.info('final 1n coverage estimate {:0.2f}'.format(smudge.brightest_smudge_n))
  logging.info('done')

  logging.info("saving " + smudge.args.o + "_smudgeplot_pythonic.png")
  smudge.plot((ymin,ymax))
  logging.info('done')

  logging.info("saving " + smudge.args.o + "_smudgematrix.tsv")
  smudge.saveMatrix()
  logging.info('done')

if __name__=='__main__':
  main()

#!/usr/bin/env python3

import argparse
import sys
import logging

class summary:
  def __init__(self, n_init):
    self.n_init = n_init;
  # def report(self):
    # report summary on the standard output

class smudge_matrix:
  def __init__(self, rel_cov, pair_cov, nbins, ymin, ymax):
    self.x = [0.1,0.2,0.3,0.4] # seq(.xlim[1], ((.nbins - 1) / .nbins) * .xlim[2], length = .nbins)
    self.y = [100,200,300,400] # c(seq(.ylim[1], ((.nbins - 1) / .nbins) * .ylim[2], length = .nbins), .ylim[2])
    # freq <-  as.data.frame(table(findInterval(.minor_variant_rel_cov, smudge_container$x, left.open = T),
    #                              findInterval(.total_pair_cov, smudge_container$y, left.open = T)),
    #                        stringsAsFactors = F)
    # freq[,1] <- as.numeric(freq[,1])
    # freq[,2] <- as.numeric(freq[,2])
    # freq <- freq[!freq$Var2 == .nbins+1,]
    #
    # smudge_container$dens <- matrix(0, .nbins, .nbins)
    # smudge_container$dens[cbind(freq[,1], freq[,2])] <- freq[,3]
    #
    # smudge_container$y <- smudge_container$y[1:.nbins]
    # smudge_container$z <- log10(smudge_container$dens)
    # smudge_container$z[is.infinite(smudge_container$z)] <- 0
  def hasDuplicitSmudges(self):
    return(True)

def main():
  parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='name of the input tsv file with covarages [default \"coverages_2.tsv\"]"')
  parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='smudgeplot')
  # parser.add_argument('-k', help='The length of the kmer.', default=21)
  # parser.add_argument('-v', '--version', action="store_true", default = FALSE, help="print the version and exit")
  args = parser.parse_args()
  kmer_pairs_file = args.infile
  logging.basicConfig(level=logging.INFO)

  minor_variant_rel_cov = []
  total_pair_cov = []

  for line in kmer_pairs_file:
    c1, c2 = line.split()
    c1 = int(c1)
    c2 = int(c2)
    sum_cov = c1 + c2
    minor_variant_rel_cov.append(c1 / sum_cov)
    total_pair_cov.append(sum_cov)

  # quantile filtering
  # -> user specified of qunatiles to be filtered out by the sum of coverages
  # -> turned off by default

  # if not specified by user
  nbins = 40
  #    L <- ifelse( length(args$L) == 0, min(total_pair_cov) / 2, args$L)
  # this is an integer division but that should be alright
  L = int(min(total_pair_cov) / 2)
  # to be replaced with actual estimate from the data
  # summary is a class gathering summary info about smudges
  smudge_summary = summary(130)

  logging.info('Running smudgeplot')
  logging.info('Processing ' + str(len(total_pair_cov)) + ' kmer pairs.')
  logging.info('Parameters:')
  logging.info('L=' + str(L))
  logging.info('k=21')

  ymax = min([10*smudge_summary.n_init, max(total_pair_cov)])
  ymin = min(total_pair_cov) - 1
  logging.info('Smudgeplot coverage range ' + str(ymin) + ' - ' + str(ymax))

  while True:
    logging.info('running with nbins: ' + str(nbins))
    smudges = smudge_matrix(minor_variant_rel_cov, total_pair_cov, nbins, ymin, ymax)
    if smudges.hasDuplicitSmudges() and nbins > 10:
      if nbins > 20 :
        nbins = nbins - 5
      else :
        nbins = nbins - 2
    else :
      break

    # smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov,
    #                                          .nbins = args$nbins, .ylim = c(ymin, ymax))
    #
    # peak_points <- peak_agregation(smudge_container)
    # peak_sizes <- get_peak_summary(peak_points, smudge_container, 0.02)
    #
    # the_smallest_n <- min(get_trinoploid_1n_est(peak_sizes), draft_n)
    # smudge_summary$n_peak_est <- estimate_1n_coverage_highest_peak(peak_sizes, minor_variant_rel_cov, total_pair_cov, the_smallest_n)
    #
    # smudge_summary$n <- ifelse(length(args$n_cov) == 0, smudge_summary$n_peak_est, args$n_cov)
    #
    # # if the organism is completely homozygous, all the detected kmer pairs are corresponding to paralogs
    # # therefore ther inference will confuse AB peak to AABB peak etc.
    # # that is recolvable just telling to guess half of the coverage instead
    # if( args$homozygous ){
    #     smudge_summary$n <- smudge_summary$n / 2
    #     smudge_summary$n_subset_est <- smudge_summary$n_subset_est / 2
    # }
    #
    # peak_sizes$structure <- apply(peak_sizes, 1,
    #                               function(x){ guess_genome_structure(x, smudge_summary$n)})
    #
    # dulpicit_structures <- any(table(peak_sizes$structure) > 1)
    # # if there are more smudges on the same location & if user have not specified nbins
    # if(dulpicit_structures & iterative_nbins){
    #     if(args$nbins > 20){
    #         args$nbins <- args$nbins - 5
    #     } else {
    #         args$nbins <- args$nbins - 2
    #     }
    #     smudge_warn(args$output, "detecting two smudges at the same positions, not enough data for this number of bins lowering number of bins to ", args$nbins)
    # } else {
    #     break
    # }


if __name__=='__main__':
  main()

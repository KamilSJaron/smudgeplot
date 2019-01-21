#!/usr/bin/env python3

from collections import defaultdict
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
# import multiprocessing as mp
# import numpy as np
import argparse
import itertools
import operator
import sys

##############
### SCRIPT ###
##############

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Calculate unique kmer pairs from a Jellyfish or KMC dump file.')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Jellyfish or KMC input dump file (stdin).')
  parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='kmerpairs')
  parser.add_argument('-k', help='The length of the kmer.', default=21)
  parser.add_argument('-t', help='Number of processes to use.', default = 4)

  args = parser.parse_args()
  dumps_file = args.infile
  output_pattern = args.o
  k = int(args.k)
  num_processes = int(args.t)

  file_one_away_pairs = open(output_pattern + '_one_away_pairs.tsv', 'w')
  file_coverages = open(output_pattern + '_coverages_2.tsv', 'w')

  #Initialize dictionaries in which the key is a kmer_half (kmer_L or kmer_R respectively), and the value is a list of (other_kmer_half, index) pairs.
  #kmer_L_to_index_family = defaultdict(list)
  kmer_R_to_index_family = defaultdict(list)

  #Get the locations for the two halves of the kmer.
  k_L = k // 2
  i_L_L = 0
  i_L_R = k_L - 1
  i_R_L = k_L
  i_R_R = k-1

  # Read each line of the input file in order to load the kmers and coverages and process the kmer halves.
  current_kmer_L = ""
  for i1, line in enumerate(dumps_file):
    kmer, coverage1 = line.split()
    coverage1 = int(coverage1)

    #coverages.append(coverage)
    new_kmer_L = kmer[i_L_L:i_L_R+1]
    kmer_R = kmer[i_R_L+1:i_R_R+1]
    if new_kmer_L == current_kmer_L:
      if kmer_R in kmer_R_to_index_family:
        for i2,coverage2 in kmer_R_to_index_family[kmer_R]:
          if coverage2 < coverage1:
            file_one_away_pairs.write(str(i2) + '\t' + str(i1) + '\n')
            file_coverages.write(str(coverage2) + '\t' + str(coverage1) + '\n')
          else:
            file_one_away_pairs.write(str(i1) + '\t' + str(i2) + '\n')
            file_coverages.write(str(coverage1) + '\t' + str(coverage2) + '\n')
    else:
      current_kmer_L = new_kmer_L
      kmer_R_to_index_family = defaultdict(list)
    kmer_R_to_index_family[kmer_R].append((i1,coverage1))

  print('Compuation done.')
  file_one_away_pairs.close()
  file_coverages.close()


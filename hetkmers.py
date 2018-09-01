#!/usr/bin/env python3

from collections import defaultdict
from multiprocessing import Pool
import itertools
import sys
import argparse

####################
# GLOBAL VARIABLES #
####################

letter_to_int = {'A':0, 'C':1, 'G':2, 'T':3}
convertString = "ACGT"

###################
###  DEFINE     ###
###   FUNCTIONS ###
###################

def kmer_to_int(kmer):
  """This function converts a kmer string into an integer."""
  k = len(kmer)
  return sum((4**(k-i-1) * letter_to_int[letter] for i, letter in enumerate(kmer)))

def int_to_kmer(n):
  """This function converts an integer into a kmer string. The variable k refers to the kmer length."""
  def int_to_truncated_kmer(n):
    if n < 4:
      return convertString[n]
    else:
      return int_to_truncated_kmer(n//4) + convertString[n%4]
  truncated_kmer = int_to_truncated_kmer(n)
  prefix = 'A'*(k - len(truncated_kmer))
  return prefix + truncated_kmer

def get_1away_pairs(kmer_ints, local_indices, i_L, i_R):
  """kmer_ints is the list of kmer_ints. local_indices is a list of indices of the kmer portions currently under consideration. i_L and i_R are the start and end locations of the kmer portions currently under consideration (0-indexed and inclusive). get_1away_pairs returns a list of pairs of indices where each pair of indices corresponds to a pair of kmers different in exactly one base."""

  #This is the base case for the recursion. Return every pair of indices where the kmer portions (between locations i_L and i_R inclusive) at those indices differ at exactly one base.
  k = i_R - i_L + 1
  if k == 1:
    return [(i,j) for (i,j) in itertools.combinations(local_indices, 2) if int_to_kmer(kmer_ints[i])[i_L:i_R+1] != int_to_kmer(kmer_ints[j])[i_L:i_R+1]]

  #Get the locations for the two halves of the kmer portion
  k_L = k // 2
  i_L_L = i_L
  i_L_R = i_L + k_L - 1
  i_R_L = i_L + k_L
  i_R_R = i_R

  #initialize dictionaries in which the key is half of the kmer portion (kmer_L or kmer_R respectively), and the value is a list of indices ("index_family") corresponding to the kmer portions with that same half
  kmer_L_to_index_family = defaultdict(list)
  kmer_R_to_index_family = defaultdict(list)

  #initialize pairs, which will be returned by get_1away_pairs
  pairs = []

  #for each index, calculate the corresponding left half and right half, then add its index to the corresponding entries of the dictionary
  for i in local_indices:
    kmer_L = int_to_kmer(kmer_ints[i])[i_L_L:i_L_R+1]
    kmer_R = int_to_kmer(kmer_ints[i])[i_R_L:i_R_R+1]
    kmer_L_to_index_family[kmer_to_int(kmer_L)] += [i]
    kmer_R_to_index_family[kmer_to_int(kmer_R)] += [i]

  #discard kmer portion halves, keep indices
  kmer_L_index_families = kmer_L_to_index_family.values()
  del kmer_L_to_index_family
  kmer_R_index_families = kmer_R_to_index_family.values()
  del kmer_R_to_index_family

  #for each left half in which there are multiple kmer portions with that left half, find the list of pairs in which the right half differs by 1. (aka, if left half matches, recurse on right half).
  for kmer_L_index_family in kmer_L_index_families: #same in left half
    if len(kmer_L_index_family) > 1:
      pairs += get_1away_pairs(kmer_ints, kmer_L_index_family, i_R_L, i_R_R) #differ by 1 in right half

  #for each right half in which there are multiple kmer portions with that same right half, find the list of pairs in which the left half differs by 1. (aka, if right half matches, recurse on left half).
  for kmer_R_index_family in kmer_R_index_families: #same in right half
    if len(kmer_R_index_family) > 1:
      pairs += get_1away_pairs(kmer_ints, kmer_R_index_family, i_L_L, i_L_R) #differ by 1 in left half
  return(pairs)

##############
### SCRIPT ###
##############
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Calculate unique kmer pairs from a Jellyfish or KMC dump file.')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Jellyfish input dump file (stdin).')
  parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='kmerpairs')
  parser.add_argument('-k', help='The length of the kmer.', default=21)
  parser.add_argument('-t', help='Number of processes to use.', default = 5)                    

  args = parser.parse_args()
  dumps_file = args.infile
  output_pattern = args.o
  k = int(args.k)
  num_processes = int(args.t)

  #Get the locations for the two halves of the kmer
  k_L = k // 2
  i_L_L = 0
  i_L_R = k_L - 1
  i_R_L = k_L
  i_R_R = k-1

  #initiate list of kmer_ints and coverages
  kmer_ints = []
  coverages = []

  #initialize dictionaries in which the key is half of the kmer (kmer_L or kmer_R respectively), and the value is a list of indices ("index_family") corresponding to the kmers with that same half
  kmer_L_to_index_family = defaultdict(list)
  kmer_R_to_index_family = defaultdict(list)

  #read each line of the input file in order to load the kmers and coverages, and process the kmer halves
  for i, line in enumerate(dumps_file):
    kmer, coverage = line.split()
    kmer_int = kmer_to_int(kmer)
    coverage = int(coverage)
    kmer_L = kmer[i_L_L:i_L_R+1]
    kmer_R = kmer[i_R_L:i_R_R+1]
    kmer_L_int = kmer_to_int(kmer_L)
    kmer_R_int = kmer_to_int(kmer_R)
    kmer_ints.append(kmer_int)
    coverages.append(coverage)
    kmer_L_to_index_family[kmer_L_int] += [i]
    kmer_R_to_index_family[kmer_R_int] += [i]

  print('Kmers and coverages loaded. Initial kmer halves processed.')

  #discard kmer halves, keep indices
  kmer_L_index_families = kmer_L_to_index_family.values()
  del kmer_L_to_index_family
  kmer_R_index_families = kmer_R_to_index_family.values()
  del kmer_R_to_index_family

  print('Starting processes.')

  #create process pool to process the kmer half index families
  with Pool(num_processes) as p:
    pairs1 = p.starmap(get_1away_pairs, ((kmer_ints, kmer_L_index_family, i_R_L, i_R_R) for kmer_L_index_family in kmer_L_index_families if len(kmer_L_index_family)>1))
    pairs2 = p.starmap(get_1away_pairs, ((kmer_ints, kmer_R_index_family, i_L_L, i_L_R) for kmer_R_index_family in kmer_R_index_families if len(kmer_R_index_family)>1))
  one_away_pairs = list(itertools.chain(itertools.chain.from_iterable(pairs1), itertools.chain.from_iterable(pairs2)))

  print('Process pool completed.')

  #save one_away_pairs to a tsv file, keep track of which kmers only occur once
  repeated = {}
  with open(output_pattern + '_one_away_pairs.tsv', 'w') as record_file:
    for (i1, i2) in one_away_pairs:
      record_file.write(str(i1) + '\t' + str(i2) + '\n')
      repeated[i1] = i1 in repeated
      repeated[i2] = i2 in repeated

  print('*_one_away_pairs.tsv file saved.')

  #save families_2 and coverages_2 tsv files, which only include one_away_pairs that don't overlap any others
  with open(output_pattern + '_families_2.tsv', 'w') as record_file1, open(output_pattern + '_coverages_2.tsv', 'w') as record_file2:
    for (i1, i2) in one_away_pairs:
      if not repeated[i1] and not repeated[i2]:
        cov1 = coverages[i1]
        cov2 = coverages[i2]
        if cov1 < cov2:
          record_file1.write(str(i1) + '\t' + str(i2) + '\n')
          record_file2.write(str(cov1) + '\t' + str(cov2) + '\n')
        else:
          record_file1.write(str(i2) + '\t' + str(i1) + '\n')
          record_file2.write(str(cov2) + '\t' + str(cov1) + '\n')

  print('*_families_2.tsv and *_coverages_2.tsv files saved.')

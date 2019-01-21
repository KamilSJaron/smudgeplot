#!/usr/bin/env python3

from collections import defaultdict
import multiprocessing as mp
import argparse
import itertools
import sys

###################
###  DEFINE     ###
###   FUNCTIONS ###
###################

def get_one_away_pairs(kmer_index_family, k):
  """kmer_index_family is a list of (kmer, index) pairs currently under consideration. k is the kmer length. get_one_away_pairs returns a list of pairs of indices where each pair of indices corresponds to a pair of kmers different in exactly one base."""

  #This is the base case for the recursion. Return every pair of indices where the kmers corresponding to those indices differ at exactly one base.
  if k == 1:
    return [(i,j) for ((kmer1,i),(kmer2,j)) in itertools.combinations(kmer_index_family, 2) if kmer1 != kmer2]

  #Initialize one_away_pairs, which will be returned by get_one_away_pairs.
  one_away_pairs = []

  #Initialize dictionaries in which the key is a kmer_half (kmer_L or kmer_R) and the value is a list of (other_kmer_half, index) pairs.
  kmer_L_to_index_family = defaultdict(list)
  kmer_R_to_index_family = defaultdict(list)

  #Get the locations for the two halves of the kmer.
  k_L = k // 2
  k_R = k-k_L
  i_L_L = 0
  i_L_R = k_L - 1
  i_R_L = k_L
  i_R_R = k-1

  #For each kmer and index calculate the corresponding left half and right half, then add the necessary (kmer_half, index) pair to the corresponding entries of the dictionary
  for kmer, i in kmer_index_family:
    kmer_L = kmer[i_L_L:i_L_R+1]
    kmer_R = kmer[i_R_L:i_R_R+1]
    kmer_L_to_index_family[kmer_L].append((kmer_R, i))
    kmer_R_to_index_family[kmer_R].append((kmer_L, i))

  #For each left half in which there are multiple kmers with that left half, find the list of pairs in which the right half differs by 1. (aka, if left half matches, recurse on right half).
  for kmer_L_index_family in kmer_L_to_index_family.values(): #same in left half
    if len(kmer_L_index_family) > 1:
      one_away_pairs.extend(get_one_away_pairs(kmer_L_index_family, k_R)) #differ by 1 in right half

  del kmer_L_to_index_family

  #For each right half in which there are multiple kmers with that same right half, find the list of pairs in which the left half differs by 1. (aka, if right half matches, recurse on left half).
  for kmer_R_index_family in kmer_R_to_index_family.values(): #same in right half
    if len(kmer_R_index_family) > 1:
      one_away_pairs.extend(get_one_away_pairs(kmer_R_index_family, k_L)) #differ by 1 in left half

  del kmer_R_to_index_family

  return(one_away_pairs)

def worker(q, results, k):
  while True:
    kmer_index_family = q.get()
    if kmer_index_family is None:
      results.put(None)
      break
    results.put(get_one_away_pairs(kmer_index_family, k))
  q.close()
  results.close()

#########################################
# WRAPPING FUNCTIONS OF THE TWO MODULES #
#########################################

def middle_one_away(args):
  dumps_file = args.infile
  output_pattern = args.o
  k = int(args.k)

  file_one_away_pairs = open(output_pattern + '_one_away_pairs.tsv', 'w')
  file_coverages = open(output_pattern + '_coverages.tsv', 'w')

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

  file_one_away_pairs.close()
  file_coverages.close()

def all_one_away(args):
  dumps_file = args.infile
  output_pattern = args.o
  k = int(args.k)
  num_processes = int(args.t)

  #Initiate coverages list.
  coverages = []

  #Initiate one_away_pairs list.
  one_away_pairs = []

  #Initialize dictionaries in which the key is a kmer_half (kmer_L or kmer_R respectively), and the value is a list of (other_kmer_half, index) pairs.
  kmer_L_to_index_family = defaultdict(list)
  kmer_R_to_index_family = defaultdict(list)

  #Get the locations for the two halves of the kmer.
  k_L = k // 2
  i_L_L = 0
  i_L_R = k_L - 1
  i_R_L = k_L
  i_R_R = k-1

  #Read each line of the input file in order to load the kmers and coverages and process the kmer halves.
  for i, line in enumerate(dumps_file):
    kmer, coverage = line.split()
    coverage = int(coverage)
    coverages.append(coverage)
    kmer_L = kmer[i_L_L:i_L_R+1]
    kmer_R = kmer[i_R_L:i_R_R+1]
    kmer_L_to_index_family[kmer_L].append((kmer_R, i))
    kmer_R_to_index_family[kmer_R].append((kmer_L, i))

  print('Kmers and coverages loaded. Initial kmer halves processed. Starting processes.')

  #Create input queue and results queue to process the kmer index families.
  mp.set_start_method('spawn')
  q = mp.Queue()
  results = mp.Queue()

  #Initiate processes.
  processes = []
  for i in range(num_processes):
    p = mp.Process(target = worker, args = (q, results, k))
    p.start()
    processes.append(p)

  print('Processes Started.')

  #Add kmer index families to the queue.
  for kmer_L_index_family in kmer_L_to_index_family.values():
    if len(kmer_L_index_family) > 1:
      q.put(kmer_L_index_family)

  del kmer_L_to_index_family

  for kmer_R_index_family in kmer_R_to_index_family.values():
    if len(kmer_R_index_family) > 1:
      q.put(kmer_R_index_family)

  del kmer_R_to_index_family

  #Add sentinels so the processes know when to terminate.
  for i in range(num_processes):
    q.put(None)

  q.close()
  print('All index families added to queue.')

  #Save one_away_pairs to a tsv file, and keep track of which indices only occur once.
  num_completed = 0
  repeated = {}
  with open(output_pattern + '_one_away_pairs.tsv', 'w') as record_file:
    while True:
      partial_one_away_pairs = results.get()
      if partial_one_away_pairs is None:
        num_completed += 1
        if num_completed == num_processes:
          break
      else:
        for (i1, i2) in partial_one_away_pairs:
          one_away_pairs.append((i1, i2))
          record_file.write(str(i1) + '\t' + str(i2) + '\n')
          repeated[i1] = i1 in repeated
          repeated[i2] = i2 in repeated

  print('*_one_away_pairs.tsv file saved.')

  #Wait for the processes to finish.
  for p in processes:
    p.join()

  results.close()
  print('Processes joined.')

  #Save families_2 and coverages_2 tsv files, which only include one_away_pairs that don't overlap any others.
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

##############
### SCRIPT ###
##############

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Calculate unique kmer pairs from a Jellyfish or KMC dump file.')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Alphabetically sorted Jellyfish or KMC dump file (stdin).')
  parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='kmerpairs')
  parser.add_argument('-k', help='The length of the kmer.', default=21)
  parser.add_argument('-t', help='Number of processes to use.', default = 4)
  parser.add_argument('--all', dest='all', action='store_const', const = True, default = False,
                      help='Get all kmer pairs one SNP away from each other (default: just the middle one).')

  args = parser.parse_args()

  if args.all:
    all_one_away(args)
  else:
    middle_one_away(args)

  print('Done!')
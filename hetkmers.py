#!/usr/bin/env python

from collections import defaultdict, Counter
from itertools import combinations
# import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pickle
import sys
import argparse

####################
# GLOBAL VARIABLES #
####################

# every time you called kmer_to_int function you created a new dictionary
# that's costly and you can do it just once
letter_to_int = {'A':1, 'C':2, 'G':3, 'T':4}

###################
###  DEFINE     ###
###   FUNCTIONS ###
###################

def kmer_to_int(kmer):
  """This function turns a kmer string into an integer."""
  k = len(kmer)
  # souradnice = strtol(zkoumaneslovo, NULL, 4);
  return sum([4**(k-i-1) * letter_to_int[letter] for i, letter in enumerate(kmer)])

def get_1away_pairs(local_index_to_kmer, k):
  """local_index_to_kmer is a dictionary where the value is a kmer portion, and the key is the index of the original kmer in which the kmer portion is found. get_1away_pairs returns a list of pairs of indices where each pair of indices correspond to a pair of kmers different in exactly one base."""

  #This is the base case for the recursion. Return every pair of indices where the kmers at those indices differ at exactly one base.
  if k == 1:
    return [(i,j) for (i,j) in combinations(local_index_to_kmer, 2) if local_index_to_kmer[i][0] != local_index_to_kmer[j][0]]

  #Get the two halves of the kmer
  k_L = k // 2
  k_R = k-k_L

  #initialize dictionaries in which the key is the hash of half of the kmer, and the value is a list of indices of the kmers with that same hash
  kmer_L_hashes = defaultdict(list)
  kmer_R_hashes = defaultdict(list)

  #initialize pairs, which will be returned by get_1away_pairs
  pairs = []

  #for each kmer, calculate its left hash and right hash, then add its index to the corresponding entries of the dictionary
  for i, kmer in local_index_to_kmer.items():
    kmer_L = kmer[:k_L]
    kmer_R = kmer[k_L:]
    kmer_L_hashes[kmer_to_int(kmer_L)] += [i]
    kmer_R_hashes[kmer_to_int(kmer_R)] += [i]

  #for each left hash in which there are multiple kmers with that left hash, find the list of pairs in which the right half differs by 1. (aka, if left half matches, recurse on right half).
  for kmer_L_hash_indices in kmer_L_hashes.values(): #same in first half
    if len(kmer_L_hash_indices) > 1:
      pairs += get_1away_pairs({kmer_L_hash_index:local_index_to_kmer[kmer_L_hash_index][k_L:] for kmer_L_hash_index in kmer_L_hash_indices}, k_R) #differ by 1 in right half

  #for each right hash in which there are multiple kmers with that right hash, find the list of pairs in which the left half differs by 1. (aka, if right half matches, recurse on left half).
  for kmer_R_hash_indices in kmer_R_hashes.values(): #same in right half
    if len(kmer_R_hash_indices) > 1:
      pairs += get_1away_pairs({kmer_R_hash_index:local_index_to_kmer[kmer_R_hash_index][:k_L] for kmer_R_hash_index in kmer_R_hash_indices}, k_L) #differ by 1 in left half
  return(pairs)

##############
### SCRIPT ###
##############

parser = argparse.ArgumentParser(description='Calculate unique kmer pairs from a Jellyfish dump file.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                     default=sys.stdin, help='Jellyfish input dump file (stdin)')
parser.add_argument('-o', help='the pattern used to name the output (kmerpairs)',
                    default='kmerpairs')

args = parser.parse_args()
dumps_file = args.infile
output_pattern = args.o
# the output string is

#create a pandas dataframe of the dumps file
df = pd.read_csv(dumps_file, sep = ' ', header = None)

#get a list of kmers and coverages
kmers = list(df[0])
coverages = list(df[1])

#delete the dataframe
del df

#get the value of k
k = len(kmers[0])

#turn the kmers and coverages into a dictionary, where the key is the index of the original kmer
test_kmers = {i:kmer for (i, kmer) in enumerate(kmers)}
test_coverages = {i:coverage for (i, coverage) in enumerate(coverages)}

#delete the kmers and coverages lists
del kmers
del coverages

#get the one away pairs and save it to a pickle file
one_away_pairs = get_1away_pairs(test_kmers, k)
pickle.dump(one_away_pairs, open(output_pattern + '_one_away_pairs.p', 'wb'))

#get the two away pairs and save it to a pickle file
#two_away_pairs = get_2away_pairs(test_kmers, k)
#pickle.dump(two_away_pairs, open('two_away_pairs.p', 'wb'))

#To determine the families, use the "connected_components" function of a networkx graph
#add each of the pairs as edges to the graph, and save it to a pickle file
G = nx.Graph()
for one_away_pair in one_away_pairs:
  G.add_edge(*one_away_pair)
#for two_away_pair in two_away_pairs:
#  G.add_edge(*two_away_pair)
# pickle.dump(G, open(output_pattern + 'G.p', 'wb'))

#get the lengths of the families, save it to a pickle file, and then look at the Counter summary
component_lengths = [len(x) for x in nx.connected_components(G)]
# pickle.dump(component_lengths, open(output_pattern + '_component_lengths.p', 'wb'))
Counter(component_lengths)

#get the families with 2-6 members
families_2 = []
families_3 = []
families_4 = []
families_5 = []
families_6 = []

for component in nx.connected_components(G):
  if len(component) == 2:
    families_2.append(list(component))
  elif len(component) == 3:
    families_3.append(list(component))
  elif len(component) == 4:
    families_4.append(list(component))
  elif len(component) == 5:
    families_5.append(list(component))
  elif len(component) == 6:
    families_6.append(list(component))

#save families_2 to a pickle file
# pickle.dump(families_2, open(output_pattern + '_families_2.p', 'wb'))
# perhaps faster way how to calculate coverages_2, probably more pythonic
# coverages_2 = [test_coverages[cov_i1] + test_coverages[cov_i2] for cov_i1, cov_i2 in families_2]
coverages_2 = [sum([test_coverages[i] for i in family]) for family in families_2]
# pickle.dump(coverages_2, open(output_pattern + '_coverages_2.p', 'wb'))

with open(output_pattern + '_coverages_2.tsv', 'w') as record_file:
  for (cov_i1, cov_i2) in families_2:
    cov1 = test_coverages[cov_i1]
    cov2 = test_coverages[cov_i2]
    if cov1 < cov2:
      record_file.write(str(cov1) + '\t' + str(cov2) + '\n')
    else:
      record_file.write(str(cov2) + '\t' + str(cov1) + '\n')

# save families_3 to a pickle file
# pickle.dump(families_3, open(output_pattern + '_families_3.p', 'wb'))
# coverages_3 = [sum([test_coverages[i] for i in family]) for family in families_3]
#
# # save families_4 to a pickle file
# pickle.dump(families_4, open(output_pattern + '_families_4.p', 'wb'))
# coverages_4 = [sum([test_coverages[i] for i in family]) for family in families_4]
#
# # save families_5 to a pickle file
# pickle.dump(families_5, open(output_pattern + '_families_5.p', 'wb'))
# coverages_5 = [sum([test_coverages[i] for i in family]) for family in families_5]
#
# # save families_6 to a pickle file
# pickle.dump(families_6, open(output_pattern + '_families_6.p', 'wb'))
# coverages_6 = [sum([test_coverages[i] for i in family]) for family in families_6]
#

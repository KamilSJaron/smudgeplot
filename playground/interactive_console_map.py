import argparse
from Bio import SeqIO
from PySAIS import sais
import numpy as np
import gzip
import pandas
import time
from bisect import bisect_right
from smudgeplot.map import mapper
# argument 1
# args.genomefile
kmer_genome_file = 'data/Tps1/1_Tps_b3v07.fa'
# kmer_fasta_file = 'data/Tps1/Tps_middle_pair_end_reads_kmers_in_smudge_1_sample_50.fasta'
output_pattern = 'data/Tps1/Tps_middle_pair_end_reads_kmers_in_smudge_1_sample_50.bed'

parser = argparse.ArgumentParser(description='whatever')
args = parser.parse_args()
args.genomefile = kmer_genome_file
args.o = 'data/Tps1/Tps_middle_pair_end_reads'
args.s = 0

kmer_map = mapper(args)
kmer_map.loadGenome()

##########################
###### SUFFIX ARRAY ######
##########################

kmer_map.constructSuffixArray()

kmer_file_name_s1 = 'data/Tps1/Tps_middle_pair_end_reads_kmers_in_smudge_1.txt'
kmer_file_name_s1
with open(kmer_file_name_s1, 'r') as s1_kmer_file:
    s1_kmers = [kmer.rstrip() for kmer in s1_kmer_file]


start = time.time()
mapping_list = [kmer_map.searchKmer(kmer) for kmer in s1_kmers]
end = time.time()
print('Done in ' + str(round(end - start, 1)) + ' s')

hit_numbers = [len(mapped_kmer) for mapped_kmer in mapping_list]
hist = Counter(hit_numbers)
fractions = [100 * round(hist[i] / len(mapping_list), 4) for i in range(0,10)]

import logging

logging.info("Missing in the assembly: " + str(fractions[0]) + "% of duplicated kmers.")
logging.info("Collapsed in the assembly: " + str(fractions[1]) + "% of duplicated kmers.")
logging.info("Correctly assembled: " + str(fractions[2]) + "% of duplicated kmers.")
logging.info("Assembled three times (once more than it should): " + str(fractions[3]) + "% of duplicated kmers.")
logging.info("Assembled more than three times: " + str(1 - sum(fractions[0:3])) + "% of duplicated kmers.")

import matplotlib.pyplot as plt
# plt.hist(, bins='auto')
bins = range(0, 10)
plt.hist(hit_numbers, bins=bins)
plt.show()

# # # # # # # # # # # # # #
# # scaffold processing # #
# # # # # # # # # # # # # #

# this will be somehow integrate in the class, but for now...

# how I can easily get scf length from scaffold name
scf2size = dict()
for i, seq in enumerate(kmer_map.sequences):
    scf2size[kmer_map.scf_names[i]] = len(seq)

### correctly assembled duplications
from collections import defaultdict
scaffold_dupl = defaultdict(int)
within_scaffold_dupl = defaultdict(int)
for mapped_kmer in mapping_list:
    if len(mapped_kmer) == 2:
        if mapped_kmer[0][0] == mapped_kmer[1][0]:
            within_scaffold_dupl[mapped_kmer[0][0]] += 1
        else:
            scaffold_dupl[mapped_kmer[0][0]] += 1
            scaffold_dupl[mapped_kmer[1][0]] += 1

### Collapsed assembled duplications
collapsed_dupl = defaultdict(int)
for mapped_kmer in mapping_list:
    if len(mapped_kmer) == 1:
        collapsed_dupl[mapped_kmer[0][0]] += 1


with open(kmer_map.output_pattern + "_scaffold_duplicates_list.txt", 'w') as dupl_file:
    for scf in kmer_map.scf_names:
        dupl_file.write(scf + "\t" + str(scaffold_dupl[scf]) + "\t" + str(scaffold_dupl[scf] / scf2size[scf]) + "\n")

with open(kmer_map.output_pattern + "_scaffold_collapsed_duplicates_list.txt", 'w') as dupl_file:
    for scf in kmer_map.scf_names:
        dupl_file.write(scf + "\t" + str(collapsed_dupl[scf]) + "\t" + str(collapsed_dupl[scf] / scf2size[scf]) + "\n")

### distributions

# ????

### Networks

dupl_network = defaultdict(int)
for mapped_kmer in mapping_list:
    if len(mapped_kmer) == 2:
        scf1 = mapped_kmer[0][0]
        scf2 = mapped_kmer[1][0]
        if scf1 < scf2:
            dupl_network[(scf1, scf2)] += 1
        else :
            dupl_network[(scf2, scf1)] += 1

## save it???

# # # # # # # # #
# # full kmer # #
# # # # # # # # #

l = 0
r = len(sa) - 1
while l <= r:
    m = (l + r) // 2
    eval_pos = sa[m]
    genome_kmer = genome[eval_pos:(eval_pos + 21)]
    if kmer == genome_kmer:
        print('jackpot at index: ', m, ' genomic position: ', eval_pos)
        break
    elif genome_kmer < kmer:
        l = m + 1
    else:
        r = m - 1
# This works
# jackpot at index:  421217008  genomic position:  480
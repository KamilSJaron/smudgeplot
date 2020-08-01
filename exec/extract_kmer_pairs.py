#!/usr/bin/env python3

#extract kmer pairs based on total coverage and minor allele frequency
"""
python extract_kmer_pairs_based_on_kmer_sum_and_ratio.py -cov <coverage_file> -seq <kmer_pair_seq_file> -minc <from_sum> -maxc <to_sum> -minr <from_ratio> -maxr <to_ratio>
python extract_kmer_pairs_based_on_kmer_sum_and_ratio.py -cov kmer_pairs_coverages.tsv -seq kmer_pairs_sequences.tsv -minc 19 -maxc 32.65 -minr 0.325 -maxr 0.5

1. what does the paired coverage input file look like?
head -n 2 kmer_pairs_coverages.tsv
28	55
15	19
2. what does the paired sequence input file look like? (1st and 2nd columns are the pairs, 3rd is the consensus)
AAAAAAAAAAAAAAAAAAGAAGAAGAAGTCT	AAAAAAAAAAAAAAAGAAGAAGAAGAAGTCT	AAAAAAAAAAAAAAANAAGAAGAAGAAGTCT
AAAAAAAAAAAAAAAGAAGGAACCAACCCTG	AAAAAAAAAAAAAAAAAAGGAACCAACCCTG	AAAAAAAAAAAAAAANAAGGAACCAACCCTG

0. what do you want output files look like?
#a)first output file -
head -n 1 good_smudge_pair_cov.tsv
15	19

#b) second output file -
head -n 1 good_smudge_pair_seq.tsv
AAAAAAAAAAAAAAAGAAGGAACCAACCCTG	AAAAAAAAAAAAAAAAAAGGAACCAACCCTG	AAAAAAAAAAAAAAANAAGGAACCAACCCTG
"""

import sys, argparse
from collections import defaultdict

#add this later when it runs on a server
ap = argparse.ArgumentParser()
ap.add_argument("-cov", "--coverageFile",required=True, help="coverage file for the kmer pairs")
ap.add_argument("-seq", "--seqFile",required=True, help="sequences of the kmer pairs")
ap.add_argument("-minc", "--countMin",required=True, help="lower bound of the summed coverage", type=int)
ap.add_argument("-maxc", "--countMax",required=True, help="upper bound of the summed coverage", type=int)
ap.add_argument("-minr", "--ratioMin",required=True, help="lower bound of minor allele ratio", type=float)
ap.add_argument("-maxr", "--ratioMax",required=True, help="upper bound of minor allele ratio", type=float)


args = ap.parse_args()

# out_cov_h = open("good_smudge_pair_cov.tsv","w")
# out_seq_h = open("good_smudge_pair_seq.tsv","w")
index2covs = defaultdict(list)
good_index_list = list()

with open(args.coverageFile,"r") as f:
	for i, row in enumerate(f):
		row_list = row.rstrip("\n").split("\t")
		countL = int(row_list[0])
		countR = int(row_list[1])

		cov_sum = countL + countR
		cov_ratio = countL / (countL + countR)
		if cov_sum <= args.countMax and cov_sum >= args.countMin and cov_ratio >= args.ratioMin and cov_ratio <= args.ratioMax:
			# out_cov_h.write(str(countL)+"\t" + str(countR) + "\n")
			#print(str(i1))
			good_index_list.append(i)
			index2covs[i] = [countL, countR]

index_seq_dict = dict()

with open(args.seqFile,"r") as f:
	for i, row in enumerate(f):
		if len(index2covs[i]) == 2:
			kmer1, kmer2 = row.rstrip('\n').split('\t')
			sys.stdout.write(">kmer_" + str(i) + "_1_cov_" + str(index2covs[i][0]) + "\n")
			sys.stdout.write(kmer1 + '\n')
			sys.stdout.write(">kmer_" + str(i) + "_2_cov_" + str(index2covs[i][1]) + "\n")
			sys.stdout.write(kmer2 + '\n')

# fraction = float(len(good_index_list) / float(i2))
# print("the fraction of pairs that meet the criteria is: " + str(fraction))

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

import sys, os, re, math, argparse
from collections import defaultdict

#add this later when it runs on a server
ap = argparse.ArgumentParser()
ap.add_argument("-cov", "--coverageFile",required=True, help="coverage file for the kmer pairs")
ap.add_argument("-seq", "--seqFile",required=True, help="sequences of the kmer pairs")
ap.add_argument("-minc", "--countMin",required=True, help="lower bound of the summed coverage")
ap.add_argument("-maxc", "--countMax",required=True, help="upper bound of the summed coverage")
ap.add_argument("-minr", "--ratioMin",required=True, help="lower bound of minor allele ratio")
ap.add_argument("-maxr", "--ratioMax",required=True, help="upper bound of minor allele ratio")


args = vars(ap.parse_args())
covFile = args["coverageFile"]
seqFile = args["seqFile"]
minC = float("{:.3f}".format(float(args["countMin"])))
maxC = float("{:.3f}".format(float(args["countMax"]))) 
minR = float("{:.3f}".format(float(args["ratioMin"]))) 
maxR = float("{:.3f}".format(float(args["ratioMax"]))) 

out_cov_h = open("good_smudge_pair_cov.tsv","w")
out_seq_h = open("good_smudge_pair_seq.tsv","w")
good_index_list = list()
i1=0
with open(covFile,"r") as f:
	for row in f:
		i1+=1
		row = row.rstrip("\n")
		row_list = re.split("\t",row)
		countL = int(row_list[0])
		countR = int(row_list[1])

		#this can be unnecessary as later I found that the kmer coverages from kmc is already sorted in that the left coverage is always smaller than right coverage
		if countL <= countR:
			countS = float("{:.3f}".format(float(countL))) 
			countB = float("{:.3f}".format(float(countR))) 
		else:
			countS = float("{:.3f}".format(float(countR))) 
			countB = float("{:.3f}".format(float(countL))) 

		if countS + countB <= maxC and countS+countB >= minC and float(countS)/(float(countS)+float(countB)) >= minR and float(countS)/(float(countS)+float(countB)) <= maxR:
			out_cov_h.write(str(countS)+"\t" + str(countB) + "\n")
			#print(str(i1))
			good_index_list.append(i1)
		else:
			pass

index_seq_dict = dict()

i2=0
with open(seqFile,"r") as f:
	for row in f:
		i2+=1
		row = row.rstrip("\n")
		index_seq_dict[i2] =  row


for good_i in good_index_list:
	out_seq = index_seq_dict[good_i]
	out_seq_h.write(out_seq+"\n")

fraction = float(len(good_index_list)/float(i2))
print("the fraction of pairs that meet the criteria is: " + str(fraction))

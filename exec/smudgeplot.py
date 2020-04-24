#!/usr/bin/env python3

import argparse
import sys
from os import system
from math import log
from math import ceil
import numpy as np
from scipy.signal import argrelextrema
# hetkmers dependencies
from collections import defaultdict
from itertools import combinations

version = '0.2.3dev'

############################
# processing of user input #
############################

class parser():
    def __init__(self):
        argparser = argparse.ArgumentParser(
            # description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
            usage='''smudgeplot <task> [options] \n
tasks: cutoff    Calculate meaningful values for lower/upper kmer histogram cutoff.
       hetkmers  Calculate unique kmer pairs from a Jellyfish or KMC dump file.
       plot      Generate 2d histogram; infere ploidy and plot a smudgeplot.\n\n''')
        argparser.add_argument('task', help='Task to execute; for task specific options execute smudgeplot <task> -h')
        argparser.add_argument('-v', '--version', action="store_true", default = False, help="print the version and exit")
        # print version is a special case
        if len(sys.argv) > 1:
            if sys.argv[1] in ['-v', '--version']:
                self.task = "version"
                return
            # the following line either prints help and die; or assign the name of task to variable task
            self.task = argparser.parse_args([sys.argv[1]]).task
        else:
            self.task = ""
        # if the task is known (i.e. defined in this file);
        if hasattr(self, self.task):
            # load arguments of that task
            getattr(self, self.task)()
        else:
            argparser.print_usage()
            sys.stderr.write('"' + self.task + '" is not a valid task name\n')
            exit(1)

    def hetkmers(self):
        '''
        Calculate unique kmer pairs from a Jellyfish or KMC dump file.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot hetkmers',
            description='Calculate unique kmer pairs from a Jellyfish or KMC dump file.')
        argparser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Alphabetically sorted Jellyfish or KMC dump file (stdin).')
        argparser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='kmerpairs')
        argparser.add_argument('--middle', dest='middle', action='store_const', const = True, default = False,
                               help='Get all kmer pairs that are exactly the same but in the middle nt. When this flag is used, the input dump must be alphabetically sorted/ (default: different by a SNP at any position).')
        self.arguments = argparser.parse_args(sys.argv[2:])

    def plot(self):
        '''
        Generate 2d histogram; infer ploidy and plot a smudgeplot.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot plot', description='Generate 2d histogram for smudgeplot')
        argparser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages (default \"coverages_2.tsv\")."')
        argparser.add_argument('-o', help='The pattern used to name the output (smudgeplot).', default='smudgeplot')
        argparser.add_argument('-q', help='Remove kmer pairs with coverage over the specified quantile; (default none).', type=float, default=1)
        argparser.add_argument('-L', help='The lower boundary used when dumping kmers (default min(total_pair_cov) / 2).', type=int, default=0)
        argparser.add_argument('-n', help='The expected haploid coverage (default estimated from data).', type=float, default=0)
        argparser.add_argument('-t', '--title', help='name printed at the top of the smudgeplot (default none).', default='')
        # argparser.add_argument('-m', '-method', help='The algorithm for annotation of smudges (default \'local_aggregation\')', default='local_aggregation')
        argparser.add_argument('-nbins', help='The number of nbins used for smudgeplot matrix (nbins x nbins) (default autodetection).', type=int, default=0)
        argparser.add_argument('-k', help='The length of the kmer.', default=21)
        # argparser.add_argument('-kmer_file', help='Name of the input files containing kmer seuqences (assuming the same order as in the coverage file)', default = "")
        argparser.add_argument('--homozygous', action="store_true", default = False, help="Assume no heterozygosity in the genome - plotting a paralog structure; (default False).")
        self.arguments = argparser.parse_args(sys.argv[2:])

    def cutoff(self):
        '''
        Calculate meaningful values for lower/upper kmer histogram cutoff.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot cutoff', description='Calculate meaningful values for lower/upper kmer histogram cutoff.')
        argparser.add_argument('infile', type=argparse.FileType('r'), help='Name of the input kmer histogram file (default \"kmer.hist\")."')
        argparser.add_argument('boundary', help='Which bounary to compute L (lower) or U (upper)')
        self.arguments = argparser.parse_args(sys.argv[2:])

###############
# task cutoff #
###############

def round_up_nice(x):
    digits = ceil(log(x, 10))
    if digits <= 1:
        multiplier = 10 ** (digits - 1)
    else:
        multiplier = 10 ** (digits - 2)
    return(ceil(x / multiplier) * multiplier)

def cutoff(args):
    # kmer_hist = open("data/Mflo2/kmer.hist","r")
    kmer_hist = args.infile
    hist = np.array([int(line.split()[1]) for line in kmer_hist])
    if args.boundary == "L":
        local_minima = argrelextrema(hist, np.less)[0][0]
        L = max(10, int(round(local_minima * 1.25)))
        sys.stdout.write(str(L))
    else:
        # take 99.8 quantile of kmers that are more than one in the read set
        hist_rel_cumsum = np.cumsum(hist[1:]) / np.sum(hist[1:])
        U = round_up_nice(np.argmax(hist_rel_cumsum > 0.998))
        sys.stdout.write(str(U))
    sys.stdout.flush()

############
# hetkmers #
############

def get_one_away_pairs(kmer_index_family, k):
    """kmer_index_family is a list of (kmer, index) pairs currently under consideration. k is the kmer length. get_one_away_pairs returns a list of pairs of indices where each pair of indices corresponds to a pair of kmers different in exactly one base."""

    #This is the base case for the recursion. Return every pair of indices where the kmers corresponding to those indices differ at exactly one base.
    if k == 1:
        return [(i,j) for ((kmer1,i),(kmer2,j)) in combinations(kmer_index_family, 2) if kmer1 != kmer2]

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

def middle_one_away(args):
    sys.stderr.write('Extracting kmer pairs that differ in the middle nt\n')

    # file_one_away_pairs = open(args.o + '_one_away_pairs.tsv', 'w')
    file_coverages = open(args.o + '_coverages.tsv', 'w')
    file_kmers = open(args.o + '_sequences.tsv', 'w')

    duplicated = set()
    filtered = set()

    #Initialize a dictionary in which the key is the right kmer_half (not including the middle nucleotide), and the value is a list of (index, coverage) tuples corresponding to kmers that have that particular right kmer_half.
    kmer_R_to_index_family = defaultdict(list)

    # read the first line to get the length of the kmer
    with open(args.infile.name) as dump_file:
        kmer, coverage = dump_file.readline().split()
        k = len(kmer)

    #Get the locations for the two halves of the kmer.
    k_middle = k // 2
    i_L_L = 0
    i_L_R = k_middle - 1
    i_R_L = k_middle + 1
    i_R_R = k - 1

    sys.stderr.write('Saving ' + args.o + '_coverages.tsv and ' + args.o + '_sequences.tsv files.\n')
    # Read each line of the input file in order to load the kmers and coverages and process the kmer halves.
    current_kmer_L = ""
    for i1, line in enumerate(args.infile):
        kmer, coverage1 = line.split()
        coverage1 = int(coverage1)

        new_kmer_L = kmer[i_L_L:i_L_R+1]
        kmer_R = kmer[i_R_L:i_R_R+1]
        if new_kmer_L == current_kmer_L:
            if kmer_R in kmer_R_to_index_family:
                if kmer_R in duplicated:
                    filtered.discard(kmer_R)
                else:
                    duplicated.add(kmer_R)
                    filtered.add(kmer_R)
        else:
            for kmer_R in filtered:
                (i1, coverage1), (i2, coverage2) = kmer_R_to_index_family[kmer_R]
                if coverage2 < coverage1:
                    file_coverages.write(str(coverage2) + '\t' + str(coverage1) + '\n')
                else:
                    file_coverages.write(str(coverage1) + '\t' + str(coverage2) + '\n')
                file_kmers.write(current_kmer_L + 'N' + kmer_R + '\n')
            duplicated = set()
            filtered = set()
            kmer_R_to_index_family = defaultdict(list)
            current_kmer_L = new_kmer_L
        kmer_R_to_index_family[kmer_R].append((i1,coverage1))

    file_coverages.close()
    file_kmers.close()

def all_one_away(args):
    #Initiate kmer and coverages lists.
    kmers = []
    coverages = []

    # Read each line of the input file in order to
    # load the kmers and coverages and process the kmer halves.
    for i, line in enumerate(args.infile):
        kmer, coverage = line.split()
        coverage = int(coverage)
        coverages.append(coverage)
        kmers.append(kmer)

    sys.stderr.write('Kmers and coverages loaded.\n')

    k = len(kmer) # all the kmers in the dump file have the same length, so I can just calc the number of nts in the last one
    # get_one_away_pairs is a recursive function that gatheres indices of all kmer 1 SNP from each other
    one_away_pairs = get_one_away_pairs([(kmer,i) for i,kmer in enumerate(kmers)], k)

    sys.stderr.write('Kmer pairs identified.\n')

    repeated = {}
    for (i1, i2) in one_away_pairs:
        repeated[i1] = i1 in repeated
        repeated[i2] = i2 in repeated

    sys.stderr.write('Kmers in unique kmer pairs identified.\n')

    with open(args.o + '_sequences.tsv', 'w') as file_seqs, open(args.o + '_coverages.tsv', 'w') as file_coverages:
        for (i1, i2) in one_away_pairs:
            if not repeated[i1] and not repeated[i2]:
                cov1 = coverages[i1]
                cov2 = coverages[i2]
                if cov1 < cov2:
                    file_coverages.write(str(cov1) + '\t' + str(cov2) + '\n')
                    file_seqs.write(kmers[i1] + '\t' + kmers[i2] + '\n')
                else:
                    file_coverages.write(str(cov2) + '\t' + str(cov1) + '\n')
                    file_seqs.write(kmers[i2] + '\t' + kmers[i1] + '\n')

    sys.stderr.write(args.o + '_families.tsv and ' + args.o + '_coverages.tsv files saved.\n')


#####################
# the script itself #
#####################

def main():
    _parser = parser()

    sys.stderr.write('Running smudgeplot v' + version + "\n")
    if _parser.task == "version":
        exit(0)

    sys.stderr.write('Task: ' + _parser.task + "\n")

    if _parser.task == "cutoff":
        cutoff(_parser.arguments)

    if _parser.task == "hetkmers":
        args = _parser.arguments
        if args.middle:
            middle_one_away(args)
        else :
            all_one_away(args)

    if _parser.task == "plot":
        # the plotting script is expected ot be installed in the system as well as the R library supporting it
        args = _parser.arguments
        plot_args = "-i \"" + args.infile + "\" -o \"" + args.o + "\" -k " + str(args.k)
        if args.q != 1:
            plot_args += " -q " + str(args.q)
        if args.L != 0:
            plot_args += " -L " + str(args.L)
        if args.n != 0:
            plot_args += " -n " + str(args.n)
        if args.title:
            plot_args += " -t \"" + args.title + "\""
        if args.nbins != 0:
            plot_args += " -nbins " + str(args.nbins)
        if args.homozygous:
            plot_args += " --homozygous"
        sys.stderr.write("Calling: smudgeplot_plot.R " + plot_args + "\n")
        system("smudgeplot_plot.R " + plot_args)

    sys.stderr.write("\nDone!\n")
    exit(0)

if __name__=='__main__':
    main()

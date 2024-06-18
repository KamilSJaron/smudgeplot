#!/usr/bin/env python3

import argparse
import sys
from os import system
from math import log
from math import ceil
# hetkmers dependencies
from collections import defaultdict
from itertools import combinations

version = '0.3.0dev'

############################
# processing of user input #
############################

class parser():
    def __init__(self):
        argparser = argparse.ArgumentParser(
            # description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
            usage='''smudgeplot <task> [options] \n
tasks: cutoff    Calculate meaningful values for lower kmer histogram cutoff.
       hetmers   Calculate unique kmer pairs from a FastK k-mer database.
       plot      Generate 2d histogram; infere ploidy and plot a smudgeplot.\n\n''')
        # removing this for now;
        #        extract   Extract kmer pairs within specified coverage sum and minor covrage ratio ranges
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

    def hetmers(self):
        '''
        Calculate unique kmer pairs from a Jellyfish or KMC dump file.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot hetkmers',
            description='Calculate unique kmer pairs from FastK k-mer database.')
        argparser.add_argument('infile', nargs='?', help='Input FastK database (.ktab) file.')
        argparser.add_argument('-L', help='Count threshold below which k-mers are considered erroneous', type=int)
        argparser.add_argument('-t', help='Number of threads (default 4)', type=int, default=4)
        argparser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='kmerpairs')
        argparser.add_argument('-tmp', help='Directory where all temporary files will be stored (default /tmp).', default='.')
        argparser.add_argument('--verbose', action="store_true", default = False, help='verbose mode')
        self.arguments = argparser.parse_args(sys.argv[2:])

    def plot(self):
        '''
        Generate 2d histogram; infer ploidy and plot a smudgeplot.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot plot', description='Generate 2d histogram for smudgeplot')
        argparser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages and frequencies (default \"coverages_2.smu\")."')
        argparser.add_argument('-o', help='The pattern used to name the output (smudgeplot).', default='smudgeplot')
        argparser.add_argument('-q', help='Remove kmer pairs with coverage over the specified quantile; (default none).', type=float, default=1)
        argparser.add_argument('-L', help='The lower boundary used when dumping kmers (default min(total_pair_cov) / 2).', type=int, default=0)
        argparser.add_argument('-c', '-cov_filter', help='Filter pairs with one of them having coverage bellow specified threshold (default 0; disables parameter L)', type=int, default=0)
        argparser.add_argument('-n', help='The expected haploid coverage (default estimated from data).', type=float, default=0)
        argparser.add_argument('-t', '--title', help='name printed at the top of the smudgeplot (default none).', default='')
        argparser.add_argument('-ylim', help='The upper limit for the coverage sum (the y axis)', type = int, default=0)
        # argparser.add_argument('-m', '-method', help='The algorithm for annotation of smudges (default \'local_aggregation\')', default='local_aggregation')
        argparser.add_argument('-nbins', help='The number of nbins used for smudgeplot matrix (nbins x nbins) (default autodetection).', type=int, default=0)
        # argparser.add_argument('-kmer_file', help='Name of the input files containing kmer seuqences (assuming the same order as in the coverage file)', default = "")
        argparser.add_argument('--homozygous', action="store_true", default = False, help="Assume no heterozygosity in the genome - plotting a paralog structure; (default False).")
        argparser.add_argument('--just_plot', action="store_true", default = False, help="Turns off the inference of coverage and annotation of smudges; simply generates smudgeplot. (default False)")
        argparser.add_argument('--alt_plot', action="store_true", default = False, help="Uses a new way to plot smudgeplots using tiling strategy, which is likely to be the default for the Oriel 0.3.0 release (default False)")

        # plotting arugments
        argparser.add_argument('-col_ramp', help='An R palette used for the plot (default "viridis", other sensible options are "magma", "mako" or "grey.colors" - recommended in combination with --invert_cols).', default='viridis')
        argparser.add_argument('--invert_cols', action="store_true", default = False, help="Revert the colour palette (default False).")
        argparser.add_argument('--plot_err_line', action="store_true", default = False, help="Add a line to the plot denoting where the error k-mers and genomic k-mers will likely pair up the most (default False).")
        
        self.arguments = argparser.parse_args(sys.argv[2:])

    def cutoff(self):
        '''
        Calculate meaningful values for lower/upper kmer histogram cutoff.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot cutoff', description='Calculate meaningful values for lower/upper kmer histogram cutoff.')
        argparser.add_argument('infile', type=argparse.FileType('r'), help='Name of the input kmer histogram file (default \"kmer.hist\")."')
        argparser.add_argument('boundary', help='Which bounary to compute L (lower) or U (upper)')
        self.arguments = argparser.parse_args(sys.argv[2:])

    def extract(self):
        '''
        Extract kmer pairs within specified coverage sum and minor covrage ratio ranges.
        '''
        argparser = argparse.ArgumentParser(prog = 'smudgeplot extract', description='Extract kmer pairs within specified coverage sum and minor covrage ratio ranges.')
        argparser.add_argument("-cov", "--coverageFile",required=True, help="coverage file for the kmer pairs")
        argparser.add_argument("-seq", "--seqFile",required=True, help="sequences of the kmer pairs")
        argparser.add_argument("-minc", "--countMin",required=True, help="lower bound of the summed coverage", type=int)
        argparser.add_argument("-maxc", "--countMax",required=True, help="upper bound of the summed coverage", type=int)
        argparser.add_argument("-minr", "--ratioMin",required=True, help="lower bound of minor allele ratio", type=float)
        argparser.add_argument("-maxr", "--ratioMax",required=True, help="upper bound of minor allele ratio", type=float)
        self.arguments = argparser.parse_args(sys.argv[2:])

###############
# task cutoff #
###############

# taken from https://stackoverflow.com/a/29614335
def local_min(ys):
    return [i for i, y in enumerate(ys)
            if ((i == 0) or (ys[i - 1] >= y))
            and ((i == len(ys) - 1) or (y < ys[i+1]))]

def round_up_nice(x):
    digits = ceil(log(x, 10))
    if digits <= 1:
        multiplier = 10 ** (digits - 1)
    else:
        multiplier = 10 ** (digits - 2)
    return(ceil(x / multiplier) * multiplier)

def cutoff(args):
    # kmer_hist = open("data/Scer/kmc_k31.hist","r")
    kmer_hist = args.infile
    hist = [int(line.split()[1]) for line in kmer_hist]
    if args.boundary == "L":
        local_minima = local_min(hist)[0]
        L = max(10, int(round(local_minima * 1.25)))
        sys.stdout.write(str(L))
    else:
        sys.stderr.write('Warning: We discourage using the original hetmer algorithm.\n\tThe updated (recommended) version does not take the argument U\n')
        # take 99.8 quantile of kmers that are more than one in the read set
        number_of_kmers = sum(hist[1:])
        hist_rel_cumsum = [sum(hist[1:i+1]) / number_of_kmers for i in range(1, len(hist))]
        min(range(len(hist_rel_cumsum))) 
        U = round_up_nice(min([i for i, q in enumerate(hist_rel_cumsum) if q > 0.998]))
        sys.stdout.write(str(U))
    sys.stdout.flush()

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

    if _parser.task == "hetmers":
        # PloidyPlot is expected ot be installed in the system as well as the R library supporting it
        args = _parser.arguments
        plot_args = " -o" + str(args.o)
        plot_args += " -e" + str(args.L)
        plot_args += " -T" + str(args.t)
        if args.verbose:
            plot_args += " -v"
        if args.tmp != '.':
            plot_args += " -P" + args.tmp
        plot_args += " " + args.infile

        sys.stderr.write("Calling: PloidyPlot " + plot_args + "\n")
        system("PloidyPlot " + plot_args)

    if _parser.task == "plot":
        # the plotting script is expected ot be installed in the system as well as the R library supporting it
        args = _parser.arguments
        plot_args = "-i \"" + args.infile + "\" -o \"" + args.o + "\""
        if args.q != 1:
            plot_args += " -q " + str(args.q)
        if args.L != 0:
            plot_args += " -L " + str(args.L)
        if args.c != 0:
            plot_args += " -c " + str(args.c)
        if args.n != 0:
            plot_args += " -n " + str(args.n)
        if args.title:
            plot_args += " -t \"" + args.title + "\""
        if args.ylim != 0:
            plot_args += " -ylim " + str(args.ylim)
        if args.col_ramp:
            plot_args += " -col_ramp \"" + args.col_ramp + "\""
        if args.nbins != 0:
            plot_args += " -nbins " + str(args.nbins)
        if args.homozygous:
            plot_args += " --homozygous"
        if args.invert_cols:
            plot_args += " --invert_cols"
        if args.plot_err_line:
            plot_args += " --plot_err_line"
        if args.just_plot:
            plot_args += " --just_plot"
        if args.alt_plot:
            plot_args += " --alt_plot"

        sys.stderr.write("Calling: smudgeplot_plot.R " + plot_args + "\n")
        system("smudgeplot_plot.R " + plot_args)

    if _parser.task == "extract":
        extract_kmer_pairs(_parser.arguments)


    sys.stderr.write("\nDone!\n")
    exit(0)

if __name__=='__main__':
    main()

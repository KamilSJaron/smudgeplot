#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from pandas import read_csv
from pandas import DataFrame
from numpy import arange
from numpy import argmin
from numpy import concatenate
from os import system
from math import log
from math import ceil
from statistics import fmean
from collections import defaultdict
# from matplotlib.pyplot import plot

version = '0.3.0dev'

############################
# processing of user input #
############################

class parser():
    def __init__(self):
        argparser = argparse.ArgumentParser(
            # description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
            usage='''smudgeplot <task> [options] \n
tasks: cutoff           Calculate meaningful values for lower kmer histogram cutoff.
       hetmers          Calculate unique kmer pairs from a FastK k-mer database.
       peak_agregation  Agregates smudges using local agregation algorithm.
       plot             Generate 2d histogram; infere ploidy and plot a smudgeplot.
       all (default)    Runs all the steps ()\n\n''')
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

    def peak_agregation(self):
        '''
        Extract kmer pairs within specified coverage sum and minor covrage ratio ranges.
        '''
        argparser = argparse.ArgumentParser()
        argparser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages and frequencies.')
        argparser.add_argument('-nf', '-noise_filter', help='Do not agregate into smudge k-mer pairs with frequency lower than this parameter', type=int, default=50)
        argparser.add_argument('-d', '-distance', help='Manthattan distance of k-mer pairs that are considered neioboring for the local agregation purposes.', type=int, default=5)
        argparser.add_argument('--mask_errors', help='instead of reporting assignments to individual smudges, just remove all monotonically decreasing points from the error line', action="store_true", default = False)
        self.arguments = argparser.parse_args(sys.argv[2:])

    def all(self):
        argparser = argparse.ArgumentParser()
        argparser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages and frequencies.')
        argparser.add_argument('-o', help='The pattern used to name the output (smudgeplot).', default='smudgeplot')
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

########################
# task peak_agregation #
########################

def load_hetmers(smufile):
    cov_tab = read_csv(smufile, names = ['covB', 'covA', 'freq'], sep='\t')
    cov_tab = cov_tab.sort_values('freq', ascending = False)
    return(cov_tab)

def local_agregation(cov_tab, distance, noise_filter, mask_errors):
    # generate a dictionary that gives us for each combination of coverages a frequency
    cov2freq = defaultdict(int)
    cov2peak = defaultdict(int)

    L = min(cov_tab['covB']) # important only when --mask_errors is on

    ### clustering
    next_peak = 1
    for idx, covB, covA, freq in cov_tab.itertuples():
        cov2freq[(covA, covB)] = freq # a make a frequency dictionary on the fly, because I don't need any value that was not processed yet
        if freq < noise_filter:
            break
        highest_neigbour_coords = (0, 0)
        highest_neigbour_freq = 0
        # for each kmer pair I will retrieve all neibours (Manhattan distance)
        for xA in range(covA - distance,covA + distance + 1):
            # for explored A coverage in neiborhood, we explore all possible B coordinates
            distanceB = distance - abs(covA - xA)
            for xB in range(covB - distanceB,covB + distanceB + 1):
                xB, xA = sorted([xA, xB]) # this is to make sure xB is smaller than xA
                # iterating only though those that were assigned already
                # and recroding only the one with highest frequency
                if cov2peak[(xA, xB)] and cov2freq[(xA, xB)] > highest_neigbour_freq:
                    highest_neigbour_coords = (xA, xB)
                    highest_neigbour_freq = cov2freq[(xA, xB)]
        if highest_neigbour_freq > 0:
            cov2peak[(covA, covB)] = cov2peak[(highest_neigbour_coords)]
        else:
            # print("new peak:", (covA, covB))
            if mask_errors:
                if covB < L + distance:
                    cov2peak[(covA, covB)] = 1 # error line
                else:
                    cov2peak[(covA, covB)] = 0 # central smudges
            else:
                cov2peak[(covA, covB)] = next_peak # if I want to keep info about all locally agregated smudges
                next_peak += 1
    return(cov2peak)

def peak_agregation(args):
    ### load data
    cov_tab = load_hetmers(args.infile)

    cov2peak = local_agregation(cov_tab, args.d, args.nf, mask_errors = False)

    cov_tab = cov_tab.sort_values(['covA', 'covB'], ascending = True)
    for idx, covB, covA, freq in cov_tab.itertuples():
        peak = cov2peak[(covA, covB)]
        sys.stdout.write(f"{covB}\t{covA}\t{freq}\t{peak}\n")
    sys.stdout.flush()

def get_smudge_container(cov_tab, cov, smudge_filter):
    smudge_container = dict()
    genomic_cov_tab = cov_tab[cov_tab['peak'] == 0] # this removed all the marked errors
    total_kmer_pairs = sum(genomic_cov_tab['freq'])

    for Bs in range(1,9):
        min_cov = 0 if Bs == 1 else cov * (Bs - 0.5)
        max_cov = cov * (Bs + 0.5)
        cov_tab_isoB = genomic_cov_tab.loc[(genomic_cov_tab["covB"] > min_cov) & (genomic_cov_tab["covB"] < max_cov)] #  

        for As in range(Bs,(17 - Bs)):
            min_cov = 0 if As == 1 else cov * (As - 0.5)
            max_cov = cov * (As + 0.5)
            cov_tab_iso_smudge = cov_tab_isoB.loc[(cov_tab_isoB["covA"] > min_cov) & (cov_tab_isoB["covA"] < max_cov)]
            if sum(cov_tab_iso_smudge['freq']) / total_kmer_pairs > smudge_filter:
                # sys.stderr.write(f"{As}A{Bs}B: {sum(cov_tab_iso_smudge['freq']) / total_kmer_pairs}\n")
                smudge_container["A" * As + "B" * Bs] = cov_tab_iso_smudge
    return(smudge_container)

def get_centrality(smudge_container, cov):
    centralities = list()
    freqs = list()
    for smudge in smudge_container.keys():
        As = smudge.count('A')
        Bs = smudge.count('B')
        smudge_tab = smudge_container[smudge]
        freqs.append(sum(smudge_tab['freq']))
        center = smudge_tab.loc[smudge_tab['freq'].idxmax()]
        ## emprical to edge
        # distA = min([abs(smudge_tab['covA'].max() - center['covA']), abs(center['covA'] - smudge_tab['covA'].min())])
        # distB = min([abs(smudge_tab['covB'].max() - center['covB']), abs(center['covB'] - smudge_tab['covB'].min())])
        ## theoretical to edge
        # distA = min(abs(center['covA'] - (cov * (As - 0.5))), abs((cov * (As + 0.5)) - center['covA']))
        # distB = min(abs(center['covB'] - (cov * (Bs - 0.5))), abs((cov * (Bs + 0.5)) - center['covB']))
        ## theoretical relative distance to the center
        distA = abs((center['covA'] - (cov * As)) / cov)
        distB = abs((center['covB'] - (cov * Bs)) / cov)

        # sys.stderr.write(f"Processing: {As}A{Bs}B; with center: {distA}, {distB}\n")
        centrality = distA + distB
        centralities.append(centrality)

    return(fmean(centralities, weights=freqs))

def test_coverage_range(cov_tab, min_c, max_c, smudge_size_cutoff = 0.02):
    # covs_to_test = range(min_c, max_c)
    covs_to_test = arange(min_c + 0.05, max_c + 0.05, 2)
    cov_centralities = list()
    for cov in covs_to_test:
        smudge_container = get_smudge_container(cov_tab, cov, smudge_size_cutoff)
        cov_centralities.append(get_centrality(smudge_container, cov))

    best_coverage = covs_to_test[argmin(cov_centralities)]

    tenths_to_test = arange(best_coverage - 1.9, best_coverage + 1.9, 0.2)
    tenths_centralities = list()
    for cov in tenths_to_test:
        smudge_container = get_smudge_container(cov_tab, cov, smudge_size_cutoff)
        tenths_centralities.append(get_centrality(smudge_container, cov))

    best_tenth = tenths_to_test[argmin(tenths_centralities)]
    sys.stderr.write(f"Best coverage to precsion of one tenth: {round(best_tenth, 2)}\n")  

    hundredths_to_test = arange(best_tenth - 0.19, best_tenth + 0.19, 0.01)
    hundredths_centralities = list()
    for cov in hundredths_to_test:
        smudge_container = get_smudge_container(cov_tab, cov, smudge_size_cutoff)
        hundredths_centralities.append(get_centrality(smudge_container, cov))

    final_cov = hundredths_to_test[argmin(hundredths_centralities)]
    sys.stderr.write(f"Best coverage to precision of one hundreth: {round(final_cov, 3)}\n")  

    all_coverages = concatenate((covs_to_test, tenths_to_test, hundredths_to_test))
    all_centralities = concatenate((cov_centralities, tenths_centralities, hundredths_centralities))

    return(DataFrame({'coverage': all_coverages, 'centrality': all_centralities}))

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

    if _parser.task == "peak_agregation":
        peak_agregation(_parser.arguments)

    if _parser.task == "all":
        args = _parser.arguments

        sys.stderr.write("\nLoading data\n")
        cov_tab = load_hetmers(args.infile)
        # cov_tab = load_hetmers("data/dicots/smu_files/daAchMill1.k31_ploidy.smu.txt")

        sys.stderr.write("\nMasking errors using local agregation algorithm\n")
        cov2peak = local_agregation(cov_tab, distance = 1, noise_filter = 1000, mask_errors = True)
        cov_tab['peak'] = [cov2peak[(covA, covB)] for idx, covB, covA, freq in cov_tab.itertuples()]

        cov_tab = cov_tab.sort_values(['covA', 'covB'], ascending = True)
        total_kmers = sum(cov_tab[cov_tab['peak'] == 0]['freq'])

        with open(args.o + "_masked_errors_smu.txt", 'w') as error_annotated_smu:
            error_annotated_smu.write("covB\tcovA\tfreq\tis_error\n")
            for idx, covB, covA, freq, is_error in cov_tab.itertuples():
                error_annotated_smu.write(f"{covB}\t{covA}\t{freq}\t{is_error}\n") # might not be needed

        sys.stderr.write("\nInfering 1n coverage using grid algorihm\n")

        smudge_size_cutoff = 0.01 # this is % of all k-mer pairs smudge needs to have to be considered a valid smudge
        centralities = test_coverage_range(cov_tab, 10, 60, smudge_size_cutoff)
        np.savetxt(args.o + "_centralities.txt", np.around(centralities, decimals=6), fmt="%.4f")
        # plot(centralities['coverage'], centralities['coverage'])

        cov = centralities['coverage'][argmin(centralities['centrality'])]
        # sys.stderr.write(f"\nInferred coverage: {cov}\n")
        final_smudges = get_smudge_container(cov_tab, cov, smudge_size_cutoff)
        
        annotated_smudges = list(final_smudges.keys())
        smudge_sizes = [round(sum(final_smudges[smudge]['freq']) / total_kmers, 4) for smudge in annotated_smudges]

        # saving smudge sizes
        # smudge_table = DataFrame({'smudge': annotated_smudges, 'size': smudge_sizes})
        # np.savetxt(args.o + "_smudge_sizes.txt", smudge_table)

        sys.stderr.write("\nPlotting\n")

        system("centrality_plot.R " + args.o + "_centralities.txt")
        # Rscript playground/alternative_fitting/alternative_plotting_testing.R -i data/dicots/peak_agregation/$ToLID.cov_tab_peaks -o data/dicots/peak_agregation/$ToLID
        system("smudgeplot_plot.R -i" + args.infile + " -o " + args.o + " -n " + str(cov) + " -ylim 300 --alt_plot") 

    sys.stderr.write("\nDone!\n")
    exit(0)

if __name__=='__main__':
    main()

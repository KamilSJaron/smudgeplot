#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from pandas import read_csv # type: ignore
from pandas import DataFrame # type: ignore
from numpy import arange
from numpy import argmin
from numpy import concatenate
from os import system
from math import log
from math import ceil
from statistics import fmean
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib.pyplot import plot

version = '0.4.0dev'

############################
# processing of user input #
############################

class Parser:
	def __init__(self):
		self.arguments=None
		argparser = argparse.ArgumentParser(
			# description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
			usage="""
			smudgeplot <task> [options] \n
			tasks: cutoff           Calculate meaningful values for lower kmer histogram cutoff.
				   hetmers          Calculate unique kmer pairs from a FastK k-mer database.
				   peak_aggregation  Agregates smudges using local aggregation algorithm.
				   plot             Generate 2d histogram; infere ploidy and plot a smudgeplot.
				   all              Runs all the steps (with default options)\n\n
				   """
		)
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
		"""
		Calculate unique kmer pairs from a Jellyfish or KMC dump file.
		"""
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
		"""
		Generate 2d histogram; infer ploidy and plot a smudgeplot.
		"""
		argparser = argparse.ArgumentParser(prog = 'smudgeplot plot', description='Generate 2d histogram for smudgeplot')
		argparser.add_argument('infile', help='name of the input tsv file with covarages and frequencies')
		argparser.add_argument('smudgefile', help='name of the input tsv file with sizes of individual smudges')
		argparser.add_argument('n', help='The expected haploid coverage.', type=float)
		argparser.add_argument('-o', help='The pattern used to name the output (smudgeplot).', default='smudgeplot')
		
		argparser = self.add_plotting_arguments(argparser)

		self.arguments = argparser.parse_args(sys.argv[2:])

	def cutoff(self):
		"""
		Calculate meaningful values for lower/upper kmer histogram cutoff.
		"""
		argparser = argparse.ArgumentParser(prog = 'smudgeplot cutoff', description='Calculate meaningful values for lower/upper kmer histogram cutoff.')
		argparser.add_argument('infile', type=argparse.FileType('r'), help='Name of the input kmer histogram file (default \"kmer.hist\")."')
		argparser.add_argument('boundary', help='Which bounary to compute L (lower) or U (upper)')
		self.arguments = argparser.parse_args(sys.argv[2:])

	def peak_aggregation(self):
		"""
		Extract kmer pairs within specified coverage sum and minor covrage ratio ranges.
		"""
		argparser = argparse.ArgumentParser()
		argparser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages and frequencies.')
		argparser.add_argument('-nf', '-noise_filter', help='Do not agregate into smudge k-mer pairs with frequency lower than this parameter', type=int, default=50)
		argparser.add_argument('-d', '-distance', help='Manthattan distance of k-mer pairs that are considered neioboring for the local aggregation purposes.', type=int, default=5)
		argparser.add_argument('--mask_errors', help='instead of reporting assignments to individual smudges, just remove all monotonically decreasing points from the error line', action="store_true", default = False)
		self.arguments = argparser.parse_args(sys.argv[2:])

	def all(self):
		argparser = argparse.ArgumentParser()
		argparser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages and frequencies.')
		argparser.add_argument('-o', help='The pattern used to name the output (smudgeplot).', default='smudgeplot')
		argparser.add_argument('-cov_min', help='Minimal coverage to explore (default 6)', default=6)
		argparser.add_argument('-cov_max', help='Maximal coverage to explore (default 50)', default=60)

		argparser = self.add_plotting_arguments(argparser)

		self.arguments = argparser.parse_args(sys.argv[2:])
	
	def add_plotting_arguments(self, argparser):
		argparser.add_argument('-c', '-cov_filter', help='Filter pairs with one of them having coverage bellow specified threshold (default 0; disables parameter L)', type=int, default=0)
		argparser.add_argument('-t', '--title', help='name printed at the top of the smudgeplot (default none).', default='')
		argparser.add_argument('-ylim', help='The upper limit for the coverage sum (the y axis)', type = int, default=0)
		argparser.add_argument('-col_ramp', help='An R palette used for the plot (default "viridis", other sensible options are "magma", "mako" or "grey.colors" - recommended in combination with --invert_cols).', default='viridis')
		argparser.add_argument('--invert_cols', action="store_true", default = False, help="Revert the colour palette (default False).")
		return argparser
	
	def format_arguments_for_R_plotting(self):
		plot_args = ""
		if self.arguments.c != 0:
			plot_args += " -c " + str(self.arguments.c)
		if self.arguments.title:
			plot_args += " -t \"" + self.arguments.title + "\""
		if self.arguments.ylim != 0:
			plot_args += " -ylim " + str(self.arguments.ylim)
		if self.arguments.col_ramp:
			plot_args += " -col_ramp \"" + self.arguments.col_ramp + "\""
		if self.arguments.invert_cols:
			plot_args += " --invert_cols"
		return plot_args

class SmudgeDataObj(object):
	def __init__(
		self,
		input_file_path
	):
		self.input_file_path = input_file_path
		self.cov_tab=None
		self.cov2peak=None
		self.final_smudges=None
		self.total_kmers=None
		self.total_genomic_kmers=None
		self.total_error_kmers=None
		self.error_fraction=None
		self.centralities=None
		self.cov=None
		self.smudge_tab=None
		self.ylim=None
		self.xlim=None
		self.error_string=None
		self.cov_string=None
		self.fig_title=None
		self.linear_plot_file=None
		self.log_plot_file=None



	def load_hetmers(self):
		self.cov_tab = read_csv(self.input_file_path, names = ['covB', 'covA', 'freq'], sep='\t')
		self.cov_tab.sort_values('freq', ascending = False, inplace=True)

	def local_aggregation(self, distance, noise_filter, mask_errors):
		# generate a dictionary that gives us for each combination of coverages a frequency
		cov2freq = defaultdict(int)
		cov2peak = defaultdict(int)

		L = min(self.cov_tab['covB']) #important only when --mask_errors is on

		### clustering
		next_peak = 1
		for idx, covB, covA, freq in self.cov_tab.itertuples():
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
					# and recording only the one with highest frequency
					if cov2peak[(xA, xB)] and cov2freq[(xA, xB)] > highest_neigbour_freq:
						highest_neigbour_coords = (xA, xB)
						highest_neigbour_freq = cov2freq[(xA, xB)]
			if highest_neigbour_freq:# > 0:
				cov2peak[covA, covB] = cov2peak[highest_neigbour_coords]
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
		self.cov2peak = cov2peak

	def peak_aggregation(self):
		self.cov_tab['peak'] = [self.cov2peak[(covA, covB)] for idx, covB, covA, freq in self.cov_tab.itertuples()]
		self.cov_tab.sort_values(['covA', 'covB'], ascending = True, inplace=True)    

	def write_peaks(self):
		self.peak_aggregation()
		for idx, covB, covA, freq, peak in self.cov_tab.itertuples():
			sys.stdout.write(f"{covB}\t{covA}\t{freq}\t{peak}\n")
		sys.stdout.flush()

	def count_kmers(self):
		self.peak_aggregation()
		
		self.total_kmers = sum(self.cov_tab['freq'])
		self.total_genomic_kmers = sum(self.cov_tab.loc[self.cov_tab['peak'] == 0]['freq'])
		self.total_error_kmers = sum(self.cov_tab.loc[self.cov_tab['peak'] == 1]['freq'])
		self.error_fraction = self.total_error_kmers / self.total_kmers

	def get_best_coverage(self, cov_list, smudge_size_cutoff = 0.02, centralities=None, last_check=False):
		if centralities is None:
			centralities = []

		if last_check:
			to_test = [cov_list[-1]]
		else:
			to_test = cov_list 
		
		for cov in to_test:
			smudge_container = self.get_smudge_container(cov, smudge_size_cutoff)
			centralities.append(get_centrality(smudge_container, cov))
		return cov_list[argmin(centralities)], centralities

	def test_coverages(self, min_c, max_c, smudge_size_cutoff = 0.02):
		grid_params = [(0.05, 0.05, 2), (-1.9, 1.9, 0.2), (-0.19, 0.19, 0.01)]
		results = []

		for i, params in enumerate(grid_params):
			cov_list = arange(min_c + params[0], max_c + params[1], params[2])
			best_cov, centralities = self.get_best_coverage(cov_list, smudge_size_cutoff)
			
			results.append({
				'covs': cov_list, 
				'centralities':centralities, 
				'best_cov':best_cov
			})

			min_c, max_c = best_cov, best_cov
			
			if i > 0:
				sys.stderr.write(f"Best coverage to precision of 1/{10**i}: {round(best_cov, 2)}\n")  
			
		# just to be sure
		results[-1]['covs'] = np.append(results[-1]['covs'], results[-1]['best_cov']/2)
		best_cov, centralities = self.get_best_coverage(
			cov_list = results[-1]['covs'], 
			smudge_size_cutoff = smudge_size_cutoff,
			centralities=results[-1]['centralities'],
			last_check=True 
		)

		sys.stderr.write(f"Best coverage to precision of 1/{10**i} (just to be sure): {round(best_cov, 3)}\n")  

		self.centralities = DataFrame({
			'coverage': concatenate([result['covs'] for result in results]), 
			'centrality': concatenate([result['centralities'] for result in results])}
		)

	def infer_coverage(self, limit=0.7):
		if self.error_fraction < limit:
			self.cov = self.centralities['coverage'][argmin(self.centralities['centrality'])]
		else:
			self.cov = 0

	def describe_smudges(self):
		annotated_smudges = list(self.final_smudges.keys())
		smudge_sizes = [round(sum(self.final_smudges[smudge]['freq']) / self.total_genomic_kmers, 4) for smudge in annotated_smudges]
		self.smudge_tab = DataFrame({'structure': annotated_smudges, 'size': smudge_sizes})
	
	def get_smudge_container(self, cov, smudge_filter):
		smudge_container = {}

		for Bs in range(1,9):
			min_cov, max_cov = get_cov_limits(Bs, cov)

			cov_tab_isoB = self.cov_tab.loc[ \
			(self.cov_tab['peak'] == 0) & \
			(self.cov_tab["covB"] > min_cov) & \
			(self.cov_tab["covB"] < max_cov)
			]  

			for As in range(Bs,(17 - Bs)):
				min_cov, max_cov = get_cov_limits(As, cov)

				cov_tab_iso_smudge = cov_tab_isoB.loc[(cov_tab_isoB["covA"] > min_cov) & (cov_tab_isoB["covA"] < max_cov)]
				if sum(cov_tab_iso_smudge['freq']) / self.total_genomic_kmers > smudge_filter:
					smudge_container["A" * As + "B" * Bs] = cov_tab_iso_smudge
					# sys.stderr.write(f"{As}A{Bs}B: {sum(cov_tab_iso_smudge['freq']) / total_kmer_pairs}\n")

		return smudge_container

	def calc_cov_columns(self):
		self.cov_tab['total_pair_cov'] = self.cov_tab[['covA','covB']].sum(axis=1)
		self.cov_tab['minor_variant_rel_cov'] = self.cov_tab['covB']/self.cov_tab['total_pair_cov']

	def filter_cov_quant(self, cov_filter=None, quant_filter=None):
		if cov_filter:
			self.cov_tab = self.cov_tab.loc[(self.cov_tab['covA'] >= cov_filter) & (self.cov_tab['covB'] >= cov_filter)]
		if quant_filter:
			# 0 < int(quant_filter) < 100
			upper_quantile = np.percentile(a=self.cov_tab['total_pair_cov'], q=quant_filter, weights = self.cov_tab['freq'], method='inverted_cdf')
			self.cov_tab = self.cov_tab.loc[self.cov_tab['total_pair_cov'] < upper_quantile]

	def get_ax_lims(self, upper_ylim=None):
		if self.cov == np.percentile(a=self.cov_tab['total_pair_cov'], q=95, weights = self.cov_tab['freq'], method='inverted_cdf'):
			self.ylim = [min(self.cov_tab['total_pair_cov']), max(self.cov_tab['total_pair_cov'])]
		else: 
			self.ylim = [min(self.cov_tab['total_pair_cov']) - 1, # or 0?
				  min(max(100, 10*self.cov), max(self.cov_tab['total_pair_cov']))]
		if upper_ylim:
			self.ylim[1] = upper_ylim
		self.xlim = [0, 0.5]

	def def_strings(self, title=None, output='smudgeplot'):
			error_fraction = np.sum(self.cov_tab[self.cov_tab['peak']==1]['freq']) / np.sum(self.cov_tab['freq']) * 100
			self.error_string = f"err = {round(error_fraction, 1)} %"
			self.cov_string = f"1n = {self.cov}"
			if title:
				self.fig_title = str(title)
			else:
				self.fig_title = '.'.join(self.input_file_path.split('/')[-1].split('.')[0:2])
			self.linear_plot_file = output+'_smudgeplot_py.pdf'
			self.log_plot_file = output+'_smudgeplot_log10_py.pdf'

	def smudgeplot(self, log=False):
	#probably needs to be a function not a method because of plot task, but is ok for now

		fig, axs = plt.subplots(nrows=2, ncols=2, width_ratios = [3, 1], height_ratios = [1, 3], figsize=(20,20))
		plt.subplots_adjust(wspace=0.05, hspace=0.05)
		fontsize=32

		main_ax = axs[1][0]
		legend_ax = axs[0][1]
		size_ax = axs[1][1]
		top_ax = axs[0][0]

		legend_ax.axis('off')
		size_ax.axis('off')
		top_ax.axis('off')

		if log:
			colour_ramp = get_col_ramp(delay=16)
			outfile = self.log_plot_file
		else:
			colour_ramp = get_col_ramp()
			outfile = self.linear_plot_file
		
		self.plot_alt(colour_ramp, log=log, fontsize=fontsize, ax=main_ax)

		if self.cov>0:
			self.plot_expected_haplotype_structure(main_ax, adjust=True, xmax = 0.49)

		plot_legend(ax=legend_ax, kmer_max=max(self.cov_tab['freq']), colour_ramp=colour_ramp, log=log)

		self.plot_smudge_sizes(ax=size_ax)

		top_ax.set_title(self.fig_title, fontsize=42, loc ='left', y=1.0, pad=-14, weight='bold')

		fig.savefig(outfile, dpi = 200)

	def plot_alt(self, colour_ramp, ax, log=False, fontsize=14):
		mask = (self.cov_tab['covA'] == self.cov_tab['covB'])
		self.cov_tab.loc[mask, 'freq'] = self.cov_tab[mask]['freq']*2

		if log:
			self.cov_tab['freq'] = np.log10(self.cov_tab['freq'])

		# might be a quicker way to do this
		self.cov_tab['col'] = [str(colour_ramp[int(i)]) for i in round((len(colour_ramp)-1) * self.cov_tab['freq'] / max(self.cov_tab['freq']))]

		ax.plot()
		ax.set_xlim(self.xlim)
		ax.set_ylim(self.ylim)
		ax.set_xlabel('Normalized minor kmer coverage: B / (A + B)', fontsize=fontsize)
		ax.set_ylabel('Total coverage of the kmer pair: A + B', fontsize=fontsize)
		ax.tick_params(axis='both', labelsize=20)
		ax.spines[['right', 'top']].set_visible(False)

		min_cov_to_plot = max(self.ylim[0],min(self.cov_tab['total_pair_cov']))

		for cov in np.arange(min_cov_to_plot, self.ylim[1]):
			self.plot_one_coverage(cov, ax)

	def plot_one_coverage(self, cov, ax):
		cov_row_to_plot = self.cov_tab[self.cov_tab['total_pair_cov'] == cov]
		width = 1/(2*cov)
		lefts = (cov_row_to_plot['minor_variant_rel_cov'] - width).to_numpy()
		rights = (cov_row_to_plot['minor_variant_rel_cov'].apply(lambda x: min(0.5, x+width))).to_numpy()
		cols = cov_row_to_plot['col'].to_numpy()
		for left, right, col in zip(lefts, rights, cols):
			plot_one_box(left, right, cov, col, ax)

	def plot_expected_haplotype_structure(self, ax, adjust=False, xmax=0.49):
		# find the slice warnings
		self.smudge_tab.loc[:,'ploidy'] = self.smudge_tab['structure'].str.len()
		self.smudge_tab = self.smudge_tab.loc[self.smudge_tab['size'] > 0.05]
		self.smudge_tab.loc[:,'corrected_minor_variant_cov'] = self.smudge_tab['structure'].str.count('B')/self.smudge_tab['ploidy']
		self.smudge_tab.loc[:,'label'] = reduce_structure_representation(self.smudge_tab['structure'])
		bordercases = np.array(self.smudge_tab['corrected_minor_variant_cov']==0.5)
		for index, row in self.smudge_tab.iterrows():
			if bordercases[index] & adjust:
				x = (xmax + 0.49) / 2
			else:
				x = row['corrected_minor_variant_cov']
			y = row['ploidy']*self.cov
			ax.text(x, y, row['label'], fontsize=28)

	def plot_smudge_sizes(self, ax):
		ax.plot()
		ax.set_title('')
		if self.cov > 0:
			labels = [smudge+"    "+str(size) for smudge, size in 
				zip(
					reduce_structure_representation(self.smudge_tab['structure']).to_list(), 
					round(self.smudge_tab['size'], 2)
					)
				]
			handles = [mpl.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',
                      visible=False) for _ in labels]
			ax.legend(loc='upper left', labels=labels, handles = handles, prop={'size': 18}, frameon=False)

		
		else:
			ax.legend(loc='lower left', handles=[self.error_string],labels=[None], prop={'size': 18}, frameon=False)

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
	return ceil(x / multiplier) * multiplier

def cutoff(kmer_hist, boundary):
	# kmer_hist = open("data/Scer/kmc_k31.hist","r")
	hist = [int(line.split()[1]) for line in kmer_hist]
	if boundary == "L":
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

#################
# general funcs #
#################

def get_centrality(smudge_container, cov):
	centralities = []
	freqs = []
	for smudge, smudge_tab in smudge_container.items():
		As = smudge.count('A')
		Bs = smudge.count('B')
		kmer_in_the_smudge = sum(smudge_tab['freq'])
		freqs.append(kmer_in_the_smudge)

		# center as a mode 
		center = smudge_tab.loc[smudge_tab['freq'].idxmax()]
		center_A = center['covA']
		center_B = center['covB']
		
		# center as a mean
		# center_A = sum((smudge_tab['freq'] * smudge_tab['covA'])) / kmer_in_the_smudge
		# center_B = sum((smudge_tab['freq'] * smudge_tab['covB'])) / kmer_in_the_smudge

		## emprical to edge
		# distA = min([abs(smudge_tab['covA'].max() - center['covA']), abs(center['covA'] - smudge_tab['covA'].min())])
		# distB = min([abs(smudge_tab['covB'].max() - center['covB']), abs(center['covB'] - smudge_tab['covB'].min())])
		
		## theoretical to edge
		# distA = min(abs(center['covA'] - (cov * (As - 0.5))), abs((cov * (As + 0.5)) - center['covA']))
		# distB = min(abs(center['covB'] - (cov * (Bs - 0.5))), abs((cov * (Bs + 0.5)) - center['covB']))
		
		## theoretical relative distance to the center
		distA = abs((center_A - (cov * As)) / cov)
		distB = abs((center_B - (cov * Bs)) / cov)

		# sys.stderr.write(f"Processing: {As}A{Bs}B; with center: {distA}, {distB}\n")
		centrality = distA + distB
		centralities.append(centrality)

	if len(centralities) == 0:
		return 1
	return fmean(centralities, weights=freqs)

def get_cov_limits(Xs, cov):
	min_cov = 0 if Xs == 1 else cov * (Xs - 0.5)
	max_cov = cov * (Xs + 0.5)
	return min_cov, max_cov

def smudgeplot_plot_R(plot_args):
	sys.stderr.write("Calling: smudgeplot_plot.R " + plot_args + "\n")
	system("smudgeplot_plot.R " + plot_args)

def smudgeplot_plot_py(dataObj, cov_filter=None, quant_filter=None, upper_ylim=None, output='smudgeplot'):
	# preplottingprep
	dataObj.calc_cov_columns()
	dataObj.filter_cov_quant(cov_filter=cov_filter, quant_filter=quant_filter)

	# plotting
	dataObj.get_ax_lims(upper_ylim=upper_ylim)
	dataObj.def_strings(output=output)
	dataObj.smudgeplot(log=False)
	#log plot and legend is currently incurrent
	dataObj.smudgeplot(log=True)

def get_col_ramp(col_ramp='viridis', delay=0, invert_cols=False):
	#is delay incorporated correctly?
	if invert_cols:
		col_ramp+='_r' 
	cmap = plt.get_cmap(col_ramp, 32-int(delay))
	return np.array([mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)])

def reduce_structure_representation(smudge_labels):
	structures_to_adjust = (smudge_labels.str.len() > 4)
	if not any(structures_to_adjust):
		return structures_to_adjust
	else:
		As = smudge_labels[structures_to_adjust].str.count('A').map(str)
		Bs = smudge_labels[structures_to_adjust].str.count('B').map(str)
		new_labels = smudge_labels.copy(deep=True)
		new_labels[structures_to_adjust] = As + 'A' + Bs + 'B'
		return new_labels

def plot_one_box(left, right, cov, colour, ax):
	width = float(right)- float(left)
	rect = mpl.patches.Rectangle((float(left), cov-0.5), width, 1, linewidth=1, edgecolor=colour, facecolor=colour)
	ax.add_patch(rect)

def plot_legend(ax, kmer_max, colour_ramp, log=False):
	title = 'kmers pairs'
	if log:
		title = 'log ' + title
	ax.set_title(title, fontsize=28, weight='bold')
	for i, colour in enumerate(colour_ramp):
		rect = mpl.patches.Rectangle((0, ((i + 1) - 0.01) / 33), 0.5, ((i + 1) + 0.99) / 33, linewidth=1,
									 edgecolor=colour, facecolor=colour)
		ax.add_patch(rect)

	for i in range(6):
		if log:
			ax.text(0.75, i / 6, str(rounding(10 ** (np.log10(kmer_max) * i / 6))), fontsize=20)
		else:
			ax.text(0.75, i / 6, str(rounding(kmer_max * i / 6)), fontsize=20)

def rounding(number):
	if number > 1000:
		return round(number / 1000) * 1000
	elif number > 100:
		return round(number / 100) * 100
	else:
		return round(number / 10) * 10

def fin():
	sys.stderr.write("\nDone!\n")
	exit(0)

#####################
# the script itself #
#####################

def main():
	_parser = Parser()

	sys.stderr.write('Running smudgeplot v' + version + "\n")
	if _parser.task == "version":
		exit(0)

	sys.stderr.write('Task: ' + _parser.task + "\n")

	args = _parser.arguments

	if _parser.task == "cutoff":
		cutoff(args.infile, args.boundary)
		fin()

	if _parser.task == "hetmers":
		# PloidyPlot is expected ot be installed in the system as well as the R library supporting it
		plot_args = " -o" + str(args.o)
		plot_args += " -e" + str(args.L)
		plot_args += " -T" + str(args.t)
		if args.verbose:
			plot_args += " -v"
		if args.tmp != '.':
			plot_args += " -P" + args.tmp
		plot_args += " " + args.infile

		sys.stderr.write("Calling: hetmers (PloidyPlot kmer pair search) " + plot_args + "\n")
		system("hetmers " + plot_args)

		fin()

	if _parser.task == "plot":
		# the plotting script is expected ot be installed in the system as well as the R library supporting it
		
		plot_args = (f'-i "{args.infile}" -s "{args.smudgefile}" -n {args.n} -o "{args.o}"'
					 + _parser.format_arguments_for_R_plotting())
		
		smudgeplot_plot_R(plot_args)

		fin()

	sys.stderr.write("\nLoading data\n")
	SmudgeData = SmudgeDataObj(args.infile)
	SmudgeData.load_hetmers()

	if _parser.task == "peak_aggregation":

		SmudgeData.local_aggregation(distance = args.d, noise_filter = args.nf, mask_errors = False)
		SmudgeData.write_peaks()

	if _parser.task == "all":
		sys.stderr.write("\nMasking errors using local aggregation algorithm\n")
		SmudgeData.local_aggregation(distance = 1, noise_filter = 1000, mask_errors = True)        
		
		SmudgeData.count_kmers()
		sys.stderr.write(f"\t\
			Total kmers: {SmudgeData.total_kmers}\n\t \
			Genomic kmers: {SmudgeData.total_genomic_kmers}\n\t \
			Sequencing errors: {SmudgeData.total_error_kmers}\n\t \
			Fraction or errors: {round(SmudgeData.total_error_kmers/SmudgeData.total_kmers, 3)}"
		)

		with open(args.o + "_masked_errors_smu.txt", 'w') as error_annotated_smu:
			error_annotated_smu.write("covB\tcovA\tfreq\tis_error\n")
			for idx, covB, covA, freq, is_error in SmudgeData.cov_tab.itertuples():
				error_annotated_smu.write(f"{covB}\t{covA}\t{freq}\t{is_error}\n") # might not be needed

		sys.stderr.write("\nInferring 1n coverage using grid algorithm\n")
		smudge_size_cutoff = 0.001 # this is % of all k-mer pairs smudge needs to have to be considered a valid smudge
		SmudgeData.test_coverages(args.cov_min, args.cov_max, smudge_size_cutoff)
		np.savetxt(args.o + "_centralities.txt", 
			np.around(SmudgeData.centralities, decimals=6), 
			fmt="%.4f", 
			delimiter = '\t'
		)
		SmudgeData.infer_coverage(limit=0.7)
		sys.stderr.write(f"\nInferred coverage: {round(SmudgeData.cov, 3)}\n")
		# plot(SmudgeData.centralities['coverage'], SmudgeData.centralities['coverage'])

		SmudgeData.final_smudges = SmudgeData.get_smudge_container(
			SmudgeData.cov,
			smudge_size_cutoff
		)
		
		SmudgeData.describe_smudges()
		sys.stderr.write(f'Detected smudges / sizes ({args.o} + "_smudge_sizes.txt):"\n')
		sys.stderr.write('\t' + str(SmudgeData.smudge_tab['structure'].to_list()) + '\n')
		sys.stderr.write('\t' + str(SmudgeData.smudge_tab['size'].to_list()) + '\n')
		np.savetxt(args.o + "_smudge_sizes.txt", 
			SmudgeData.smudge_tab, 
			fmt='%s', 
			delimiter = '\t'
		)

		sys.stderr.write("\nPlotting\n")
		system("centrality_plot.R " + args.o + "_centralities.txt")
		
		# Rscript playground/alternative_fitting/alternative_plotting_testing.R -i data/dicots/peak_aggregation/$ToLID.cov_tab_peaks -o data/dicots/peak_aggregation/$ToLID
		
		plot_args = (f' -i "{args.o}_masked_errors_smu.txt" -s "{args.o}_smudge_sizes.txt" -n {round(SmudgeData.cov, 3)} -o "{args.o}"'
					 + _parser.format_arguments_for_R_plotting())
		
		#smudgeplot_plot_R(plot_args)
		
		smudgeplot_plot_py(SmudgeData)

	fin()

if __name__=='__main__':
	main()
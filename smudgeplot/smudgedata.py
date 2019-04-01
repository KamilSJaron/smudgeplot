import numpy as np
import matplotlib as mpl
mpl.use('Agg') # https://stackoverflow.com/a/33873802/2962344
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from os.path import isfile
from collections import defaultdict
from bisect import bisect_left
import logging

# smudgedata - a class for easy control over the data related to smudgeplot
#
# data:
#  - args (output of argparse)
#    * infile (open file for reading)
#  - rel_cov (relative coverage of kmer pairs; B / (A + B), where B <= A )
#  - sum_cov (sum of coverages; A + B)
#  - x, y (borders of 2d histogram bins - x corresponds to rel_cov; y to sum_cov)
#  - hist (matrix of kmers in bins)
#  - sorted_hist_indices
#  - smudge_assignment
#  - smudge_centers
#  - brightest_smudge_n

class smudgedata:
    def __init__(self, user_args):
        self.args = user_args
        if self.args.nbins == 0:
            self.nbins = 40
        else:
            self.nbins = self.args.nbins

    def validateArguments(self):
        # check if kmer file is valid
        if self.args.kmer_file and not isfile(self.args.kmer_file):
            # TODO: throw and exception
            logging.error(self.args.kmer_file + ' is not a readable file.')
            exit(1)


    def loadData(self):
        # add a protection against erasing already loaded data
        self.rel_cov = []
        self.sum_cov = []
        for line in self.args.infile:
            c1, c2 = line.split()
            c1 = int(c1)
            c2 = int(c2)
            sum_cov = c1 + c2
            self.rel_cov.append(c1 / sum_cov)
            self.sum_cov.append(sum_cov)
        self.args.infile.close()
        self.rel_cov = np.array(self.rel_cov)
        self.sum_cov = np.array(self.sum_cov)
        # check if L was specified
        if self.args.L == 0:
            self.L = int(min(self.sum_cov) / 2)
        else:
            self.L = self.args.L

    def quantile_filt(self):
        if self.args.q < 1:
            quantile = np.percentile(self.sum_cov, self.args.q * 100)
            self.rel_cov = self.rel_cov[self.sum_cov < quantile]
            self.sum_cov = self.sum_cov[self.sum_cov < quantile]

    def get1dSmudges(self, coverage_data, number_of_peaks, bandwidth = 0, depth = 1):
        if bandwidth == 0:
            coverage_kernel = gaussian_kde(coverage_data)
            bandwidth = coverage_kernel.factor
        else:
            coverage_kernel = gaussian_kde(coverage_data, bandwidth)
        coverage_range = range(max(coverage_data))
        coverage_hist = coverage_kernel.evaluate(coverage_range)
        peaks, _ = find_peaks(coverage_hist)
        if len(peaks) > number_of_peaks and depth < 10:
            new_bandwidth = bandwidth * 1.4
            # print("Nah, I need bigger bandwidth: " + str(new_bandwidth))
            peaks, _ = self.get1dSmudges(coverage_data, number_of_peaks, new_bandwidth, depth + 1)
        return peaks, coverage_hist[peaks]

# originally I used a weighted average including tetra and pendaploid smudges
# but I think diploid and tirploid will be sufficient
    def initialNEstimate(self):
        # don't estimate 1n coverage if user have speciefied it
        if self.args.n > 0:
            self.n_init = self.args.n
            return

        even_coverage_kmer_pairs = np.array([self.sum_cov[i] for i in range(len(self.rel_cov)) if self.rel_cov[i] > 0.47])
        triploid_coverage_kmer_pairs = np.array([self.sum_cov[i] for i in range(len(self.rel_cov)) if self.rel_cov[i] > 0.31 and  self.rel_cov[i] < 0.35])

        di_peaks, di_heights = self.get1dSmudges(even_coverage_kmer_pairs, 3)
        tri_peaks, tri_heights = self.get1dSmudges(triploid_coverage_kmer_pairs, 2)

        di_peaks_genome_cov = [peak / round(2 * peak / di_peaks[0]) for peak in di_peaks]
        tri_peaks_genome_cov = [peak / round(3 * peak / tri_peaks[0]) for peak in tri_peaks]

        self.n_init = np.average(np.concatenate((di_peaks_genome_cov, tri_peaks_genome_cov)),
                                 weights = np.concatenate((di_heights, tri_heights)))

    def calculateHist(self, ymin, ymax):
        self.x = np.linspace(0, 0.5, self.nbins + 1)
        self.y = np.linspace(ymin, ymax, self.nbins + 1)
        self.hist, self.y, self.x = np.histogram2d(self.sum_cov, self.rel_cov, bins=(self.y, self.x))
        # self.loghist = np.log10(self.hist)
        # self.loghist[self.loghist == -np.inf] = 0

    def agregateSmudges(self):
        #nogo_x = findInterval(.L / .smudge_container$y, .smudge_container$x, left.open = T)
        self.sorted_hist_indices = np.dstack(np.unravel_index(np.argsort(self.hist.ravel()), (self.nbins, self.nbins)))[0][::-1]
        self.smudge_assignment = []
        max_smudge = 1
        self.smudge_centers = dict()
        # x and y are coordinates in the array
        for i_to_assign, xy in enumerate(self.sorted_hist_indices):
            # I don't want to assign tiles with no kmers
            if self.hist[xy[0], xy[1]] == 0:
                self.smudge_assignment.append(0)
                continue
            # for non-zero tiles
            i_assigned = 0
            self.smudge_assignment.append(max_smudge)
            # search if there is a neibouring tile that was parsed already (== carries more kmers)
            while True:
                if (abs(xy[0] - self.sorted_hist_indices[i_assigned][0]) <= 1) and (abs(xy[1] - self.sorted_hist_indices[i_assigned][1]) <= 1):
                    self.smudge_assignment[i_to_assign] = self.smudge_assignment[i_assigned]
                    break
                i_assigned += 1
            if self.smudge_assignment[i_to_assign] == max_smudge:
                self.smudge_centers[max_smudge] = [xy]
                max_smudge += 1

    def countSmudgeSize(self, treshold = 0.02):
        all_kmers = sum(sum(self.hist))

        kmers_in_smudges = dict()
        for i, smudge_index in enumerate(self.smudge_assignment):
            kmer_pairs = self.hist[self.sorted_hist_indices[i][0], self.sorted_hist_indices[i][1]]
            try:
                kmers_in_smudges[smudge_index] += kmer_pairs
            except KeyError:
                kmers_in_smudges[smudge_index] = kmer_pairs

        for smudge_index in list(self.smudge_centers):
            rel_smudge_size = kmers_in_smudges[smudge_index] / all_kmers
            if rel_smudge_size > treshold:
                # number of kmers in the center
                self.smudge_centers[smudge_index].append(self.hist[self.smudge_centers[smudge_index][0][0], self.smudge_centers[smudge_index][0][1]])
                # absolute number kmers per smudge
                self.smudge_centers[smudge_index].append(kmers_in_smudges[smudge_index])
                # relative number kmers per smudge
                self.smudge_centers[smudge_index].append(rel_smudge_size)
            else:
                self.smudge_assignment = [i if i != smudge_index else 0 for i in self.smudge_assignment]
                del self.smudge_centers[smudge_index]

    def brightestSmudgeNEstimate(self, margin = 0.02):
        # extract smudge intensities from self.smudge_centers - data structure carring the information about assigned smudges
        smudge_intensities = self.getSmudgeStats(3, sorted = False)
        # find the smudge with the highest coverage
        bighest_smudge = list(self.smudge_centers)[smudge_intensities.index(max(smudge_intensities))]
        rec_cov_coordinate = self.smudge_centers[bighest_smudge][0][1]
        # get the relative coverage of the posisiotn of the brighest peak (likely 0.5, 0.33 or 0.25)
        rel_cov = (self.x[rec_cov_coordinate + 1] + self.x[rec_cov_coordinate]) / 2
        min_cov = rel_cov - margin
        max_cov = rel_cov + margin
        # extract kmers that are contain the birghest smudge for one more kernel smoothing fit
        subset_sum_cov = np.array([self.sum_cov[i] for i in range(len(self.rel_cov)) if self.rel_cov[i] > min_cov and self.rel_cov[i] < max_cov])
        smudge_centers, smudge_brightness = self.get1dSmudges(subset_sum_cov, 10)
        brighest_coverage = smudge_centers[np.argmax(smudge_brightness)]
        # round(brighest_coverage / self.n_init)) gives estimate of ploidy of the smudge
        # brighest_coverage / smudge_ploidy
        self.ploidy_est = int(brighest_coverage / self.n_init)
        self.brightest_smudge_n = brighest_coverage / self.ploidy_est

    def guessGenomeStructure(self):
        # -> get AB; AAB... annotations of smudges
        genome_struct_template = np.array(["A", "B"])
        for smudge_index in self.smudge_centers:
            sum_cov_index = self.smudge_centers[smudge_index][0][0]
            rel_cov_index = self.smudge_centers[smudge_index][0][1]
            sum_cov = ((self.y[sum_cov_index] + self.y[sum_cov_index + 1]) / 2)
            rel_cov = ((self.x[rel_cov_index] + self.x[rel_cov_index + 1]) / 2)
            genome_count = round(sum_cov / self.brightest_smudge_n)
            As = round((1 - rel_cov) * genome_count)
            Bs = round(rel_cov * genome_count)
            struct_list = np.repeat(genome_struct_template, [As, Bs], axis=0)
            structure = "".join(struct_list)
            self.smudge_centers[smudge_index].append(structure)

    def hasDuplicitSmudges(self):
        structures = np.array(self.getSmudgeStats(4, sorted = False))
        uniq_sctruc, cunts = np.unique(structures, return_counts=True)
        return(any(cunts > 1))

    def saveKmersFromSmudges(self):
        # build a dictionary of position > assigned smudge
        coor_smudge_dict = defaultdict(int)
        for i, assigned_smudge in enumerate(self.smudge_assignment):
             coords = self.sorted_hist_indices[i]
             coor_smudge_dict[tuple(coords)] = assigned_smudge

        # build a dictionary of smudges > list of kmers they carry
        with open(self.args.kmer_file, 'r') as kmer_file:
            smudge_kmers = defaultdict(list)
            for i, kmer in enumerate(kmer_file):
                x_coord = bisect_left(self.x, self.rel_cov[i]) - 1
                y_coord = bisect_left(self.y, self.sum_cov[i]) - 1
                assigned_smudge = coor_smudge_dict[(y_coord, x_coord)]
                smudge_kmers[assigned_smudge].append(kmer)

        # save the kmers in individual files
        for processed_smudge in smudge_kmers.keys():
            if processed_smudge == 0:
                logging.info("Skipping " + str(len(smudge_kmers[processed_smudge])) + " kmer pairs without assigned smudge (close to the error line, dispersed in space, etc.)")
                continue
            with open(self.args.o + "_kmers_in_smudge_" + str(processed_smudge) + ".txt", 'w') as kmer_smudge_file:
                for kmer in smudge_kmers[processed_smudge]:
                    kmer_smudge_file.write(kmer)

    def plot(self, ylim):
        fig = plt.figure(figsize=(8, 8))
        mpl.rcParams.update({'font.size': 14})

        ax_joint = plt.subplot2grid((8, 8), (2, 0), colspan=6, rowspan=6)
        ax_marg_x = plt.subplot2grid((8, 8), (0, 0), colspan=6, rowspan=2, sharex=ax_joint, frameon=False)
        ax_marg_y = plt.subplot2grid((8, 8), (2, 6), colspan=2, rowspan=6, sharey=ax_joint, frameon=False)
        ax_legend = plt.subplot2grid((8, 8), (0, 6), colspan=1, rowspan=2)

        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_x.yaxis.get_majorticklines(), visible=False)
        plt.setp(ax_marg_x.yaxis.get_minorticklines(), visible=False)
        # to remove even the x ticks
        # plt.setp(ax_marg_x.xaxis.get_majorticklines(), visible=False)
        # plt.setp(ax_marg_x.xaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_x.get_yticklabels(), visible=False)
        ax_marg_x.yaxis.grid(False)
        ax_marg_x.hist(self.rel_cov, self.nbins, orientation = 'vertical', color = 'darkred')
        # plot title at top left of the plot
        if self.args.title:
            ax_marg_x.text(0, 0.7, '$\it{' + self.args.title + '}$', size = 20, horizontalalignment='left',
                           verticalalignment='bottom', transform=ax_marg_x.transAxes)
        # estimated ploidy just bellow
        verbose_ploidy = self.getVerbosePloidy()
        ax_marg_x.text(0.05, 0.55, 'estimated ' + verbose_ploidy, horizontalalignment='left',
                       verticalalignment='bottom', transform=ax_marg_x.transAxes)
        ax_marg_x.text(0.05, 0.4, '1n coverage: ' + str(round(self.brightest_smudge_n, 1)), horizontalalignment='left',
                       verticalalignment='bottom', transform=ax_marg_x.transAxes)

        plt.setp(ax_marg_y.get_yticklabels(), visible=False)
        # to remove even the y ticks
        # plt.setp(ax_marg_y.yaxis.get_majorticklines(), visible=False)
        # plt.setp(ax_marg_y.yaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_y.xaxis.get_majorticklines(), visible=False)
        plt.setp(ax_marg_y.xaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_y.get_xticklabels(), visible=False)
        ax_marg_y.xaxis.grid(False)
        ax_marg_y.hist(self.sum_cov, self.nbins, orientation = 'horizontal', range = ylim, color = 'darkred')
        # ax_marg_y.text(1, 0, '1n = ' + str(round(self.brightest_smudge_n, 1)), horizontalalignment='right',
        #                verticalalignment='bottom', transform=ax_marg_y.transAxes)
        #### print summary
        het_structures = self.getSmudgeStats(4)
        het_structure_sizes = [round(x, 2) for x in self.getSmudgeStats(3)]
        self.plotTextVector(ax_marg_y, het_structures, 0.1)
        self.plotTextVector(ax_marg_y, het_structure_sizes, 0.8)

        ax_joint.set_xlabel("Normalized minor kmer coverage: B / (A + B)")
        ax_joint.set_ylabel("Total coverage of the kmer pair: A + B")
        ax_joint.set_xlim((0,0.5))
        ax_joint.xaxis.set_major_locator(plt.FixedLocator([0.5, 0.4, 0.33, 0.25, 0.2]))
        ax_joint.xaxis.set_major_formatter(plt.FixedFormatter(["1/2", "2/5", "1/3", "1/4", "1/5"]))
        ax_joint.yaxis.set_major_locator(plt.MultipleLocator(self.brightest_smudge_n))
        #### Annotate smudges
        self.annotateSmudges(ax_joint, ylim)
        ax_joint.set_ylim(ylim)
        cmap = plt.get_cmap('viridis')
        ax_joint.pcolormesh(self.x, self.y, self.hist, cmap = cmap)

        # ax_legend.cla()
        # fig.colorbar(im, ax=ax_legend)
        # fig = plt.figure()
        # ax_legend = plt.subplot2grid((4, 4), (0, 3), colspan=1, rowspan=1)
        bounds = np.linspace(0,self.hist.max(), 64)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        mpl.colorbar.ColorbarBase(ax_legend, cmap, norm, ticks = np.linspace(0,self.hist.max(), 6))
        # fig.tight_layout()
        plt.subplots_adjust(left = 0.12, bottom = 0.12, right = 0.95, top = 0.95, wspace = 0.2, hspace = 0.2)

        plt.savefig(self.args.o + "_smudgeplot_pythonic.png")

    def annotateSmudges(self, ax, ylim):
        for smudge_index in self.smudge_centers:
            sum_cov_index = self.smudge_centers[smudge_index][0][0]
            rel_cov_index = self.smudge_centers[smudge_index][0][1]
            y = (self.y[sum_cov_index] + self.y[sum_cov_index + 1]) / 2
            x = self.x[rel_cov_index] + self.x[rel_cov_index + 1]
            label = self.smudge_centers[smudge_index][4]
            if x > 0.90:
                # the labels of AB, AABB, AAABBB ... would not fit on the pic if they were centered on the smudge
                # I will plot it rightmost instead, because it should be very close to 0.5 anyway
                aln = 'right'
                x = 1
            else :
                aln = 'center'
            # print("x: " + str(x) + " y: " + str(y) + " as " + label)
            ax.text(x, (y - ylim[0]) / (ylim[1] - ylim[0]), label,
                    horizontalalignment = aln,
                    verticalalignment = 'center', transform = ax.transAxes)

    def getVerbosePloidy(self):
        ploidy_map = {
            2: "diploid",
            3: "triploid",
            4: 'tetraploid',
            5: 'pentaploid',
            6: 'hexaploid',
            7: 'heptaploid',
            8: 'octoploid'
        }
        return(ploidy_map.get(self.ploidy_est, "unknown"))

    def getSmudgeStats(self, index, sorted = True):
        if sorted:
            sizes = [self.smudge_centers[smudge_index][3] for smudge_index in self.smudge_centers]
            order_of_sizes = np.argsort(sizes)[::-1]
            return([self.smudge_centers[list(self.smudge_centers.keys())[smudge_order]][index] for smudge_order in order_of_sizes])
        else:
            return([self.smudge_centers[smudge_index][index] for smudge_index in self.smudge_centers])

    def plotTextVector(self, ax, vector, x, y = 1, aln = 'left'):
        lining = 0.045
        for i, value in enumerate(vector):
            ax.text(x, y - lining * i, value, horizontalalignment=aln,
                    verticalalignment='top', transform=ax.transAxes)

    def saveMatrix(self):
        np.savetxt(self.args.o + "_smudgematrix.tsv", self.hist, delimiter="\t", fmt='%i')

    def lowerNbins(self):
        if self.nbins > 20 :
            self.nbins = self.nbins - 5
        else :
            self.nbins = self.nbins - 2

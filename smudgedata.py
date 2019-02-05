import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import scipy.signal

# smudgedata - a class for easy control over the data related to smudgeplot
#
# data:
#  - args (output of argparse)
#    * infile (open file for reading)
#  - rel_cov (relative coverage of kmer pairs; B / (A + B), where B <= A )
#  - sum_cov (sum of coverages; A + B)
#  - x, y (borders of 2d histogram bins - x corresponds to rel_cov; y to sum_cov)
#  - hist (matrix of kmers in bins)


class smudgedata:
  def __init__(self, user_args):
    self.args = user_args
    if self.args.nbins == 0:
      self.args.iterative_bins = True
      self.args.nbins = 40
    else:
      self.args.iterative_bins = False

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
    if self.args.L == 0:
      self.args.L = int(min(self.sum_cov) / 2)

  def get1dSmudges(self, coverage_data, number_of_peaks, bandwidth = 0, depth = 1):
    if bandwidth == 0:
      coverage_kernel = gaussian_kde(coverage_data)
      bandwidth = coverage_kernel.factor
    else:
      coverage_kernel = gaussian_kde(coverage_data, bandwidth)
    coverage_range = range(len(coverage_data))
    coverage_hist = coverage_kernel.evaluate(coverage_range)
    peaks, _ = find_peaks(coverage_hist)
    if len(peaks) > number_of_peaks and depth < 10:
      new_bandwidth = bandwidth * 1.4
      peaks, _ = self.get1dSmudges(coverage_data, number_of_peaks, new_bandwidth, depth + 1)
    return peaks, coverage_hist[peaks]

# originally I used a weighted average including tetra and pendaploid smudges
# but I think diploid and tirploid will be sufficient
  def initialNEstimate(self):
    even_coverage_kmer_pairs = np.array([self.sum_cov[i] for i in range(len(self.rel_cov)) if self.rel_cov[i] > 0.47])
    triploid_coverage_kmer_pairs = np.array([self.sum_cov[i] for i in range(len(self.rel_cov)) if self.rel_cov[i] > 0.31 and  self.rel_cov[i] < 0.35])

    di_peaks, di_heights = self.get1dSmudges(even_coverage_kmer_pairs, 3)
    tri_peaks, tri_heights = self.get1dSmudges(triploid_coverage_kmer_pairs, 2)

    di_peaks_genome_cov = [peak / round(2 * peak / di_peaks[0]) for peak in di_peaks]
    tri_peaks_genome_cov = [peak / round(3 * peak / tri_peaks[0]) for peak in tri_peaks]

    self.n_init = np.average(np.concatenate((di_peaks_genome_cov, tri_peaks_genome_cov)),
                             weights = np.concatenate((di_heights, tri_heights)))

  def calculateHist(self, ymin, ymax):
    self.x = np.linspace(0, 0.5, self.args.nbins + 1)
    self.y = np.linspace(ymin, ymax, self.args.nbins + 1)
    self.hist, self.y, self.x = np.histogram2d(self.sum_cov, self.rel_cov, bins=(self.y, self.x))
    # self.loghist = np.log10(self.hist)
    # self.loghist[self.loghist == -np.inf] = 0

  def agregateSmudges(self):
    #nogo_x = findInterval(.L / .smudge_container$y, .smudge_container$x, left.open = T)
    self.sorted_hist_indices = np.dstack(np.unravel_index(np.argsort(self.hist.ravel()), (self.args.nbins, self.args.nbins)))[0][::-1]
    self.smudge_assignment = []
    max_smudge = 0
    self.smudge_centers = dict()
    # x and y are coordinates in the array
    for i_to_assign, xy in enumerate(self.sorted_hist_indices):
      i_assigned = 0
      self.smudge_assignment.append(max_smudge)
      while True:
        if (abs(xy[0] - self.sorted_hist_indices[i_assigned][0]) <= 1) or (abs(xy[1] - self.sorted_hist_indices[i_assigned][1]) <= 1):
          self.smudge_assignment[i_to_assign] = self.smudge_assignment[i_assigned]
          break
        i_assigned += 1
      if self.smudge_assignment[i_to_assign] == max_smudge:
        self.smudge_centers[max_smudge] = xy
        max_smudge += 1

  def plot(self):
    # , plot_log = False
    # if plot_log:
    #   plt.pcolormesh(self.x, self.y, self.loghist)
    # else:
    plt.pcolormesh(self.x, self.y, self.hist)
    # plt.show()
    plt.savefig(self.args.o + "_smudgeplot_pythonic.png")

  def saveMatrix(self):
    np.savetxt(self.args.o + "_smudgematrix.tsv", self.hist, delimiter="\t", fmt='%i')

  def hasDuplicitSmudges(self):
    return(True)

  def lowerNbins(self):
    if self.args.nbins > 20 :
      self.args.nbins = self.args.nbins - 5
    else :
      self.args.nbins = self.args.nbins - 2
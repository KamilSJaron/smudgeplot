import numpy as np
import matplotlib.pyplot as plt

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

  def initialNEstimate(self):
    self.n_init = 206

  def calculateHist(self, nbins, ymin, ymax):
    self.x = np.linspace(0, 0.5, nbins+1)
    self.y = np.linspace(ymin, ymax, nbins+1)
    self.hist, self.y, self.x = np.histogram2d(self.sum_cov, self.rel_cov, bins=(self.y, self.x))
    # self.loghist = np.log10(self.hist)
    # self.loghist[self.loghist == -np.inf] = 0

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
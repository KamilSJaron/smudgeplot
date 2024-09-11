
# cov2freq = defualtdict(covA, covB) -> freq
# cov2peak = dict(covA, covB) -> peak
# dict(peak) -> summit (if relevant)
# import numpy as np

import argparse
from pandas import read_csv # type: ignore
from collections import defaultdict
from statistics import mean
import matplotlib.pyplot as plt

####

parser = argparse.ArgumentParser()
parser.add_argument('infile', nargs='?', help='name of the input tsv file with covarages and frequencies.')
parser.add_argument('-nf', '-noise_filter', help='Do not agregate into smudge k-mer pairs with frequency lower than this parameter', type=int, default=50)
parser.add_argument('-d', '-distance', help='Manthattan distance of k-mer pairs that are considered neioboring for the local agregation purposes.', type=int, default=5)
parser.add_argument('--mask_errors', help='instead of reporting assignments to individual smudges, just remove all monotonically decreasing points from the error line', action="store_true", default = False)
args = parser.parse_args()

### what should be aguments at some point
# smu_file = 'data/ddSalArbu1/ddSalArbu1.k31_ploidy_converted.smu.txt'
# distance = 5
# noise_filter = 100

smu_file = args.infile
distance = args.d
noise_filter = args.nf

### load data
# cov_tab = np.loadtxt(smu_file, dtype=int)
cov_tab = read_csv(smu_file, names = ['covB', 'covA', 'freq'], sep='\t')
cov_tab = cov_tab.sort_values('freq', ascending = False)
L = min(cov_tab['covB']) # important only when --mask_errors is on

# generate a dictionary that gives us for each combination of coverages a frequency
cov2freq = defaultdict(int)
cov2peak = defaultdict(int)
# for idx, covB, covA, freq in cov_tab.itertuples():
#     cov2freq[(covA, covB)] = freq
# I can create this one when I iterate though the data though

# plt.hist(means)
# plt.hist([x for x in means if x < 100 and x > -100])
# plt.show() 

### clustering
next_peak = 1
for idx, covB, covA, freq in cov_tab.itertuples():
    cov2freq[(covA, covB)] = freq # with this I can get rid of lines 23 24 that pre-makes this table
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
        if args.mask_errors:
            if covB < L + args.d:
                cov2peak[(covA, covB)] = 1 # error line
            else:
                cov2peak[(covA, covB)] = 0 # central smudges
        else:
            cov2peak[(covA, covB)] = next_peak # if I want to keep info about all locally agregated smudges
            next_peak += 1

cov_tab = cov_tab.sort_values(['covA', 'covB'], ascending = True)
for idx, covB, covA, freq in cov_tab.itertuples():
    print(covB, covA, freq, cov2peak[(covA, covB)])
    # if idx > 20:
    #     break

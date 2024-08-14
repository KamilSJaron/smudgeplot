
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

# generate a dictionary that gives us for each combination of coverages a frequency
cov2freq = defaultdict(int)
cov2peak = defaultdict(int)
# for idx, covB, covA, freq in cov_tab.itertuples():
#     cov2freq[(covA, covB)] = freq
# I can create this one when I iterate though the data though

# SOME EXPLORATION THAT HELPED ME TO DIAGNOSE WHAT COULD BE DONE TO MAKE CLUSTERING MORE ROBUST
# mins = []
# maxs = []
# means = []
# freqs = []
# ### just calculating some stats
# for idx, covB, covA, freq in cov_tab.itertuples():
#     # print("Processing coverage pair: ", covB, covA)
#     # if idx > 10:
#     #     break
#     diffs = []
#     for xA in range(covA - distance,covA + distance + 1):
#         # for explored A coverage in neiborhood, we explore all possible B coordinates
#         distanceB = distance - abs(covA - xA)
#         for xB in range(covB - distanceB,covB + distanceB + 1):
#             xB, xA = sorted([xA, xB]) # this is to make sure xB is smaller than xA
#             if (xB == covB) and (xA == covA):
#                 continue # this does not matter in clustering, but it does for the stats
#             diff = cov2freq[(covA, covB)] - cov2freq[(xA, xB)]
#             # print('\t', xA, xB, ' diff: ', diff)
#             diffs.append(diff)
#     mins.append(min(diffs))
#     maxs.append(max(diffs))
#     means.append(mean(diffs))
#     freqs.append(freq)

# max(means)
# [(i, freqs[i], value) for i,value in enumerate(mins) if value > 0][0:40]
# Notes:
# the slope of points around could be potentially used to assess the quality of peak initiating, 
# but in the end of the day, it's really the proximity to the other peak that is the most important
# increasing the explored distance is actually ok speedwise, and makes the most sense in reducing noise.
# outstanding question is if to add blockers to stop peak initiation for small / edge peaks
# or if to sort this in post processing

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
        cov2peak[(covA, covB)] = next_peak
        next_peak += 1

cov_tab = cov_tab.sort_values(['covA', 'covB'], ascending = True)
for idx, covB, covA, freq in cov_tab.itertuples():
    print(covB, covA, freq, cov2peak[(covA, covB)])
    # if idx > 20:
    #     break

### distances
# covA = 5
# covB = 10
# distance = 2
# 3   10  2
# 4   9   2
# 4   10  1
# 4   11  2
# 5   8   2
# 5   9   1
# 5   10  0
# 5   11  1
# 5   12  2
# 6   9   2
# 6   10  1
# 6   11  2
# 7   10  2

# covA = 5
# covB = 10
# distance = 3

# covA_neigh  covB_neigh distance
# 2   10  3

# 3   9  3
# 3   10  2
# 3   11  3

# 4   8   3
# 4   9   3
# ...
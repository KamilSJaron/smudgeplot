from math import log
from math import ceil
import numpy as np
from scipy.signal import argrelextrema

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
        print(L, end = '')
    else:
        # take 99.8 quantile of kmers that are more than one in the read set
        hist_rel_cumsum = np.cumsum(hist[1:]) / np.sum(hist[1:])
        U = round_up_nice(np.argmax(hist_rel_cumsum > 0.998))
        print(U, end = '')
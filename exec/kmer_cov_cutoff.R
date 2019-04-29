#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# kmer.hist
# L or R
# load kmer histogram specified by user
kmer_hist <- read.table(args[1])$V2
# print L, U or both cutoffs
which_cutoff <- ifelse(length(args) > 1, args[2], 'both')

# rounding funciton taken from https://stackoverflow.com/a/6463946/2962344
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# U is representing 99.8% of the kmers in the dataset
U <- roundUpNice(min(which(cumsum(as.numeric(kmer_hist)) / sum(as.numeric(kmer_hist)) > 0.998)))

# find local minima
selected_points <- which(diff(sign(diff(kmer_hist))) == 2)
# Here I take the first local minimum as the treshold, but 10 is hard floor
L <- max(10, round(min(selected_points) * 1.25))

if ( which_cutoff == 'L'){
    cat(L)
} else if ( which_cutoff == 'U' ){
    cat(U)
} else {
    cat(paste(L,U, '\n'))
}
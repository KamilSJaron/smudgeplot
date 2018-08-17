#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# kmer.hist
# L or R

kmer_hist <- read.table(args[1])
which_cutoff <- args[2]

if ( which_cutoff = 'L'){
    cat(20)
}

if ( which_cutoff = 'R' ){
    cat(800)
}
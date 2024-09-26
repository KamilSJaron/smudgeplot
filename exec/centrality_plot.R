#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
input_name <- args[1]
# input_name = '~/test/daAchMill1_centralities.txt'

output_name <- gsub('.txt', '.pdf', input_name)
tested_covs <- read.table(input_name, col.names = c('cov', 'centrality'))

pdf(output_name)
    plot(tested_covs[, 'cov'], tested_covs[, 'centrality'], xlab = 'Coverage', ylab = 'Centrality [(theoretical_center - actual_center) / coverage ]', pch = 20)
dev.off()
# Sploidyplot

The goal is to have a smudge inference based on an explicit model. I hoped for a model that would make a lot of sense - based on negative binomials.

Gene - made his own EM algorithm, I think. I could not decipher it
Richard and Tianyi - made me an EM algorithm that work on simply coverage A and B; also using normal distributions


## alternative plotting

A minimalist attempt

```R
minidata <- daArtCamp1[daArtCamp1[, 'total_pair_cov'] < 20, ]
coverages_to_plot <- unique(minidata[, 'total_pair_cov'])
number_of_coverages_to_plot <- length(coverages_to_plot)
mini_ylim <- c(5, 20)
L <- 4
cols <- c(rgb(1,0,0, 0.5), rgb(1,1,0, 0.5), rgb(0,1,1, 0.5), rgb(1,0,1, 0.5), rgb(0,1,0, 0.5), rgb(0,0,1, 0.5))
# plot(1:6, pch = 20, cex = 5, col = cols)

plot_dot_smudgeplot(minidata, rep('black', 32), xlim, mini_ylim, cex = 3)
points((L - 1) / coverages_to_plot, coverages_to_plot, cex = 3, pch = 20, col = 'blue')

for( cov in 8:19){
    rect(0, cov - 0.5, 0.5, cov + 0.5, col = NA, border = 'black')
    width <- 1 / (2 * cov)
    min_ratio <- L / cov
    rect(min_ratio - width, cov - 0.5, min(0.5, min_ratio + width), cov + 0.5, col = sample(cols, 1))    
}

text(rep(0.05, number_of_coverages_to_plot), 8:19, 8:19)
```

This is more serious attempt that does not really work.

Alternative local aggregation

```bash
for ToLID in daAchMill1 daAchPtar1 daAdoMosc1 daAjuCham1 daAjuRept1 daArcMinu1 daArtCamp1 daArtMari1 daArtVulg1 daAtrBela1; do
    python3 playground/alternative_fitting/pair_clustering.py data/dicots/smu_files/$ToLID.k31_ploidy.smu.txt --mask_errors > data/dicots/peak_agregation/$ToLID.cov_tab_peaks
    Rscript playground/alternative_fitting/alternative_plotting_testing.R -i data/dicots/peak_agregation/$ToLID.cov_tab_peaks -o data/dicots/peak_agregation/$ToLID
done
```

This worked well. The agregation produced beautiful blocks, and vastly of the same shape as noticed by Richard. He suggested we should fix their shape and fit only a single parameter - coverage.

```R
smudge_tab <- read.table('data/dicots/peak_agregation/daArtMari1.cov_tab_peaks', col.names = c('covB', 'covA', 'freq', 'smudge'))

all_smudges <- unique(smudge_tab[, 'smudge'])
all_smudge_sizes <- sapply(all_smudges, function(x){ sum(smudge_tab[smudge_tab[, 'smudge'] == x, 'freq']) })

# plot(sort(all_smudge_sizes, decreasing = T) / sum(all_smudge_sizes), ylim = c(0, 1))
# sort(all_smudge_sizes, decreasing = T) / sum(all_smudge_sizes) > 0.02
# 2% of data soiunds reasonable

smudges <- all_smudges[all_smudge_sizes / sum(all_smudge_sizes) > 0.02 & all_smudges != 0]
smudge_sizes <- all_smudge_sizes[all_smudge_sizes / sum(all_smudge_sizes) > 0.02 & all_smudges != 0]

smudge_tab[, 'total_pair_cov'] <- smudge_tab[, 1] + smudge_tab[, 2]
smudge_tab[, 'minor_variant_rel_cov'] <- smudge_tab[, 1] / smudge_tab[, 'total_pair_cov']


smudge_tab_reduced <- smudge_tab[smudge_tab[, 'smudge'] %in% smudges, ]
source('playground/alternative_fitting/alternative_plotting_functions.R')

per_smudge_cov_list <- lapply(smudges, function(x){ smudge_tab_reduced[smudge_tab_reduced$smudge == x, ] })
names(per_smudge_cov_list) <- smudges

cov_sum_summary <- sapply(per_smudge_cov_list, function(x){ summary(x[, 'total_pair_cov']) } )
rel_cov_summary <- sapply(per_smudge_cov_list, function(x){ summary(x[, 'minor_variant_rel_cov']) } )

colnames(cov_sum_summary) <- colnames(rel_cov_summary) <- smudges

data.frame(smudges = smudges, total_pair_cov = round(cov_sum_summary[4, ], 1), minor_variant_rel_cov = round(rel_cov_summary[4, ], 3))

head(per_smudge_cov_list[['2']])

one_smudge <- per_smudge_cov_list[['2']]
one_smudge[one_smudge[, 'minor_variant_rel_cov'] == 0.5, ]

table(one_smudge[, 'covB'])
(one_smudge[, 'minor_variant_rel_cov'])



    # cov_range <- seq((2 * .L) - 2, max_cov_pair, length = 500)
    # lines((.L - 1)/cov_range, cov_range, lwd = 2.5, lty = 2, 

plot_peakmap(smudge_tab_reduced, xlim = c(0, 0.5), ylim = c(0, 300))

plot_seq_error_line(smudge_tab, 4)
plot_seq_error_line(smudge_tab, 13)
plot_seq_error_line(smudge_tab, 48)
plot_seq_error_line(smudge_tab, 80)

one_smudge <- per_smudge_cov_list[['4']]
min(one_smudge[ ,'total_pair_cov'])

one_smudge[one_smudge[ ,'total_pair_cov'] == 61, ]
one_smudge <- one_smudge[order(one_smudge[, 'minor_variant_rel_cov']), ]

right_part_of_the_smudge <- one_smudge[one_smudge[ ,'minor_variant_rel_cov'] > 0.2131147, ]

all_minor_var_rel_covs <- sort(unique(round(right_part_of_the_smudge[, 'minor_variant_rel_cov'], 2)))
corresponding_min_cov_sums <- sapply(all_minor_var_rel_covs, function(x){ min(right_part_of_the_smudge[round(right_part_of_the_smudge[, 'minor_variant_rel_cov'], 2) == x, 'total_pair_cov']) } )

lines(all_minor_var_rel_covs, corresponding_min_cov_sums, lwd = 3, lty = 3, col = 'red')

subtract_line <- function(rel_cov, cov_tab){
    approx_rel_cov = round(rel_cov, 2)
    band_covs = round(cov_tab[, 'minor_variant_rel_cov'], 2) == approx_rel_cov
    cov_tab[band_covs, ][which.min(cov_tab[band_covs, 'total_pair_cov']), ]
}

edge_points <- t(sapply(all_minor_var_rel_covs, subtract_line, right_part_of_the_smudge))
total_pair_cov <- sapply(1:29, function(x){edge_points[[x,5]]})
minor_variant_rel_cov <- sapply(1:29, function(x){edge_points[[x,6]]})
lm(total_pair_cov ~ minor_variant_rel_cov + I(minor_variant_rel_cov^2))

plot_isoA_line <- function (.covA, .cov_tab, .col = "black") {
    min_covB <- min(.cov_tab[, 'covB']) # should be L really
    max_covB <- .covA
    B_covs <- seq(min_covB, max_covB, length = 500)
    lines(B_covs/ (B_covs + .covA), B_covs + .covA, lwd = 2.5, lty = 2, 
        col = .col)
}

plot_isoA_line(48, smudge_tab, 'blue')
plot_isoA_line(79, smudge_tab, 'blue')
plot_isoA_line(110, smudge_tab, 'blue')
plot_isoA_line(141, smudge_tab, 'blue')
plot_isoA_line(172, smudge_tab, 'blue')

```

HA, looks great! Let's plot it on the background...

```R
library(smudgeplot)
source('playground/alternative_fitting/alternative_plotting_functions.R')
colour_ramp <- viridis(32)


smudge_tab <- read.table('data/dicots/peak_agregation/daAchMill1.cov_tab_errors', col.names = c('covB', 'covA', 'freq', 'is_error'))
# smudge_tab <- read.table('data/dicots/peak_agregation/daArtMari1.cov_tab_peaks', col.names = c('covB', 'covA', 'freq', 'smudge'))
smudge_tab[, 'total_pair_cov'] <- smudge_tab[, 1] + smudge_tab[, 2]
smudge_tab[, 'minor_variant_rel_cov'] <- smudge_tab[, 1] / smudge_tab[, 'total_pair_cov']
cov = 31.1 # this is from GenomeScope this time

plot_alt(smudge_tab[smudge_tab[, 'is_error'] != 0, ], c(0, 100), colour_ramp, T) 
plot_alt(smudge_tab[smudge_tab[, 'is_error'] != 1, ], c(0, 100), colour_ramp, T)
plot_alt(smudge_tab, c(0, 100), colour_ramp, T)
plot_iso_grid(31.1, 4, 100)
plot_smudge_labels(18.1, 100)
# .peak_points, .peak_sizes, .min_kmerpair_cov, .max_kmerpair_cov, col = "red"
dev.off()

plot_iso_grid()

plot_smudge_labels(cov, 240)
text(0.49, cov / 2, "2err", offset = 0, cex = 1.3, xpd = T, pos = 2)
```

Say we will test ploidy up to 16 (capturing up to octoploid paralogs). That makes

```R
smudge_tab_with_err <- read.table('data/dicots/peak_agregation/daAchMill1.cov_tab_errors', col.names = c('covB', 'covA', 'freq', 'is_error'))

smudge_filtering_threshold <- 0.01 # at least 1% of genomic kmers
colour_ramp <- viridis(32)

# # error band, done on non filtered data
# smudge_tab[, 'edgepoint'] <- F
# smudge_tab[smudge_tab[, 'covB'] < L + 3, 'edgepoint'] <- T
# plot_alt(smudge_tab[smudge_tab[, 'edgepoint'], ], c(0, 500), colour_ramp, T)

cov <- 19.55 # this will be unknown
L <- min(smudge_tab_with_err[, 'covB'])
smudge_tab <- smudge_tab_with_err[smudge_tab_with_err[, 'is_error'] == 0, ]
genomic_kmer_pairs <- sum(smudge_tab[ ,'freq'])

plot_alt(smudge_tab, c(0, 300), colour_ramp, T)

smudge_tab[, 'total_pair_cov'] <- smudge_tab[, 1] + smudge_tab[, 2]
smudge_tab[, 'minor_variant_rel_cov'] <- smudge_tab[, 1] / smudge_tab[, 'total_pair_cov']

plot_alt(smudge_tab, c(0, 300), colour_ramp) 
plot_all_smudge_labels(cov, 300)
dev.off()

#### isolating all smudges given cov

# total_genomic_kmer_pairs <- sum(smudge_tab[, 'freq'])

# plot_alt(smudge_container[[1]], c(0, 300), colour_ramp, T) 
# looks good!

# two functions need to be sources from the smudgeplot package here
covs_to_test <- seq(10.05, 60.05, by = 0.1)
centrality_grid <- sapply(covs_to_test, run_replicate, smudge_tab, smudge_filtering_threshold)
covs_to_test[which.max(centrality_grid)]

sapply(c(21.71, 21.72, 21.73), run_replicate, smudge_tab, smudge_filtering_threshold)
# 21.72 is our winner!

tested_covs <- test_grid_of_coverages(smudge_tab, smudge_filtering_threshold, min_to_explore, max_to_explore)
plot(tested_covs[, 'cov'], tested_covs[, 'centrality'])

```

Fixing the main package

```bash
for smu_file in  data/dicots/smu_files/*.k31_ploidy.smu.txt; do 
    ToLID=$(basename $smu_file .k31_ploidy.smu.txt); 
    smudgeplot.py all $smu_file -o data/dicots/grid_fits/$ToLID
done


```


## Homopolymer compressed testing

Datasets with lots of errors. Sacharomyces will do.

```
FastK -v -c -t4 -k31 -M16 -T4 data/Scer/SRR3265401_[12].fastq.gz -Ndata/Scer/FastK_Table_hc
hetmers -e4 -k -v -T4 -odata/Scer/kmerpairs_hc data/Scer/FastK_Table_hc

smudgeplot.py hetmers -L 4 -t 4 -o data/Scer/kmerpairs_default_e --verbose data/Scer/FastK_Table

smudgeplot.py all -o data/Scer/homopolymer_e4_wo data/Scer/kmerpairs_default_e_text.smu

smudgeplot.py all -o data/Scer/homopolymer_e4_with data/Scer/kmerpairs_hc_text.smu
```

## Other


### .smu to smu.txt

For the legacy `.smu` files, we have a convertor for to flat files.

```bash
gcc src_ploidyplot/smu2text_smu.c -o exec/smu2text_smu
exec/smu2text_smu data/ddSalArbu1/ddSalArbu1.k31_ploidy.smu | less
```
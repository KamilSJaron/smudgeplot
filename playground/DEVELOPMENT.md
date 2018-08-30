# STANDARDS

## R

- spaces around operators
- verbose, snakecase variable and function names
- arguments of functions have variable names starting by `.`

# FUNCTIONS

notes about R functions

## APPROACH ONE
estimate 1n coverage by smothing stripes of points at minor_variant_rel_cov 0.5, 0.33, 0.25, 0.2
and assuming that the firt 0.5 is the diploid peak (which is the problem I guess)
get local maxima


## APPROACH TWO
propagation of clusters from summits of the 2d histogram

LATEST peak_sizes function
LATEST peak agregation function


## VISUALIZATION
 - smudgeplot: weapper

 - plot_smudgeplot: core function for 2d histogram
 - annotate_peaks:
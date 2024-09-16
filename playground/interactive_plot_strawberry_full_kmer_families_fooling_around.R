library("methods")
library("argparse")
library("smudgeplot")
# library("hexbin")

# preprocessing
# to get simply number of memebers / family (exploration)
# cat data/strawberry_iinumae/kmer_counts_L109_U.tsv  | cut -f 1 > data/strawberry_iinumae/kmer_counts_L109_U_family_members.tsv
# awk '{row_sum = 0; row_max = 0; row_min = 10000; for (i=2; i <= NF; i++){ row_sum += $i; if ($i > row_max){row_max = $i} if ($i < row_min){row_min = $i} } print row_sum "\t" row_min "\t" row_max }' data/strawberry_iinumae/kmer_counts_L109_U.tsv > data/strawberry_iinumae/kmer_counts_L109_U_sums_min_max.tsv
# (exploration)
#
#

args <- ArgumentParser()$parse_args()
args$homozygous <- F
args$input <- 'data/Fiin/kmerpairs_k51_text.smu'
args$output = './data/Fiin/testrun'
args$title = 'F. iinumae'
args$nbins <- 40
args$L <- NULL
args$n_cov <- NULL
args$k <- 21


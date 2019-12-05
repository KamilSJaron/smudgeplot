# Smudgeplot

This tool extracts heterozygous kmer pairs from kmer dump files and performs gymnastics with them. We are able to disentangle genome structure by comparing the sum of kmer pair coverages (CovA + CovB) to their relative coverage (CovA / (CovA + CovB)). Such an approach also allows us to analyze obscure genomes with duplications, various ploidy levels, etc.

Smudgeplots are computed from raw or even better from trimmed reads and show the haplotype structure using heterozygous kmer pairs. For example:

![smudgeexample](https://user-images.githubusercontent.com/8181573/45959760-f1032d00-c01a-11e8-8576-ff0512c33da9.png)

Every haplotype structure has a unique smudge on the graph and the heat of the smudge indicates how frequently the haplotype structure is represented in the genome compared to the other structures. The image above is an ideal case, where the sequencing coverage is sufficient to beautifully separate all the smudges, providing very strong and clear evidence of triploidy.

This tool is planned to be a part of [GenomeScope](https://github.com/schatzlab/genomescope) in the near future.

## Installation

You need a program for counting kmers installed, such as [KMC](https://github.com/refresh-bio/KMC) and you should definitely run as well [GenomeScope](https://github.com/tbenavi1/genomescope2.0) (a classical kmer spectra analysis). It's not rare that both GenomeScope and smudgeplot are needed to make a sense out of the sequencing data.

To run smudgpelot, you will need `python3` with couple of pretty standard packages and `R`.

1. Download this repository

```
git clone https://github.com/KamilSJaron/smudgeplot
cd smudgeplot
```

2. Install the R package that is needed for plotting. For the installation you will need R package `devtools`, the rest of the dependencies is installed automatically.

```
Rscript install.R
```

3. install two scripts

```
install -C exec/smudgeplot.py /usr/local/bin
install -C exec/smudgeplot_plot.R /usr/local/bin
```

Congratulations, you should have `smudgeplot` operational.
Just to sure, check that the smudgeplot script works

```
smudgeplot.py --version
```

something like `Running smudgeplot v0.2.0` is expected to be printed.

If the installation procedure does not work, if you encounter any other problem, or if you would like to get help with the interpretation of your smudgeplot, please open an [issue](https://github.com/KamilSJaron/smudgeplot/issues/new).

## Usage

The input is a set of whole genome sequencing reads, the more coverage the better. The method is designed to process big datasets, don't hesitate to pull all single-end/pair-end libraries together.

The workflow is automatic, but it's not fool-proof. It requires some decisions. The usage is shown here is using [KMC](https://github.com/refresh-bio/KMC) and [GenomeScope](https://github.com/schatzlab/genomescope). We provide also a real data tutorial on our wiki: [strawberry genome analysis](https://github.com/KamilSJaron/smudgeplot/wiki/strawberry-tutorial). If you are interested in running Jellyfish instead of KMC, look at the [manual of smudgeplot with Jellyfish](https://github.com/KamilSJaron/smudgeplot/wiki/manual-of-smudgeplot-with-jellyfish).

Give KMC all the files with trimmed reads to calculate kmer frequencies and then generate a histogram of kmers:

```
mkdir tmp
ls *.fastq.gz > FILES
# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES kmer_counts tmp
kmc_tools transform kmer_counts histogram kmer_k21.hist -cx10000
```

where `-k` is the kmer length, `-m` is the approximate amount of RAM to use in GB (1 to 1024), `-ci<value>` excludes kmers occurring less than \<value\> times, `-cs` is the maximum value of a counter, `FILES` is a file name with a list of input files, `kmer_counts` is the output file name prefix, `tmp` is a temporary directory, and `-cx<value>` is the maximum value of counter to be stored in the histogram file.

The next step is to extract genomic kmers using reasonable coverage thresholds. You can either inspect the kmer spectra and choose the L (lower) and U (upper) coverage thresholds via visual inspection, or you can estimate them using command `smudgeplot.py cutoff <kmer.hist> <L/U>`. Then, extract kmers in the coverage range from `L` to `U` using `kmc_tools`. Then run `smudge_families` on the reduced file to compute the set of kmer families.

`smudge_families` is available at [tbenavi1/KMC](https://github.com/tbenavi1/KMC). Specifically, after compiling this forked version of KMC, `smudge_families` will be in the bin directory.

```
L=$(smudgeplot.py cutoff kmer_k21.hist L)
U=$(smudgeplot.py cutoff kmer_k21.hist U)
echo $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000
kmc_tools transform kmer_counts -ci$L -cx$U reduce kmer_counts_L$L\_U$U
smudge_families kmer_counts_L$L\_U$U kmer_counts_L$L\_U$U_coverages.tsv
```

now you can finally generate the smudgeplot using the coverages of the identified kmer families (`*_coverages.tsv` file). You can either supply the haploid kmer coverage (reported by GenomeScope) or let it be estimated directly from the data. Something like this

```
smudgeplot.py plot kmer_pairs_coverages.tsv
```

will generate a basic smugeplot. To see all the parameters of `smudgeplot.py plot` you can run `smudgeplot.py plot --help`.

Smudgeplot generates two plots, one with coloration on a log scale and on a linear scale. The legend indicates approximate kmer pairs per tile densities. Note that a single polymorphism generates multiple heterozygous kmers. As such, the reported numbers do not directly correspond to the number of variants. Instead, the actual number is approximately 1/k times the reported numbers, where k is the kmer size (in summary already recalculated). It's important to note that this process does not exhaustively attempt to find all of the heterozygous kmers from the genome. Instead, only a sufficient sample is obtained in order to identify relative genome structure. You can also report the minimal number of loci that are heterozygous if the inference is correct.

### GenomeScope for estimation of L/U

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use the [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmer_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

This script estimates the size, heterozygosity, and repetitive fraction of the genome. By inspecting the fitted model you can determine the location of the smallest peak after the error tail. Then, you can decide the low end cutoff below which all kmers will be discarded as errors (cca 1/2 of the haploid kmer coverage), and the high end cutoff above which all kmers will be discarded (cca 8x of the haploid kmer coverage).

## Frequently Asked Questions

Are collected on [our wiki](https://github.com/KamilSJaron/smudgeplot/wiki/FAQ). Smudgeplot is not much demanding on computational resources, but make sure you check [memory requirements](https://github.com/KamilSJaron/smudgeplot/wiki/smudgeplot-hetkmers#memory-requirements) before you extract kmer pairs (`hetkmers` task). If you won't find an answer for your question in FAQ, open an [issue](https://github.com/KamilSJaron/smudgeplot/issues/new/choose) or drop us an email.

Check [projects](https://github.com/KamilSJaron/smudgeplot/projects) to see how the development goes.

## Contributions

This is definitely an open project, contributions are welcome. You can check some of the ideas for future in [projects](https://github.com/KamilSJaron/smudgeplot/projects) and the development in [dev](https://github.com/KamilSJaron/smudgeplot/tree/dev) branch.

The file [playground/DEVELOPMENT.md](playground/DEVELOPMENT.md) contains some development notes. The directory [playground](playground) contains some snippets, attempts, and other items of interest.

## Reference

GenomeScope 2.0 and Smudgeplots: Reference-free profiling of polyploid genomes Timothy Rhyker Ranallo-Benavidez, Kamil S. Jaron, Michael C. Schatz bioRxiv 747568; doi: https://doi.org/10.1101/747568

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison has largely inspired the visual output of smudgeplots. The colourblind friendly colour theme was suggested by @ggrimes. Grateful for helpful comments of beta testers and pre-release chatters!

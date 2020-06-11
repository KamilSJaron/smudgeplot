# Smudgeplot

This tool extracts heterozygous kmer pairs from kmer count databases and performs gymnastics with them. We are able to disentangle genome structure by comparing the sum of kmer pair coverages (CovA + CovB) to their relative coverage (CovB / (CovA + CovB)). Such an approach also allows us to analyze obscure genomes with duplications, various ploidy levels, etc.

Smudgeplots are computed from raw or even better from trimmed reads and show the haplotype structure using heterozygous kmer pairs. For example:

![smudgeexample](https://user-images.githubusercontent.com/8181573/45959760-f1032d00-c01a-11e8-8576-ff0512c33da9.png)

Every haplotype structure has a unique smudge on the graph and the heat of the smudge indicates how frequently the haplotype structure is represented in the genome compared to the other structures. The image above is an ideal case, where the sequencing coverage is sufficient to beautifully separate all the smudges, providing very strong and clear evidence of triploidy.

This tool is planned to be a part of [GenomeScope](https://github.com/tbenavi1/genomescope2.0) in the near future.

## Installation

You will need a program for counting kmers such as [KMC](https://github.com/refresh-bio/KMC) installed, and you should definitely also run [GenomeScope](https://github.com/tbenavi1/genomescope2.0) (a classical kmer spectra analysis). It's not rare that both GenomeScope and Smudgeplot are needed to make sense out of the sequencing data. We recommend using [tbenavi1/KMC](https://github.com/tbenavi1/KMC). This version of KMC has an additional `smudge_pairs` program which finds heterozygous kmer pairs more quickly and using less memory than the python version of hetkmers.

To run Smudgeplot, you will need `python3` with a couple of standard packages and `R`.

1. Download this repository

```
git clone https://github.com/KamilSJaron/smudgeplot
cd smudgeplot
```

2. Install the R package that is needed for plotting. For the installation you will need the R package `devtools`, and the rest of the dependencies are installed automatically.

```
Rscript install.R
```

3. install two scripts

```
install -C exec/smudgeplot.py /usr/local/bin
install -C exec/smudgeplot_plot.R /usr/local/bin
```

Congratulations, `smudgeplot` should be operational.
Just to be sure, check that the smudgeplot script works:

```
smudgeplot.py --version
```

Something like `Running smudgeplot v0.2.2` is expected to be printed.

If the installation procedure does not work, if you encounter any other problem, or if you would like to get help with the interpretation of your smudgeplot, please open an [issue](https://github.com/KamilSJaron/smudgeplot/issues/new).

## Usage

The input is a set of whole genome sequencing reads, the more coverage the better. The method is designed to process big datasets, don't hesitate to pull all single-end/pair-end libraries together.

The workflow is automatic, but it's not fool-proof. It requires some decisions. The usage is shown here using [tbenavi1/KMC](https://github.com/tbenavi1/KMC) and [GenomeScope](https://github.com/tbenavi1/genomescope2.0). We also provide a real data tutorial on our wiki: [strawberry genome analysis](https://github.com/KamilSJaron/smudgeplot/wiki/strawberry-tutorial). If you are interested in running Jellyfish instead of KMC, look at the [manual of smudgeplot with Jellyfish](https://github.com/KamilSJaron/smudgeplot/wiki/manual-of-smudgeplot-with-jellyfish).

Give KMC all the files with trimmed reads to calculate kmer frequencies and then generate a histogram of kmers:

```
mkdir tmp
ls *.fastq.gz > FILES
# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES kmcdb tmp
kmc_tools transform kmcdb histogram kmcdb_k21.hist -cx10000
```

where `-k` is the kmer length, `-m` is the approximate amount of RAM to use in GB (1 to 1024), `-ci<value>` excludes kmers occurring less than \<value\> times, `-cs` is the maximum value of a counter, `FILES` is a file name with a list of input files, `kmcdb` is the output file name prefix for the KMC database, `tmp` is a temporary directory, and `-cx<value>` is the maximum value of counter to be stored in the histogram file.

The next step is to extract genomic kmers using reasonable coverage thresholds. You can either inspect the kmer spectra and choose the L (lower) and U (upper) coverage thresholds via visual inspection, or you can estimate them using command `smudgeplot.py cutoff <kmcdb_k21.hist> <L/U>`.
```
L=$(smudgeplot.py cutoff kmcdb_k21.hist L)
U=$(smudgeplot.py cutoff kmcdb_k21.hist U)
echo $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000
```

Then, extract kmers in the coverage range from `L` to `U` using `kmc_tools`. Then run `smudge_pairs` on the reduced file to compute the set of kmer pairs.

`smudge_pairs` is available at [tbenavi1/KMC](https://github.com/tbenavi1/KMC). Specifically, after compiling this version of KMC, `smudge_pairs` will be in the bin directory.

```
kmc_tools transform kmcdb -ci$L -cx$U reduce kmcdb_L$L\_U$U
smudge_pairs kmcdb_L$L\_U$U kmcdb_L$L\_U$U_coverages.tsv kmcdb_L$L\_U$U_pairs.tsv > kmcdb_L$L\_U$U_familysizes.tsv
```

Alternatively, if you don't have [tbenavi1/KMC](https://github.com/tbenavi1/KMC) installed, you can extract kmers in the coverage range from `L` to `U` using `kmc_dump`. Then run `smudgeplot.py hetkmers` on the dump of kmers the compute the set of kmer pairs.
```
kmc_tools transform kmcdb -ci$L -cx$U dump -s kmcdb_L$L\_U$U.dump
smudgeplot.py hetkmers -o kmcdb_L$L\_U$U < kmcdb_L$L\_U$U.dump
```

Now you can finally generate the smudgeplot using the coverages of the identified kmer pairs (`*_coverages.tsv` file). You can either supply the haploid kmer coverage (reported by GenomeScope) or let it be estimated directly from the data. Something like this

```
smudgeplot.py plot kmcdb_L$L\_U$U_coverages.tsv
```

will generate a basic smudgeplot. To see all the parameters of `smudgeplot.py plot` you can run `smudgeplot.py plot --help`.

Smudgeplot generates two plots, one with coloration on a log scale and the other on a linear scale. The legend indicates approximate kmer pairs per tile densities. Note that a single polymorphism generates multiple heterozygous kmers. As such, the reported numbers do not directly correspond to the number of variants. Instead, the actual number is approximately 1/k times the reported numbers, where k is the kmer size (in summary already recalculated). It's important to note that this process does not exhaustively attempt to find all of the heterozygous kmers from the genome. Instead, only a sufficient sample is obtained in order to identify relative genome structure. You can also report the minimal number of loci that are heterozygous if the inference is correct.

### GenomeScope for estimation of L/U

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use the [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmcdb_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

This script estimates the size, heterozygosity, and repetitive fraction of the genome. By inspecting the fitted model you can determine the location of the smallest peak after the error tail. Then, you can decide the low end cutoff below which all kmers will be discarded as errors (cca 0.5 times the haploid kmer coverage), and the high end cutoff above which all kmers will be discarded (cca 8.5 times the haploid kmer coverage).

## Frequently Asked Questions

Are collected on [our wiki](https://github.com/KamilSJaron/smudgeplot/wiki/FAQ). Smudgeplot does not demand much on computational resources, but make sure you check [memory requirements](https://github.com/KamilSJaron/smudgeplot/wiki/smudgeplot-hetkmers#memory-requirements) before you extract kmer pairs (`hetkmers` task). If you don't find an answer for your question in FAQ, open an [issue](https://github.com/KamilSJaron/smudgeplot/issues/new/choose) or drop us an email.

Check [projects](https://github.com/KamilSJaron/smudgeplot/projects) to see how the development goes.

## Contributions

This is definitely an open project, contributions are welcome. You can check some of the ideas for the future in [projects](https://github.com/KamilSJaron/smudgeplot/projects) and in the development [dev](https://github.com/KamilSJaron/smudgeplot/tree/dev) branch.

The file [playground/DEVELOPMENT.md](playground/DEVELOPMENT.md) contains some development notes. The directory [playground](playground) contains some snippets, attempts, and other items of interest.

## Reference

Ranallo-Benavidez, T.R., Jaron, K.S. & Schatz, M.C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* **11**, 1432 (2020). https://doi.org/10.1038/s41467-020-14998-3

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison has largely inspired the visual output of smudgeplots. The colourblind friendly colour theme was suggested by @ggrimes. Grateful for helpful comments of beta testers and pre-release chatters!

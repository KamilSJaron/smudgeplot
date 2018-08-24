# Sibling kmer analysis

This tool extract heterozygous kmer pairs from jellyfish dump files and does a gymnastics with it. We are able to disentangle genome structure by comparison of the sum of kmer pair coverages and their relative coverage. Such approach allows to analyze also obscure genomes with duplications, various ploidy levels and stuff. it's work in progress but it's already great.

## Installation

You need [jellyfish](https://github.com/gmarcais/Jellyfish) install and to get most of your kmer spectra, also [GenomeScope](https://github.com/schatzlab/genomescope).

TODO: add required python libraries and R packages (you can find them in the scripts)

To install you can run

```
git clone https://github.com/tbenavi1/sibling_kmers
cd sibling_kmers
make install
```

This will install R package `smudgeplot` to place where all the other R packages are installed and two scripts `smudgeplot.R` and `hetkmers.py` to `/usr/local/bin`. If you if you wish to specify directory where to install you can define location as `make install INSTAL_PREFIX=/home/slim/` (installs to `/home/slim/bin`).

If the installation procedure does not work, if you encounter any other problem or if you would like to get a help with interpretation of your smudgeplot, please open an [issue](https://github.com/tbenavi1/sibling_kmers/issues/new).

## Usage

Right now the workflow is semi-automatic. It requires some decisions. The first step is to run [jellyfish](https://github.com/gmarcais/Jellyfish) and [GenomeScope](https://github.com/schatzlab/genomescope).

Give jellyfish all the files with trimmed reads to calculate kmer frequencies and then generate histogram of kmers:

```
jellyfish count -C -m 21 -s 1000000000 -t 8 *.fastq -o kmer_counts.jf
jellyfish histo kmer_counts.jf > kmer_k21.hist
```

Now feed the histogram to GenomeScope (either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use their [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmer_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

You just made estimates of genome size, heterozygosity and repetitive fraction of the genome. Look at the fitted model and figure out what is the smallest peak after the error tail. Now you need to decide the low end cutoff that will discard all the kmers that have that small coverage and look like errors (cca 1/2 of the haploid kmer coverage), also check the right side of the graph to find reasonable upped threshold (cca 8x of the haploid kmer coverage). Dump kmers in the coverage range you have defined using jellyfish and pipe it directly to the python script that will compute set of heterozygous kmers.

```
jellyfish dump -c -L 40 -U 600 kmer_counts.jf | hetkmers.py -o kmer_pairs
```

Now finally generate smudgeplot using coverages of the kmer pairs. You need to supply the haploid kmer coverage (reported by GenomeScope). If GenomeScope correctly identified the peak of haplod kmers, the expected positions of haplotype structures will overlap with high density smudges on the smudgeplot. If the overlap is not great you might consider to adjust both GenomeScope model and redo the plod with the better estimate of the haploid coverage.

```
smudgeplot.R [input.tsv] [output.png] [plot_title] [haplod_cov]
```

The smudgeplot uses colouration on squared scale, the legend indicate approximate kmer pairs per tile densities. Note that single polymorphism generates multiple heterozygous kmers, these numbers do not directly correspond to variants, but should be fairly correlated.

## Content

(Moreless) working code:
- [hetkmers.py](hetkmers.py): python script to compute coverages heterozygous kmer pairs (`*_coverages_2.tsv` file)
- `smudgeplot`: an R package for plotting smudgeplot (`DESCRIPTION` file, `R` and `tests` directories)
- [smudgeplot.R](smudgeplot.R): R script wraps `smudgeplot` package functions. The smudgeplot is generated from `*_coverages_2.tsv` file - a two column file with kmer sorted kmer coverages.

Under development:
- [DEVELOPMENT.md](DEVELOPMENT.md): some development notes
- [playground.R](playground.R): snippets of R code for testing that should be maintained.
- [playground.py](playground.py): snippets of python code for testing that should be maintained.
- [kmer_cov_cutoff.R](kmer_cov_cutoff.R): a script that will be used in future to create L and U cutoffs for kmer spectra dumping
- [more_away_pairs.py](more_away_pairs.py): functions that allow to extract kmers that are further away (but for a great ocmputational cost)

## Contributions

TODO

## License

TODO GPL3? MIT? CC-0?
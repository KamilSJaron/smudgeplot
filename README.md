# Sibling kmer analysis

This tool extract heterozygous kmer pairs from jellyfish dump files and does a gymnastics with it. We are able to disentangle genome stricture by comparison of the sum of kmer pair coverages and their relative coverage. Such approach allows to analyze also obscure genomes with duplications, various ploidy levels and stuff. it's work in progress but it's going to be great.

## Content

- [hetkmers.py](hetkmers.py): python script to compute heterozygous kmer pairs
- [smudgeplot.R](smudgeplot.R): R script that creates so-called smudgeplot from `coverages_2.tsv` file - a two column file with kmer sorted kmer coverages.

- [playground.py](playground.py): snippets of python code for testing that should be maintained.

## Usage

Right now the workflow is semi-automatic. It requires some decisions. The first step is to run [jellyfish](https://github.com/gmarcais/Jellyfish) and [GenomeScope](https://github.com/schatzlab/genomescope).

Give jellyfish all the files with trimmed reads to calculate kmer frequencies and then generate histogram of kmers:

```
jellyfish count -C -m 21 -s 1000000000 -t 8 *.fastq -o kmer_counts.jf
jellyfish histo kmer_counts.jf > kmer_k21.hist
```

Now feed the histogram to GenomeScope (either run the [script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use their [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmer_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

You just made estimates of genome size, heterozygosity and repetitive fraction of the genome. Look at the fitted model and figure out what is the smallest peak after the error tail. Now you need to decide the low end cutoff that will discard all the kmers that have that small coverage and look like errors (cca 1/2 of the haplod kmer coverage), also check the right side of the graph to find reasonable upped threshold (cca 8x of the haploid kmer coverage). Dump kmers in the coverage range you have defined using jellyfish and pipe it directly to the python script that will compute set of heterozygous kmers.

```
jellyfish dump -c -L 40 -U 600 kmer_counts.jf | hetkmers.py -o kmer_pairs
```

Now finally generate smudgeplot using coverages of the kmer pairs. You need to supply the haploid kmer coverage (reported by GenomeScope). If GenomeScope correctly identified the peak of haplod kmers, the expected positions of haplotype structures will overlap with high density smudges on the smudgeplot. If the overlap is not great you might consider to adjust both GenomeScope model and redo the plod with the better estimate of the haploid coverage.

```
smudgeplot.R [input.tsv] [output.png] [plot_title] [haplod_cov]
```

The smudgeplot uses colouration on squared scale, the legend indicate approximate kmer pairs per tile densities. Note that single polymorphism generates multiple heterozygous kmers, these numbers do not directly correspond to variants, but should be fairly correlated.

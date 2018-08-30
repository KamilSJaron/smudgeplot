# Smudgeplot

This tool extract heterozygous kmer pairs from jellyfish dump files and perform gymnastics with it. We are able to disentangle genome structure by comparison of the sum of kmer pair coverages and their relative coverage. Such approach allows to analyze also obscure genomes with duplications, various ploidy levels and stuff. it's a work in progress but it's already great.

This tool is planned to be a part of [GenomeScope](https://github.com/schatzlab/genomescope) in future.

## Installation

You need [jellyfish](https://github.com/gmarcais/Jellyfish) installed and if you want to get most of your kmer spectra, also consider running [GenomeScope](https://github.com/schatzlab/genomescope).

You will need python libraries `networkx` and `pandas`. Both of them are available in `pip`. Required R libraries are installed together with the R library `smudgeplot`.

Get this repository

```
git clone https://github.com/tbenavi1/sibling_kmers
cd sibling_kmers
```

Install the R package `smudgeplot`

```
Rscript install.R    # installs R library
```

Finally copy this three scripts somewhere where you system will see it (places in `$PATH`)

```
# copy scripts somewhere where your shell will see them
install -C hetkmers.py /usr/local/bin
install -C smudgeplot.R /usr/local/bin
install -C kmer_cov_cutoff.R /usr/local/bin
```

If you don't have rights to write at these directories you can alternatively add this directory to your PATH

```
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

If the installation procedure does not work, if you encounter any other problem or if you would like to get a help with interpretation of your smudgeplot, please open an [issue](https://github.com/tbenavi1/sibling_kmers/issues/new).

## Usage

Right now the workflow is automatic, but it's not fool-proof. It requires some decisions. The first step is to run [jellyfish](https://github.com/gmarcais/Jellyfish) and [GenomeScope](https://github.com/schatzlab/genomescope).

Give jellyfish all the files with trimmed reads to calculate kmer frequencies and then generate histogram of kmers:

```
jellyfish count -C -m 21 -s 1000000000 -t 8 *.fastq -o kmer_counts.jf
jellyfish histo kmer_counts.jf > kmer_k21.hist
```

The next step is to extract genomic kmers using reasonable coverage threshold. You can either inspect the kmer spectra and chose the L (lower) and U (upper) coverage thresholds via visual inspection, or you can estimate them using script `kmer_cov_cutoff.R <kmer.hist> <L/U>`. Extract kmers in the coverage range from `L` to `U` using `Jellyfish` and pipe it directly to the python script `hetkmers.py`  that will compute set of heterozygous kmers.

```
L=$(kmer_cov_cutoff.R kmer_k21.hist L)
U=$(kmer_cov_cutoff.R kmer_k21.hist U)
echo $L $U # these need to be sane values like 30 800 or so
jellyfish dump -c -L $L -U $U kmer_counts.jf | hetkmers.py -o kmer_pairs
```

Now finally generate smudgeplot using coverages of the kmer pairs (`*_coverages_2.tsv` file). You can either supply the haploid kmer coverage (reported by GenomeScope) or you can let it to be estimated form directly from the data and compare it afterwards. If GenomeScope correctly identified the peak of haploid kmers, the expected positions of haplotype structures will overlap with high density smudges on the smudgeplot. If the overlap is not great you might consider to adjust both GenomeScope model and redo the plot with the better estimate of the haploid coverage. Usage of `smudgeplot.R` script is

```
usage: ./smudgeplot.R [-h] [-v] [-i INPUT] [-o OUTPUT] [-t TITLE] [-n N_COV]
                      [-L LOW_CUTOFF] [-nbins NBINS] [-k KMER_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print the version and exit
  -i INPUT, --input INPUT
                        name of the input tsv file with covarages [default
                        "coverages_2.tsv"]
  -o OUTPUT, --output OUTPUT
                        name pattern used for the output files
                        (OUTPUT_smudgeplot.png, OUTPUT_summary.txt,
                        OUTPUT_warrnings.txt) [default "smudgeplot"]
  -t TITLE, --title TITLE
                        name printed at the top of the smudgeplot [default
                        none]
  -n N_COV, --n_cov N_COV
                        the haploid coverage of the sequencing data [default
                        inference from data]
  -L LOW_CUTOFF, --low_cutoff LOW_CUTOFF
                        the lower boundary used when dumping kmers from
                        jellyfish [default min(total_pair_cov) / 2]
  -nbins NBINS          the number of nbins used for smudgeplot matrix (nbins
                        x nbins) [default 40]
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        The kmer size used to calculate kmer spectra [default
                        21]
```

The smudgeplot uses colouration on squared scale, the legend indicate approximate kmer pairs per tile densities. Note that single polymorphism generates multiple heterozygous kmers, these numbers do not directly correspond to variants, it's approximately kmer size times less. It's important to note that we are obviouslty not fishing out all the heterozygous kmers from the genome, we just get enough to identify relative genome structure. We are also able as well report minimal number of loci that are heterozygous if the inference was correct.

### GenomeScope for estimation of L/U

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use their [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmer_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

You just made estimates of genome size, heterozygosity and repetitive fraction of the genome. Look at the fitted model and figure out what is the smallest peak after the error tail. Now you need to decide the low end cutoff that will discard all the kmers that have that small coverage and look like errors (cca 1/2 of the haploid kmer coverage), also check the right side of the graph to find reasonable upped threshold (cca 8x of the haploid kmer coverage).

## Obvious near future

There is one amazing thing we can potentially do with smudges (not implemented yet). We could use smudges to quantify error in the genome assembly. The arch nemesis of variant calling is missmapping mostly because of collapsed paralogs or separately assembled alleles. We can not fix those mistakes, but we can quantify them.

Imagine you have diploid organism and supposedly haploid assembly. If we extract kmer pairs with inferred four genomic copies, we expect to find them exactly twice (twice because haplotypes are expected to be collapsed). Therefore we can quantify what is the proportion of four-copy heterozygous kmers was assembled in corrected copy number. Analogically we can extract heterozygous kmers in two copies (those that should be collapsed into one sequence in the assembly) and map them to genome. Those that will be identified two

## Contributions

In file [DEVELOPMENT.md](playground/DEVELOPMENT.md) there are some development notes. Directory [playground](playground) contains some snippets, attempts and stuff that some of us at some point wanted to store.

TODO

## License

TODO GPL3? MIT? CC-0?

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison have largely inspired visual of smudgeplots (including colour theme in fact).
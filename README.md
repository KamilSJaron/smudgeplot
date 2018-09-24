# Smudgeplot

This tool extracts heterozygous kmer pairs from kmer dump files (from jellyfish or KMC) and performs gymnastics with them. We are able to disentangle genome structure by comparing the sum of kmer pair coverages (CovA + CovB) to their relative coverage (CovA / (CovA + CovB)). Such an approach also allows us to analyze obscure genomes with duplications, various ploidy levels, etc. It's a work in progress but already provides great insights.

Smudgeplots are computed from raw/trimmed reads and show the haplotype structure using heterozygous kmer pairs. For example:

![smudgeexample](https://user-images.githubusercontent.com/8181573/44874203-ac45d980-ac9a-11e8-9943-889acbba81cd.png)

Every haplotype structure has a unique smudge on the graph and the heat of the smudge indicates how frequently the haplotype structure is represented in the genome compared to the other structures. The image above is an ideal case, where the sequencing coverage is sufficient to beautifully separate all the smudges, providing very strong and clear evidence of triploidy.

This tool is planned to be a part of [GenomeScope](https://github.com/schatzlab/genomescope) in the near future.

## Installation

You need [jellyfish](https://github.com/gmarcais/Jellyfish) or [KMC](https://github.com/refresh-bio/KMC) installed and if you want to get the most out of your kmer spectra, also consider running [GenomeScope](https://github.com/schatzlab/genomescope).

Required R libraries are installed together with the R library `smudgeplot`.

Get this repository

```
git clone https://github.com/tbenavi1/sibling_kmers
cd sibling_kmers
```

Install the R package `smudgeplot`

```
Rscript install.R    # installs R library
```

Finally copy these three scripts somewhere where your system will see them (places in `$PATH`)

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

If the installation procedure does not work, if you encounter any other problem, or if you would like to get help with the interpretation of your smudgeplot, please open an [issue](https://github.com/tbenavi1/sibling_kmers/issues/new).

## Usage

The input is a set of whole genome sequencing reads, the more coverage the better. The method is designed to process big datasets, don't hesitate to pull all single-end/pair-end libraries together.

Right now the workflow is automatic, but it's not fool-proof. It requires some decisions. The first step is to run [jellyfish](https://github.com/gmarcais/Jellyfish) or [KMC](https://github.com/refresh-bio/KMC) and optionally also [GenomeScope](https://github.com/schatzlab/genomescope).

If using jellyfish, give jellyfish all the files with trimmed reads to calculate kmer frequencies and then generate a histogram of kmers:

```
jellyfish count -C -m 21 -s 1000000000 -t 8 *.fastq -o kmer_counts.jf
jellyfish histo kmer_counts.jf > kmer_k21.hist
```

If using KMC, give KMC all the files with trimmed reads to calculate kmer frequencies and then generate a histogram of kmers:

```
mkdir tmp
ls *.fastq > FILES
kmc -k21 -m64 -ci1 -cs10000 @FILES kmer_counts tmp
kmc_tools transform kmer_counts histogram kmer_k21.hist
```

where `-k` is the kmer length, `-m` is the max amount of RAM to use in GB (1 to 1024), `-ci<value>` excludes kmers occurring less than \<value\> times, `-cs` is the maximum value of a counter, `FILES` is a file name with a list of input files, `kmer_counts` is the output file name prefix, and `tmp` is a temporary directory.

The next step is to extract genomic kmers using reasonable coverage thresholds. You can either inspect the kmer spectra and choose the L (lower) and U (upper) coverage thresholds via visual inspection, or you can estimate them using the script `kmer_cov_cutoff.R <kmer.hist> <L/U>`. Then, extract kmers in the coverage range from `L` to `U` using `Jellyfish` or `KMC` and pipe them directly to the python script `hetkmers.py` to compute the set of heterozygous kmers.

If using jellyfish

```
L=$(kmer_cov_cutoff.R kmer_k21.hist L)
U=$(kmer_cov_cutoff.R kmer_k21.hist U)
echo $L $U # these need to be sane values like 30 800 or so
jellyfish dump -c -L $L -U $U kmer_counts.jf | hetkmers.py -k 21 -t 8 -o kmer_pairs
```

If using KMC

```
L=$(kmer_cov_cutoff.R kmer_k21.hist L)
U=$(kmer_cov_cutoff.R kmer_k21.hist U)
echo $L $U # these need to be sane values like 30 800 or so
kmc_dump -ci$L -cx$U kmer_counts kmer_k21.dump
hetkmers.py -k 21 -t 8 -o kmer_pairs < kmer_k21.dump
```

After using jellyfish or KMC, generate the smudgeplot using the coverages of the kmer pairs (`*_coverages_2.tsv` file). You can either supply the haploid kmer coverage (reported by GenomeScope) or let it be estimated directly from the data and compare it afterwards. If GenomeScope correctly identifies the peak of haploid kmers, the expected positions of the haplotype structures will overlap with high density smudges on the smudgeplot. If the overlap is not great you might consider adjusting the GenomeScope model and redoing the plot with a better estimate of the haploid coverage. Usage of `smudgeplot.R` script is

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

The smudgeplot uses coloration on a log scale. The legend indicates approximate kmer pairs per tile densities. Note that a single polymorphism generates multiple heterozygous kmers. As such, the reported numbers do not directly correspond to the number of variants. Instead, the actual number is approximately 1/k times the reported numbers, where k is the kmer size. It's important to note that this process does not exhaustively attempt to find all of the heterozygous kmers from the genome. Instead, only a sufficient sample is obtained in order to identify relative genome structure. You can also report the minimal number of loci that are heterozygous if the inference is correct.

### GenomeScope for estimation of L/U

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use the [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmer_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

This script estimates the size, heterozygosity, and repetitive fraction of the genome. By inspecting the fitted model you can determine the location of the smallest peak after the error tail. Then, you can decide the low end cutoff below which all kmers will be discarded as errors (cca 1/2 of the haploid kmer coverage), and the high end cutoff above which all kmers will be discarded (cca 8x of the haploid kmer coverage).

## Obvious near future

One important potential feature of smudgeplots (not implemented yet) is to quantify error in the genome assembly. As we know, the arch nemesis of variant calling is missmapping mostly because of collapsed paralogs or separately assembled alleles. Smudgeplots can be used to quantify such mistakes.

For example, imagine you have a diploid organism and a supposedly haploid assembly. If you extract kmer pairs with inferred four genomic copies, you would expect to find them exactly twice (because the haplotypes are expected to be collapsed). Therefore you can quantify the proportion of four-copy heterozygous kmers was assembled in corrected copy number. Analogically we can extract heterozygous kmers in two copies (those that should be collapsed into one sequence in the assembly) and map them to genome. Those that will be identified two

## Contributions

The file [DEVELOPMENT.md](playground/DEVELOPMENT.md) contains some development notes. The directory [playground](playground) contains some snippets, attempts, and other items of interest.

TODO

## License

TODO GPL3? MIT? CC-0?

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison has largely inspired the visual output of smudgeplots (including the color theme).

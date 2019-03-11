# Warning

you are stepping on developers' ground. Thank you for that, feedback before release is always welcome.

Recently we undergone a switch form R to python code (most for versatility and performance purposes), but that might mean that some of the stuff is not working as it should and it's possible that I am not aware of it. If you encounter any unexpected behavior, open an issue or send us a mail.

Cheers,
Kamil

# Smudgeplot

This tool extracts heterozygous kmer pairs from kmer dump files and performs gymnastics with them. We are able to disentangle genome structure by comparing the sum of kmer pair coverages (CovA + CovB) to their relative coverage (CovA / (CovA + CovB)). Such an approach also allows us to analyze obscure genomes with duplications, various ploidy levels, etc. It's a work in progress but already provides great insights.

Smudgeplots are computed from raw/trimmed reads and show the haplotype structure using heterozygous kmer pairs. For example:

![smudgeexample](https://user-images.githubusercontent.com/8181573/45959760-f1032d00-c01a-11e8-8576-ff0512c33da9.png)

Every haplotype structure has a unique smudge on the graph and the heat of the smudge indicates how frequently the haplotype structure is represented in the genome compared to the other structures. The image above is an ideal case, where the sequencing coverage is sufficient to beautifully separate all the smudges, providing very strong and clear evidence of triploidy.

This tool is planned to be a part of [GenomeScope](https://github.com/schatzlab/genomescope) in the near future.

## Installation

You need [jellyfish](https://github.com/gmarcais/Jellyfish) or [KMC](https://github.com/refresh-bio/KMC) installed and we heavily recommend running [GenomeScope](https://github.com/schatzlab/genomescope) as well, sometimes both GenomeScope and smudgeplot are needed to make a sense out of the sequencing data.

Get this repository

```
git clone https://github.com/tbenavi1/smudgeplot
cd smudgeplot
```

Install the python package `smudgeplot`

```
python setup.py install --user    # installs the library
```

Congratulations, you should have `smudgeplot` operational.
Just to sure, check that the smudgeplot script works

```
smudgeplot -v
```

something like `INFO:root:Running smudgeplot v0.1.3 beta3-development` is expected to be printed.

If the installation procedure does not work, if you encounter any other problem, or if you would like to get help with the interpretation of your smudgeplot, please open an [issue](https://github.com/tbenavi1/smudgeplot/issues/new).

## Usage

The input is a set of whole genome sequencing reads, the more coverage the better. The method is designed to process big datasets, don't hesitate to pull all single-end/pair-end libraries together.

The workflow is automatic, but it's not fool-proof. It requires some decisions. The usage is shown using [KMC](https://github.com/refresh-bio/KMC) and [GenomeScope](https://github.com/schatzlab/genomescope). We provide also a real data tutorial on our wiki: [strawberry genome analysis](https://github.com/tbenavi1/smudgeplot/wiki/strawberry-tutorial). If you are interested in running Jellyfish instead of KMC, look at the [manual of smudgeplot with Jellyfish](https://github.com/tbenavi1/smudgeplot/wiki/manual-of-smudgeplot-with-jellyfish).

Give KMC all the files with trimmed reads to calculate kmer frequencies and then generate a histogram of kmers:

```
mkdir tmp
ls *.fastq.gz > FILES
# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES kmer_counts tmp
kmc_tools transform kmer_counts histogram kmer_k21.hist -cx10000
```

where `-k` is the kmer length, `-m` is the max amount of RAM to use in GB (1 to 1024), `-ci<value>` excludes kmers occurring less than \<value\> times, `-cs` is the maximum value of a counter, `FILES` is a file name with a list of input files, `kmer_counts` is the output file name prefix, `tmp` is a temporary directory, and `-cx<value>` is the maximum value of counter to be stored in the histogram file.

The next step is to extract genomic kmers using reasonable coverage thresholds. You can either inspect the kmer spectra and choose the L (lower) and U (upper) coverage thresholds via visual inspection, or you can estimate them using the script `kmer_cov_cutoff.R <kmer.hist> <L/U>`. Then, extract kmers in the coverage range from `L` to `U` using `kmc_dump`. Then give the dump of kmers to the python script `hetkmers.py` to compute the set of heterozygous kmers.

```
L=$(smudgeplot cutoff kmer_k21.hist L)
U=$(smudgeplot cutoff kmer_k21.hist U)
echo $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000
kmc_tools transform kmer_counts -ci$L -cx$U dump -s kmer_k21.dump
smudgeplot hetkmers -k 21 -o kmer_pairs < kmer_k21.dump
```

After using KMC, generate the smudgeplot using the coverages of the kmer pairs (`*_coverages_2.tsv` file). You can either supply the haploid kmer coverage (reported by GenomeScope) or let it be estimated directly from the data and compare it afterwards. If GenomeScope correctly identifies the peak of haploid kmers, the expected positions of the haplotype structures will overlap with high density smudges on the smudgeplot. If the overlap is not great you might consider adjusting the GenomeScope model and redoing the plot with a better estimate of the haploid coverage. Something like

```
smudgeplot plot -i kmer_pairs_coverages_2.tsv
```

will generate a basic smugeplot, the full usage of `smudgeplot.R` script is

TODO: fix following section
```
usage: smudgeplot.R [-h] [-v] [--homozygous] [-i INPUT]
                    [-o OUTPUT] [-t TITLE] [-q QUANTILE_FILT]
                    [-n N_COV] [-L LOW_CUTOFF] [-nbins NBINS]
                    [-k KMER_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print the version and exit
  --homozygous          Assume no heterozygosity in the genome - plotting a
                        paralog structure; [default FALSE]
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
  -q QUANTILE_FILT, --quantile_filt QUANTILE_FILT
                        Remove kmer pairs with coverage over the specified
                        quantile; [default none]
  -n N_COV, --n_cov N_COV
                        the haploid coverage of the sequencing data [default
                        inference from data]
  -L LOW_CUTOFF, --low_cutoff LOW_CUTOFF
                        the lower boundary used when dumping kmers [default
                        min(total_pair_cov) / 2]
  -nbins NBINS          the number of nbins used for smudgeplot matrix (nbins
                        x nbins) [default autodetection]
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        The kmer size used to calculate kmer spectra [default
                        21]
```

Smudgeplot generates two plots, one with coloration on a log scale and on a linear scale. The legend indicates approximate kmer pairs per tile densities. Note that a single polymorphism generates multiple heterozygous kmers. As such, the reported numbers do not directly correspond to the number of variants. Instead, the actual number is approximately 1/k times the reported numbers, where k is the kmer size (in summary already recalculated). It's important to note that this process does not exhaustively attempt to find all of the heterozygous kmers from the genome. Instead, only a sufficient sample is obtained in order to identify relative genome structure. You can also report the minimal number of loci that are heterozygous if the inference is correct.

### GenomeScope for estimation of L/U

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use the [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmer_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

This script estimates the size, heterozygosity, and repetitive fraction of the genome. By inspecting the fitted model you can determine the location of the smallest peak after the error tail. Then, you can decide the low end cutoff below which all kmers will be discarded as errors (cca 1/2 of the haploid kmer coverage), and the high end cutoff above which all kmers will be discarded (cca 8x of the haploid kmer coverage).

## Frequently Asked Questions

Are collected on [our wiki](https://github.com/tbenavi1/smudgeplot/wiki/FAQ). If you won't find an answer for your question there, open an [issue](https://github.com/tbenavi1/smudgeplot/issues/new/choose) or drop us an email.

## Obvious near future

One important potential feature of smudgeplots (not implemented yet) is to quantify errors in a genome assembly. The arch nemesis of variant calling is missmapping mostly because of flaws of a genome assembly - collapsed paralogs or separately assembled alleles. Smudgeplot could be potentially used to quantify these mistakes.

For example, imagine you have a diploid organism and a supposedly haploid assembly. If you extract kmer pairs with inferred four genomic copies, you would expect to find them exactly twice (because the haplotypes are expected to be collapsed). Therefore you can quantify the proportion of four-copy heterozygous kmers was assembled in corrected copy number. Analogically we can extract heterozygous kmers in two copies (those that should be collapsed into one sequence in the assembly) and map them to genome. Those that will be identified two copies are separately assembled alleles, those with only one matching assembly region are those assembled correctly.

There will be always some errors in a genome assembly, the best we can do is to try to quantify it!

## Computational requirements

The memory required scale linearly with the number of kmers and it is approximately 15x higher than the size of the dump file
(for 20Gb dump file you will need approx ~250Gb of RAM). Alternatively, you can estimate requerements for RAM by number of dumped kmers. It's approximately 350x higher than number of kmers in the dump file. If your file has too many kmers you can decrease computational requirement by reruning the kmer spectra with a smaller kmer size or by more strict filtering of the dumped kmers (higher L and smaller U).

We have not calculated the complexity of the algorithm yet. Usually for smaller genomes (<250Gb) it's couple of hours, the longest computation took bit more than one day.

The biggest genome we analyzed so far was a triplod genome with a haploid size 3.5Gbp. We have processed 1.5e9 genomic kmers and it have required 520GB of memory and two days of computation on eight cores.

## Contributions

The file [DEVELOPMENT.md](playground/DEVELOPMENT.md) contains some development notes. The directory [playground](playground) contains some snippets, attempts, and other items of interest.

This is definitely an open project, but we have not decided about how to recognize contributions yet.

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison has largely inspired the visual output of smudgeplots. The colourblind friendly colour theme was suggested by @ggrimes. Grateful for helpful comments of beta testers and pre-release chatters!

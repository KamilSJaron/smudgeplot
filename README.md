# Smudgeplot 

<font size ="4">**_Version: 0.5.3 Skylight_**</font>

<font size ="4">**_Authors: Sam Ebdon, [Gene W Myers](https://github.com/thegenemyers) and [Kamil S. Jaron](https://github.com/KamilSJaron), Tianyi Ma._**</font>

## Installation
 
This version of smudgeplot operates on FastK k-mer databases. The smudgeplot installation consists of a python package and C-backend to search for all the k-mer pairs (hetmers) and extract sequences of k-mer pairs (extract_kmer_pairs).

We recommend installing smudgeplot within a [conda](https://conda-forge.org/download/) environment.

```
#optional conda environment setup
conda create -n smudgeplot && conda activate smudgeplot
conda install pip

# install via pypi
pip install smudgeplot

# or download and install directly. See below if you need to compile the C dependencies.
git clone https://github.com/KamilSJaron/smudgeplot.git
cd smudgeplot && pip install .
smudgeplot -h # check installation succeeded
```

That should do everything necesarry to make smudgeplot fully operational. You can run `smudgeplot --help` to see if it worked. If you activated a virtual environment prior to installation (either `conda` or any other) then smudgeplot installed within the environment.

Note the smudgeplot version downloadable from conda itself is not currently up to date.

### Compiling the C code

The process above should install everything including compilation of the C backend. If you need or would like to know how to compile the code yourself you can simply run

```
make
```

This will not, however, install the smudgeplot python package. 

### Pypi installation [EXPERIMENTAL] 

We are working on packaging smudgeplot for pypi. You are welcome to try installing from pypi if you are interested and please open an issue if you have problems. If it fails please follow the main instructions above to install for now.

```
pip install smudgeplot
```

## Example run on Saccharomyces data

Requires ~2.1GB of space and `FastK` and `smudgeplot` installed.

```bash
# download data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR326/001/SRR3265401/SRR3265401_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR326/001/SRR3265401/SRR3265401_2.fastq.gz

# sort them in a reasonable place
mkdir data/Scer
mv *fastq.gz data/Scer/

# run FastK to create a k-mer database
FastK -v -t4 -k31 -M16 -T4 data/Scer/SRR3265401_[12].fastq.gz -Ndata/Scer/FastK_Table

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot hetmers -L 12 -t 4 -o data/Scer/kmerpairs --verbose data/Scer/FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot all -o data/Scer/trial_run data/Scer/kmerpairs.smu

# check that bunch files are generated (3 pdfs; some summary tables and logs)
ls data/Scer/trial_run_*
```

The y-axis scaling is by default 100, one can spcify argument `ylim` to scale it differently

```bash
smudgeplot all -o data/Scer/trial_run_ylim70 data/Scer/kmerpairs.smu -ylim 70
```

There is also a plotting module that requires the coverage, a list of smudges, and the smudge sizes listed in a tabular file. This plotting module does not perform the inference and should be used only if you know the right answers already. 

## How smudgeplot works

This tool extracts heterozygous kmer pairs from kmer count databases and performs gymnastics with them. We are able to disentangle genome structure by comparing the sum of kmer pair coverages (CovA + CovB) to their relative coverage (CovB / (CovA + CovB)). Such an approach also allows us to analyze obscure genomes with duplications, various ploidy levels, etc.

Smudgeplots are computed from raw or even better from trimmed reads and show the haplotype structure using heterozygous kmer pairs. For example (of an older version):

![smudgeexample](https://user-images.githubusercontent.com/8181573/45959760-f1032d00-c01a-11e8-8576-ff0512c33da9.png)

Every haplotype structure has a unique smudge on the graph and the heat of the smudge indicates how frequently the haplotype structure is represented in the genome compared to the other structures. The image above is an ideal case, where the sequencing coverage is sufficient to beautifully separate all the smudges, providing very strong and clear evidence of triploidy.

This tool is planned to be a part of [GenomeScope](https://github.com/tbenavi1/genomescope2.0) in the near future.

### More usage information

The input is a set of whole genome sequencing reads, the higher the coverage the better. The method is designed to process big datasets, don't hesitate to pull all single-end/pair-end libraries together.

The workflow is automatic, but it's not fool-proof. It requires some decisions. Use this tool jointly with [GenomeScope](https://github.com/tbenavi1/genomescope2.0).

Smudgeplot generates two plots, one with coloration on a log scale and the other on a linear scale. The legend indicates approximate kmer pairs per tile densities. Note that a single polymorphism generates multiple heterozygous kmers. As such, the reported numbers do not directly correspond to the number of variants. Instead, the actual number is approximately 1/k times the reported numbers, where k is the kmer size (in summary already recalculated). It's important to note that this process does not exhaustively attempt to find all of the heterozygous kmers from the genome. Instead, only a sufficient sample is obtained in order to identify relative genome structure. You can also report the minimal number of loci that are heterozygous if the inference is correct.

### GenomeScope

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope](https://github.com/tbenavi1/genomescope2.0) or use the [web server](http://genomescope.org/genomescope2.0))

This tool estimates the size, heterozygosity, and repetitive fraction of the genome. By inspecting the fitted model you can determine the location of the smallest peak after the error tail. Then, you can decide the low end cutoff below which all kmers will be discarded as errors (cca 0.5 times the haploid kmer coverage), and the high end cutoff above which all kmers will be discarded (cca 8.5 times the haploid kmer coverage).

## Frequently Asked Questions

Are collected on [our wiki](https://github.com/KamilSJaron/smudgeplot/wiki/FAQ). Smudgeplot does not demand much computational resources, but make sure you check [memory requirements](https://github.com/KamilSJaron/smudgeplot/wiki/smudgeplot-hetkmers#memory-requirements) before you extract kmer pairs (`hetkmers` task). If you don't find an answer for your question in FAQ, open an [issue](https://github.com/KamilSJaron/smudgeplot/issues/new/choose) or drop us an email.

Check [projects](https://github.com/KamilSJaron/smudgeplot/projects) to see how the development goes.

## Contributions

This is definitely an open project, contributions are welcome. You can check some of the ideas for the future in [projects](https://github.com/KamilSJaron/smudgeplot/projects) and in the development [dev](https://github.com/KamilSJaron/smudgeplot/tree/dev) branch. The file [playground/DEVELOPMENT.md](playground/DEVELOPMENT.md) contains some development notes. The directory [playground](playground) contains some snippets, attempts, and other items of interest.

## Reference

Ranallo-Benavidez, T.R., Jaron, K.S. & Schatz, M.C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* **11**, 1432 (2020). https://doi.org/10.1038/s41467-020-14998-3

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison has largely inspired the visual output of smudgeplots. The colourblind friendly colour theme was suggested by @ggrimes. Grateful for helpful comments of beta testers and pre-release chatters!

## Changelog

#### 0.5

 + experimental feature to extract sequences of the kmers in the pair; this functionality will hopefully restore at some point together with functionality to assess the quality of assembly.
 + histograms are back

#### 0.4

 + the search for the kmer pair will be within both canonical and non-canonical k-mer sets (Gene demonstrated it makes a difference)
 + the tool will be supporting FastK kmer counter only
 + the backend by Gene is paralelized and massively faster
 + the intermediate file will be a flat file with the 2d histogram with cov1, cov2, freq columns (as opposed to list of coverages of pairs cov1 cov2);
 + completelly revamped plot showing how all individual kmer pairs insead of agregating them into squares
 + new smudge detection algorithm based on grid projection on the smudge plane (working, but under revisions at the moment)
 + R package smudgeplot was retired and is no longer used

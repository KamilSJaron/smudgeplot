# Smudgeplot 

<font size ="4">**_Version: 0.4.0 Arched_**</font>

<font size ="4">**_Authors: [Gene W Myers](https://github.com/thegenemyers) and [Kamil S. Jaron](https://github.com/KamilSJaron), Tianyi Ma._**</font>

### Install the whole thing
 
This version of smudgeplot operates on FastK k-mer databases. So, before installing smudgeplot, please install [FastK](https://github.com/thegenemyers/FASTK). The smudgeplot installation consist of one python, two R scripts and the C-backend to search for all the k-mer pairs (hetmers) that needs to be compilet.

#### Quick

Assuming you have admin right / can write to `/usr/local/bin`, you can simply run

```bash
sudo make
```
That should do everything necesarry to make smudgeplot fully operational. You can run `smudgeplot.py --help` to see if it worked.

#### Custom installation location

If there is a different directory where you store your executables, you can specify `INSTALL_PREFIX` variable to make. The binaries are then added to `$INSTALL_PREFIX/bin`. For example

```bash
make -s INSTALL_PREFIX=~
```

will install smudgeplot to `~/bin/`.

#### Manual installation

Compiling the `C` executable

```
make exec/hetmers # this will compile hetmers (kmer pair searching engine of PloidyPlot) backend
```

Now you can move all three files from the `exec` directory somewhere your system will see it (or alternativelly, you can add that directory to `$PATH` variable).

```
install -C exec/smudgeplot.py /usr/local/bin
install -C exec/hetmers /usr/local/bin
install -C exec/smudgeplot_plot.R /usr/local/bin
install -C exec/centrality_plot.R /usr/local/bin
```

### Runing this version on Sacharomyces data
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
smudgeplot.py hetmers -L 12 -t 4 -o data/Scer/kmerpairs --verbose data/Scer/FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o data/Scer/trial_run data/Scer/kmerpairs_text.smu

# check that bunch files are generated (3 pdfs; some summary tables and logs)
ls data/Scer/trial_run_*
```

The y-axis scaling is by default 100, one can spcify argument `ylim` to scale it differently

```bash
smudgeplot.py all -o data/Scer/trial_run_ylim70 data/Scer/kmerpairs_text.smu -ylim 70
```

There is also a plotting module that requires the coverage and a list of smudges and their respective sizes listed in a tabular file. This plotting module does not inference and should be used only if you know the right answers already. 

### How smudgeplot works

This tool extracts heterozygous kmer pairs from kmer count databases and performs gymnastics with them. We are able to disentangle genome structure by comparing the sum of kmer pair coverages (CovA + CovB) to their relative coverage (CovB / (CovA + CovB)). Such an approach also allows us to analyze obscure genomes with duplications, various ploidy levels, etc.

Smudgeplots are computed from raw or even better from trimmed reads and show the haplotype structure using heterozygous kmer pairs. For example (of an older version):

![smudgeexample](https://user-images.githubusercontent.com/8181573/45959760-f1032d00-c01a-11e8-8576-ff0512c33da9.png)

Every haplotype structure has a unique smudge on the graph and the heat of the smudge indicates how frequently the haplotype structure is represented in the genome compared to the other structures. The image above is an ideal case, where the sequencing coverage is sufficient to beautifully separate all the smudges, providing very strong and clear evidence of triploidy.

This tool is planned to be a part of [GenomeScope](https://github.com/tbenavi1/genomescope2.0) in the near future.

### More about the use

The input is a set of whole genome sequencing reads, the more coverage the better. The method is designed to process big datasets, don't hesitate to pull all single-end/pair-end libraries together.

The workflow is automatic, but it's not fool-proof. It requires some decisions. Use this tool joinlty with [GenomeScope](https://github.com/tbenavi1/genomescope2.0). The tutorials on our wiki are currently outdated (build for version 0.2.5), and will be updated by 18th of October. 

Smudgeplot generates two plots, one with coloration on a log scale and the other on a linear scale. The legend indicates approximate kmer pairs per tile densities. Note that a single polymorphism generates multiple heterozygous kmers. As such, the reported numbers do not directly correspond to the number of variants. Instead, the actual number is approximately 1/k times the reported numbers, where k is the kmer size (in summary already recalculated). It's important to note that this process does not exhaustively attempt to find all of the heterozygous kmers from the genome. Instead, only a sufficient sample is obtained in order to identify relative genome structure. You can also report the minimal number of loci that are heterozygous if the inference is correct.

### GenomeScope

You can feed the kmer coverage histogram to GenomeScope. (Either run the [genomescope script](https://github.com/schatzlab/genomescope/blob/master/genomescope.R) or use the [web server](http://qb.cshl.edu/genomescope/))

```
Rscript genomescope.R kmcdb_k21.hist <k-mer_length> <read_length> <output_dir> [kmer_max] [verbose]
```

This script estimates the size, heterozygosity, and repetitive fraction of the genome. By inspecting the fitted model you can determine the location of the smallest peak after the error tail. Then, you can decide the low end cutoff below which all kmers will be discarded as errors (cca 0.5 times the haploid kmer coverage), and the high end cutoff above which all kmers will be discarded (cca 8.5 times the haploid kmer coverage).

## Frequently Asked Questions

Are collected on [our wiki](https://github.com/KamilSJaron/smudgeplot/wiki/FAQ). Smudgeplot does not demand much on computational resources, but make sure you check [memory requirements](https://github.com/KamilSJaron/smudgeplot/wiki/smudgeplot-hetkmers#memory-requirements) before you extract kmer pairs (`hetkmers` task). If you don't find an answer for your question in FAQ, open an [issue](https://github.com/KamilSJaron/smudgeplot/issues/new/choose) or drop us an email.

Check [projects](https://github.com/KamilSJaron/smudgeplot/projects) to see how the development goes.

## Contributions

This is definitely an open project, contributions are welcome. You can check some of the ideas for the future in [projects](https://github.com/KamilSJaron/smudgeplot/projects) and in the development [dev](https://github.com/KamilSJaron/smudgeplot/tree/dev) branch. The file [playground/DEVELOPMENT.md](playground/DEVELOPMENT.md) contains some development notes. The directory [playground](playground) contains some snippets, attempts, and other items of interest.

## Reference

Ranallo-Benavidez, T.R., Jaron, K.S. & Schatz, M.C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* **11**, 1432 (2020). https://doi.org/10.1038/s41467-020-14998-3

## Acknowledgements

This [blogpost](http://www.everydayanalytics.ca/2014/09/5-ways-to-do-2d-histograms-in-r.html) by Myles Harrison has largely inspired the visual output of smudgeplots. The colourblind friendly colour theme was suggested by @ggrimes. Grateful for helpful comments of beta testers and pre-release chatters!

## Smudgeplot

Have you ever sequenced something not-well studied? Something that might show strange genomic signatures? Smudgeplot is a visualisation technique for whole-genome sequencing reads from a single individual. The visualisation techique is based on the idea of het-mers. Het-mers are k-mer pairs that are exactly one nucleotide pair away from each other, while forming a unique pair in the sequencing dataset. These k-mers are assumed to be mostly representing two alleles of a heterozygous, but potentially can also show pairing of imperfect paralogs, or sequencing errors paired up with a homozygous genomic k-mer. Nevertheless, the predicted ploidy by smudgeplot is simply the ploidy with the highest number of k-mer pairs (if a reasonable estimate must be evaluated for each individual case!).



### Installing the software


```
mkdir src bin && cd src # create directories for source code & binaries
git clone -b sploidyplot https://github.com/KamilSJaron/smudgeplot
git clone https://github.com/thegenemyers/FastK

cd smudgeplot && make -s INSTALL_PREFIX=/workspace && cd ..
cd FastK && make
ln -s `pwd`/FastK /workspace/bin/

export PATH="/workspace/bin:"$PATH
```

R packages

```
install.packages('argparse')
install.packages('devtools')

```

Finally, if your session is running; Start downloading the data;

### Constructing a database

The whole process operates with raw, or trimmed sequencing reads. From those we generate a k-mer database using [FastK](https://github.com/thegenemyers/FASTK). FastK is currently the fastest k-mer counter out there and the only supported by the lastest version of smudgeplot*. This database contains an index of all the k-mers and their coverages in the sequencing readset. Within this set user must chose a theshold for excluding low frequencing k-mers that will be considered errors. That choice is not too difficult to make by looking at the k-mer spectra. Of all the retained k-mers we find all the het-mers. Then we plot a 2d histogram. 



*Note: The previous versions of smudgeplot (up to 2.5.0) were operating on k-mer "dumps" flat files you can generate with any counter you like. You can imagine that text files are very inefficient to operate on. The new version is operating directly on the optimised k-mer database instead.

```
FastK
```

### Constructing a database

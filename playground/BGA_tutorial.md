## Smudgeplot

Have you ever sequenced something not-well studied? Something that might show strange genomic signatures? Smudgeplot is a visualisation technique for whole-genome sequencing reads from a single individual. The visualisation techique is based on the idea of het-mers. Het-mers are k-mer pairs that are exactly one nucleotide pair away from each other, while forming a unique pair in the sequencing dataset. These k-mers are assumed to be mostly representing two alleles of a heterozygous, but potentially can also show pairing of imperfect paralogs, or sequencing errors paired up with a homozygous genomic k-mer. Nevertheless, the predicted ploidy by smudgeplot is simply the ploidy with the highest number of k-mer pairs (if a reasonable estimate must be evaluated for each individual case!).



### Installing the software

Open gitpod. And install the development version of smudgeplot (branch sploidyplot) & FastK. 

```
mkdir src bin && cd src # create directories for source code & binaries
git clone -b sploidyplot https://github.com/KamilSJaron/smudgeplot
git clone https://github.com/thegenemyers/FastK
```

Now smudgeplot make install smudgeplot R package, compiles the C kernel for searching for k-mer pairs and copy all the executables to `workspace/bin/` (which will be our dedicated spot for executables). 

```
cd smudgeplot && make -s INSTALL_PREFIX=/workspace && cd ..
cd FastK && make FastK Histex
install -c Histex FastK /workspace/bin/
```


** 8 Datasets **

Species name	SRA/ENA ID
Pseudoloma neurophilia	SRR926312
Tubulinosema ratisbonensis	ERR3154977
Nosema ceranae	SRR17317293
Nematocida ausubeli	SRR058692
Nematocida ausubeli	SRR350188
Hamiltosporidium tvaerminnensis	SRR16954898
Encephalitozoon hellem	SRR14017862
Agmasoma penaei	SRR926341

TODO: get them urls

Finally, if your session is running; Start downloading the data; For example:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR926/SRR926341/SRR926341_[12].fastq.gz
```

### Constructing a database

The whole process operates with raw, or trimmed sequencing reads. From those we generate a k-mer database using [FastK](https://github.com/thegenemyers/FASTK). FastK is currently the fastest k-mer counter out there and the only supported by the lastest version of smudgeplot*. This database contains an index of all the k-mers and their coverages in the sequencing readset. Within this set user must chose a theshold for excluding low frequencing k-mers that will be considered errors. That choice is not too difficult to make by looking at the k-mer spectra. Of all the retained k-mers we find all the het-mers. Then we plot a 2d histogram. 


*Note: The previous versions of smudgeplot (up to 2.5.0) were operating on k-mer "dumps" flat files you can generate with any counter you like. You can imagine that text files are very inefficient to operate on. The new version is operating directly on the optimised k-mer database instead.

```
FastK -v -t4 -k31 -M16 -T4 *.fastq.gz -NSRR8495097
```

20'

23:24

one file is also ~20', it's mostly function of the number of k-mers, we could speed it up by chosing higher t maybe?

### Getting the k-mer spectra out of it

```
Histex -G SRR8495097 > SRR8495097_k31.hist

| GeneScopeFK -o data/Pvir1/GenomeScopeFK/ -k 17


# Run PloidyPlot to find all k-mer pairs in the dataset
PloidyPlot -e12 -k -v -T4 -odata/Scer/kmerpairs data/Scer/FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py plot -n 15 -t Sacharomyces -o data/Scer/trial_run data/Scer/kmerpairs_text.smu 
```

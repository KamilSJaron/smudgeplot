### Manual tests

This is a place for manual tests that have not been automated (yet?).
Don't forget to re-install the package/script before execution. Somehting like

```
Rscript install.R
install -C exec/smudgeplot /usr/local/bin
install -C exec/smudgeplot_plot.R /usr/local/bin
```

should do the job.

#### interface tests

##### smudgeplot plot

```
./exec/smudgeplot plot tests/data/toy_middle_coverages.tsv -o tests/data/toy_middle_R -n 100 -nbins 20
```

generates a smudgeplot from the toy data.

##### smudgeplot hetkmers

two different methods to extract homologous kmers

###### All

```
./exec/smudgeplot hetkmers -o tests/data/toy_all tests/data/toy_kmer_k21.dump
```

###### Middle

```
./exec/smudgeplot hetkmers -o tests/data/toy_middle tests/data/toy_kmer_k21.dump --middle
```

generates two files `tests/data/toy_middle_coverages.tsv` and `tests/data/toy_middle_sequences.tsv` with coverages and sequences of kmer pairs.


##### smudgeplot kmer extraction

```
exec/smudgeplot.py extract -cov tests/data/toy_all_coverages.tsv -seq tests/data/toy_all_sequences.tsv   -minc 400 -maxc 410 -minr 0.49 -maxr 0.5 | head
```

would extract 80 lines made of 20 kmers pairs, but showing only the first ten.

##### Dicots

This is a large dataset of the first 540 dicot genomes sequenced by the Tree of Life. Some of them are completed, some of them are with insufficient coverage or otherwise qc failed data. The idea here is to be able to tell those apart, get reasonable defaults so the generated plot is meaningful for a reasonable number (i.e. nearly all) of them.

```
time ./exec/smudgeplot.py plot data/dicots/smu_files/daAchMill1.k31_ploidy.smu.txt -o data/dicots/alt_plots/daAchMill1 --alt_plot -q 0.9
```

```bash
for smu in data/dicots/smu_files/*.smu.txt; do
    species=$(basename $smu)
    echo $species $smu
    time ./exec/smudgeplot.py plot $smu -o data/dicots/alt_plots/$species --alt_plot -q 0.9
done

for smu in data/dicots/smu_files/*.smu.txt; do
    species=$(basename $smu .smu.txt)
    echo $species $smu
    time ./exec/smudgeplot.py plot $smu -c 10 -o data/dicots/alt_plots_c10/$species --alt_plot -q 0.9
done
```
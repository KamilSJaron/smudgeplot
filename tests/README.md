### Manual tests

This is a place for manual tests that have not been automated (yet?).
Don't forget to re-install the package/script before execution. Somehting like

```
make install INSTALL_PREFIX=~
```

should do the job.

#### interface tests

##### data prep

Download `SRR3265401` - nice teteraploid Sacharomyces run I use often for testing.

##### smudgeplot plot

Defaults:

```
smudgeplot.py all data/Scer/kmerpairs_default_e2_text.smu -o data/Scer/240918_trial
```

Testing parameters:

```
smudgeplot.py all data/Scer/kmerpairs_default_e2_text.smu -o data/Scer/240918_trial_params -t "Species 1" -c 20 -ylim 80 -col_ramp magma --invert_cols
```

##### smudgeplot hetkmers

two different methods to extract homologous kmers

TODO

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

for smu in $(ls data/dicots/smu_files/*.smu.txt | head -20); do
    species=$(basename $smu .smu.txt)
    echo $species $smu
    smudgeplot.py all $smu -t $species -o data/dicots/automated_smudgeplots/$species
done
```
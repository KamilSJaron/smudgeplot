# Testing

There are three levels we can test on. We can test if the installation procedure works correctly, we can test individual component of the

## Instalation

not tested now.

## Units

Tests of functionalities of classes.

When a need for testing units appear, we have `unittest` environment set up. I don't have a capacity to write tests for everything, but it's useful to write new tests for individual cases that will be latter identified as problematic.

These tests can be executed by

```
make test
```

## Integration tests

Tests of the scripts. Simulating use cases of the software. I hope to hook these to Travis at some point. Right now I will just make here a collection of commands that will do the testing.

Generate a toy dumpfile `tests/data/toy_kmer_k21.dump`.

```
python3 tests/generatedDumpfileForTesting.py
```

there is a bit of random generation involved, but very little, should not have any consequences for the downstream tests (adding bit of sequencing noise to make pictures nicer).

### API of cutoff

### API of smudgeplot

Runs, but wrong inference:

```
smudgeplot plot tests/data/toy_middle_coverages.tsv -t "fake\ data" -nbins 20 -o tests/data/toy_middle_plain
```

```
smudgeplot plot tests/data/toy_middle_coverages.tsv -t "fake\ data" -nbins 20 -o tests/data/toy_middle_plain -n 100
```

### API of hetkmers

Works:

```
cat data/Minc1/kmer_k21.dump | smudgeplot hetkmers -o data/Minc_test/files -k 21
```

```
smudgeplot hetkmers -o tests/data/toy_middle -k 21 tests/data/toy_kmer_k21.dump
```

Does not:

```
smudgeplot hetkmers -o tests/data/minc_sample_all -k 21 tests/data/minc_k21_sample.dump --all
```

```
smudgeplot hetkmers -o tests/data/toy_middle -k 21 tests/data/toy_kmer_k21.dump --all
```

### Virtual environment testing

Testing of virtual env if it works

```
virtualenv -p python3 venv
source venv/bin/activate
pip3 install -r requirements.txt
pip3 install PySAIS==1.0.4
python3 setup.py install
smudgeplot --version
# heck yeah
rm -r venv
```

### Test that is failing

```
python3 tests/generatedDumpfileForTesting.py
smudgeplot hetkmers -o tests/data/toy_middle -k 21 tests/data/toy_kmer_k21.dump
python3 tests/kmers_in_genome.py
# line 4
```

Basically, I made this toy dataset, where I know exactly what variants I am implanting in. And one of the heterozygous loci is not identified although they represent a unique pair. Look at `tests/kmers_in_genome.py` for more details.
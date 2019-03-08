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


### API of cutoff

### API of smudgeplot

### API of hetkmers

```
cat data/Minc1/kmer_k21.dump | smudgeplot hetkmers -o data/Minc_test/files -k 21
```
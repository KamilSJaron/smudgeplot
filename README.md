# Sibling kmer analysis

This tool extract heterozygous kmer pairs from jellyfish dump files and does a gymnastics with it. We are able to disentangle genome stricture by comparison of the sum of kmer pair coverages and their relative coverage. Such approach allows to analyze also obscure genomes with duplications, various ploidy levels and stuff. it's work in progress but it's going to be great.

Content :

- [hetkmers.py](hetkmers.py): python script to compute heterozygous kmer pairs
- [blurplot.R](blurplot.R): R script that creates so-called blur plot from `coverages_2.tsv` file - a two column file with kmer sorted kmer coverages.

- [playground.py](playground.py): snippets of python code for testing that should be maintained.
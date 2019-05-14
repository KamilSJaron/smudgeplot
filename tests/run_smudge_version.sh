#!/bin/bash

version=$(smudgeplot -v | cut -f 3 -d ' ')
echo "testing $version"

outdir=figures/$version
mkdir -p $outdir
rm $outdir/*

Rscript install.R
install -C exec/smudgeplot /usr/local/bin
install -C exec/smudgeplot_plot.R /usr/local/bin

for sp in "Aric1" "Avag1" "Mflo2" "Rmag1"; do
	smudgeplot plot -o $outdir/"$sp"_smudgeplot_"$version" -t "$sp $version" data/$sp/*coverages_2.tsv
done

for sp in "Ps791" "Rvar1"; do
	smudgeplot plot -o $outdir/"$sp"_smudgeplot_"$version" -t "$sp $version" data/$sp/*coverages_2.tsv -nbins 15
done

sp="Lcla1"

smudgeplot plot -o $outdir/"$sp"_smudgeplot_"$version"_homozyg -t "$sp $version" --homozygous data/$sp/*coverages_2.tsv

smudgeplot plot -o $outdir/"$sp"_smudgeplot_"$version" -t "$sp $version" data/$sp/*coverages_2.tsv

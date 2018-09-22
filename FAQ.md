
##### 1. What are heterozygous kmers? Are these kmer pairs with a close but not perfect sequence match?


Yes. Right now heterozygous kmers are thost that:
- are exactly one SNP from each other ( for instance ATGTCA ATGTTA)
- form a unique pair (i.e. there are no other kmers one SNP away from them. for instance ATGTCA ATGTTA ATGTTT would be discarded - they three not two)

Like this we heavily subsample the genome, but so far this was very sufficient to sample enough heterozygous kmers to see the genome structure.


##### 2. Are you equating kmer with haplotype? If not, how do you infer haplotypes?


I assume that at least some of the kmers will be heterozygous - i.e. one of the kmer is from one haplotype and the other kmers is the other. If the ploidy is higher, like triploid, we expect to find locations that have the same kmers in two of the haplotypes one one SNP difference in the third haplotype. If the genome is tetraploid we will detect some pairs where two haplotypes are similar and two haplotypes are too diverged (looking like AB), some cases where one haplytype will be diverged but the three other will be "triploid like" and therefore it will look like AAB and finally majority of heterozygosity will be carried by either AABB or AAAB structures which will also tell us whether the genome structure is AA'BB' or AA'A''B (i.e. what is the branching of haplotypes)

##### 3. In the example plot, I suspect that the red dot at x = 1/3 indicates triploidy, but only if two of the subgenomes are more similar to each other relative to the third (which I presume is meant by the AAB label).

Yes, because we search for unique pairs, if the pair comes from three haplotypes it must be twice same (A) and once one SNP away (B).


##### 4. If the three genomes of a triploid organism were equidistant to each other, then that blob would move to and fuse with the one labelled 'AB' at x = 0.5 because you are considering pairs. Is that correct?

This could happen, but under slightly different circumstances. Only if one haplotype would be very diverged and the genome would be ABC (where AB are similar, C is distant), then we would not be able to identify C as the corresponding haplotype and smudgeplot would indeed look like diploid (because the only kmer pairs one SNP away from each other would be kmers heterozygous kmers between A and B). If all thee haplotypes are evenly distant but still close than AAB smudge will be mixture of genomic positions where 1. 2. haplotypes are the same and 3. has a SNP, 1. 3. are the same and 2. is different and 2. 3. are same and 1. is different. Smudgeplots are not phasing haplotypes, so we can not make strong claims about haplotype divergence for triploids (although I am trying to make some guesses).

##### 5. What are the assumptions of smudgeplots?

We require:
- single individual (or clonal population)
- not too high coverage variance (this can be compensated with genome coverage, but for instance a 100x genome with whole genome amplification was not enough)
- low sequencing error rates (<1%)

However, the method should be robust to:
- contamination (to cetain extend)
- genome subsampling (it should work on exon capture etc)
- quite some coverage variance if the coverage is sufficient (I see a signal in tardigrades that was a messy sequencing too)
- presence of adapters (because they are high frequency kmers, but this was not verified yet)

And mainly we do NOT use reference genome, de-novo genome assembly, mapping or base quality scores for the inference. Therefore we are free of all the downstream biases that are coming from the steps above.

##### 6. Can I use mate-pairs?

Theoretically you can. If you have sufficient coverage from pair-end / single end reads (>50x), I would probably try without mate pairs first.

The reason is that mate pair library is constructed using circularization of fragments using linker adapters and subsequent fragmentation and pair-end sequencing.
This process is way more complicated than simple pair-end sequencing and results in lot of adapters sequenced. For the same reason mate-pairs are rarely used for genome assembly (they are rather used for scaffolding).

However, we have never tried it with mate pairs. Maybe my worry is not in place and maybe it will work just find. If you try, definitely let me know :-).

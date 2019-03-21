from Bio import SeqIO
from PySAIS import sais
from bisect import bisect_right
import numpy as np
import gzip
import logging

def evaluated2assembly_position(eval_pos):
    _scf = bisect_right(scf_sizes, eval_pos)
    eval_pos - scf_sizes[_scf - 1]
    if _scf > 0:
        return((scf_names[_scf], eval_pos - scf_sizes[_scf - 1]))
    else :
        return((scf_names[_scf], eval_pos))

def searchKmer(kmer):
    alignments = []
    L_kmer = kmer[0:10]
    l = 0
    r = len(sa) - 1
    while l <= r:
        m = (l + r) // 2
        eval_pos = sa[m]
        L_genome_kmer = genome[eval_pos:(eval_pos + 10)]
        if L_kmer == L_genome_kmer:
            r_L_converged = r
            for N in ['A', 'C', 'G', 'T']:
                kmer = kmer[0:10] + N + kmer[11:21]
                while l <= r:
                    m = (l + r) // 2
                    eval_pos = sa[m]
                    genome_kmer = genome[eval_pos:(eval_pos + 21)]
                    if kmer == genome_kmer:
                        m_hit = m
                        # following ms
                        while kmer == genome_kmer:
                            scf, pos = evaluated2assembly_position(eval_pos)
                            alignments.append([scf, pos, N])
                            m += 1
                            eval_pos = sa[m]
                            genome_kmer = genome[eval_pos:(eval_pos + 21)]
                        # previous ms
                        m = m_hit - 1
                        eval_pos = sa[m]
                        genome_kmer = genome[eval_pos:(eval_pos + 21)]
                        while kmer == genome_kmer:
                            scf, pos = evaluated2assembly_position(eval_pos)
                            alignments.append([scf, pos, N])
                            m -= 1
                            eval_pos = sa[m]
                            genome_kmer = genome[eval_pos:(eval_pos + 21)]
                        break
                    elif genome_kmer < kmer:
                        l = m + 1
                    else:
                        r = m - 1
                r = r_L_converged # recover the r after converging on the left sub-kmer
            break
        elif L_genome_kmer < L_kmer:
            l = m + 1
        else:
            r = m - 1
    return(alignments)

def map(args):
    # consider also zipped genomes
    # parser.add_argument('genomefile', type=argparse.FileType('r'), help='The reference genome.')
    kmer_genome_file = args.genomefile
    output_pattern = args.o
    # assume args.s == 0

    # LOAD GENOME; scf_names are names of scaffolds; sequences is a list of sequences
    logging.info('Loading genome')
    if kmer_genome_file[-2:] == "gz":
        ffile = SeqIO.parse(gzip.open(kmer_genome_file, "rt"), "fasta")
    else :
        ffile = SeqIO.parse(kmer_genome_file, "fasta")

    sequences = []
    scf_names = []
    for seq_record in ffile:
        scf_names.append(seq_record.name)
        sequences.append(str(seq_record.seq).upper())

    ffile.close()

    scf_sizes = np.cumsum([len(scf) + 1 for i,scf in enumerate(sequences)])
    # scf = bisect_right(scf_sizes, eval_pos) + 1

    logging.info('Building suffix array (might take several minutes)')
    genome = '$'.join(sequences)
    sa = sais(genome)

    logging.info('Loading kmers')

    # output_pattern -> kmer_files
    # kmer_file_name = 'data/Tps1/Tps_middle_pair_end_reads_kmers_in_smudge_2.txt'
    kmer_file_name = output_pattern + '_kmers_in_smudge_2.txt'
    with open(kmer_file_name, 'r') as kmer_file:
        kmers = [kmer.rstrip() for kmer in kmer_file]

    logging.info('Mapping kmers')
    start = time. time()
    mapping_list = [searchKmer(kmer) for kmer in s1_kmers]
    end = time. time()
    print(end - start)

    logging.info('Generating stats')

    hit_numbers = [len(mapped_kmer) for mapped_kmer in mapping_list]
    hit_numbers.value_counts()
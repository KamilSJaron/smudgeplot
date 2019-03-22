from Bio import SeqIO
from PySAIS import sais
from bisect import bisect_right
import numpy as np
import gzip
import logging

class mapper:
    def __init__(self, user_args):
        self.genome_file_name = user_args.genomefile
        self.output_pattern = user_args.o
        self.to_plot = user_args.s

    def evaluated2assembly_position(self, eval_pos):
        scf_index = bisect_right(self.scf_sizes, eval_pos)
        eval_pos - self.scf_sizes[scf_index - 1]
        if scf_index > 0:
            return((self.scf_names[scf_index], eval_pos - self.scf_sizes[scf_index - 1]))
        else :
            return((self.scf_names[scf_index], eval_pos))

    def searchKmer(self, kmer):
        alignments = []
        L_kmer = kmer[0:10]
        l = 0
        r = len(self.sa) - 1
        while l <= r:
            m = (l + r) // 2
            eval_pos = self.sa[m]
            L_genome_kmer = self.genome[eval_pos:(eval_pos + 10)]
            if L_kmer == L_genome_kmer:
                r_L_converged = r
                for N in ['A', 'C', 'G', 'T']:
                    kmer = kmer[0:10] + N + kmer[11:21]
                    while l <= r:
                        m = (l + r) // 2
                        eval_pos = self.sa[m]
                        genome_kmer = self.genome[eval_pos:(eval_pos + 21)]
                        if kmer == genome_kmer:
                            m_hit = m
                            # following ms
                            while kmer == genome_kmer:
                                scf, pos = self.evaluated2assembly_position(eval_pos)
                                alignments.append([scf, pos, N])
                                m += 1
                                eval_pos = self.sa[m]
                                genome_kmer = self.genome[eval_pos:(eval_pos + 21)]
                            # previous ms
                            m = m_hit - 1
                            eval_pos = self.sa[m]
                            genome_kmer = self.genome[eval_pos:(eval_pos + 21)]
                            while kmer == genome_kmer:
                                scf, pos = self.evaluated2assembly_position(eval_pos)
                                alignments.append([scf, pos, N])
                                m -= 1
                                eval_pos = self.sa[m]
                                genome_kmer = self.genome[eval_pos:(eval_pos + 21)]
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

    def loadGenome(self):
        # LOAD GENOME; scf_names are names of scaffolds; sequences is a list of sequences
        logging.info('Loading genome')
        # autodetect gzipped files
        if self.genome_file_name[-2:] == "gz":
            ffile = SeqIO.parse(gzip.open(self.genome_file_name, "rt"), "fasta")
        else :
            ffile = SeqIO.parse(self.genome_file_name, "fasta")
        # load scaffold names and sequences
        self.scf_names = []
        self.sequences = []
        for seq_record in ffile:
            self.scf_names.append(seq_record.name)
            self.sequences.append(str(seq_record.seq).upper())
        ffile.close()
        # calculate scaffold sizes
        self.scf_sizes = np.cumsum([len(scf) + 1 for scf in self.sequences])
        # create a string representing the gneome
        self.genome = '$'.join(self.sequences)

    def constructSuffixArray(self):
        logging.info('Building suffix array (might take several minutes)')
        self.sa = sais(self.genome)

    def map(self, args):

        logging.info('Loading kmers')

        # output_pattern -> kmer_files
        # kmer_file_name = 'data/Tps1/Tps_middle_pair_end_reads_kmers_in_smudge_2.txt'
        kmer_file_name = self.output_pattern + '_kmers_in_smudge_1.txt'
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
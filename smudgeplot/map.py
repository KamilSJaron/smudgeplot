from Bio import SeqIO
from Bio import bgzf
from Bio.Seq import Seq
from PySAIS import sais
from bisect import bisect_right
from time import time
import numpy as np
import gzip
import logging
import csv
import os
import time

class mapper:
    def __init__(self, user_args):
        self.genome_file_name = user_args.genomefile
        self.output_pattern = user_args.o
        self.to_plot = user_args.s

    def evaluated2assembly_position(self, eval_pos):
        scf_index = bisect_right(self.scf_cum_sizes, eval_pos)
        eval_pos - self.scf_cum_sizes[scf_index - 1]
        if scf_index > 0:
            return((self.scf_names[scf_index], eval_pos - self.scf_cum_sizes[scf_index - 1]))
        else :
            return((self.scf_names[scf_index], eval_pos))

    def searchKmer(self, kmer):
        reverse_complementary_kmer = str(Seq(kmer).reverse_complement())
        alignments = []
        for strand in ['+', '-']:
            if strand == '-':
                kmer = reverse_complementary_kmer
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
                                    alignments.append([scf, pos, strand, N])
                                    m += 1
                                    eval_pos = self.sa[m]
                                    genome_kmer = self.genome[eval_pos:(eval_pos + 21)]
                                # previous ms
                                m = m_hit - 1
                                eval_pos = self.sa[m]
                                genome_kmer = self.genome[eval_pos:(eval_pos + 21)]
                                while kmer == genome_kmer:
                                    scf, pos = self.evaluated2assembly_position(eval_pos)
                                    alignments.append([scf, pos, strand, N])
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
        self.scf_sizes = [len(scf) for scf in self.sequences]
        self.scf_cum_sizes = np.cumsum([len(scf) + 1 for scf in self.sequences])
        # create a string representing the gneome
        self.genome = '$'.join(self.sequences)

    def constructSuffixArray(self):
        logging.info('Building suffix array (might take several minutes)')
        self.sa = sais(self.genome)

    def map(self):
        # if user specifies a smudge to process
        if self.to_plot == 0:
            proc_smudge = 1
        else:
            proc_smudge = self.to_plot

        while os.path.isfile(self.output_pattern + '_kmers_in_smudge_' + str(proc_smudge) + '.txt'):
            logging.info('Processig smudge ' + str(proc_smudge))

            logging.info('Loading kmers')
            kmer_file_name = self.output_pattern + '_kmers_in_smudge_' + str(proc_smudge) + '.txt'
            with open(kmer_file_name, 'r') as kmer_file:
                kmers = [kmer.rstrip() for kmer in kmer_file]

            logging.info('Mapping kmers')
            start = time.time()
            mapping_list = [self.searchKmer(kmer) for kmer in kmers]
            end = time.time()
            logging.info('Kmers mapped in ' + str(round(end - start, 4)) + ' seconds')

            # logging.info('Generating stats')
            # hit_numbers = [len(mapped_kmer) for mapped_kmer in mapping_list]
            # hit_numbers.value_counts()

            logging.info('Generating bam file.')
            hlaf_of_kmer = str(int((len(kmers[0]) - 1) / 2))
            cigar = hlaf_of_kmer + '=1X' + hlaf_of_kmer + '='
            bamfile_name = self.output_pattern + "_" + str(proc_smudge) + "_mapped.bam"
            with bgzf.BgzfWriter(bamfile_name, 'w') as bamfile:
                writer = csv.writer(bamfile, delimiter='\t', lineterminator='\n')
                writer.writerow(['@HD', 'VN:1.3', 'SO:coordinate'])
                for i, scf in enumerate(self.scf_names):
                    writer.writerow(['@SQ', 'SN:' + scf, 'LN:' + str(self.scf_sizes[i])])
                for i, mapped_kmer in enumerate(mapping_list):
                    for entry in mapped_kmer:
                        if entry[2] == '+':
                            flag = 0
                        else:
                            flag = 16
                        writer.writerow(['k' + str(i + 1), str(flag), entry[0], entry[1], 255, cigar, '*', '0', '0', kmers[i], '*'])
                    if not mapped_kmer:
                        flag = 4
                        writer.writerow(['k' + str(i + 1), str(flag), '*', 0, 255, '*', '*', '0', '0', kmers[i], '*'])
            logging.info('bam file saved.')
            if self.to_plot != 0:
                break
            proc_smudge += 1


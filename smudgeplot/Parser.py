import argparse
import sys
import os

class Parser():

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            # description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
            usage='''smudgeplot <task> [options]
tasks: hetkmers  Calculate unique kmer pairs from a Jellyfish or KMC dump file.
       plot      Generate 2d histogram; infere ploidy and plot a smudgeplot\n\n''')
        self.parser.add_argument('task', help='Task to execute; for task specific options execute smudgeplot <task> -h')
        self.parser.add_argument('-v', '--version', action="store_true", default = False, help="print the version and exit")
        for arg in sys.argv:
            print(arg)
        self.arguments = self.parser.parse_args([''.join(sys.argv[1])])
        if self.arguments.version:
            exit(1)
        if not hasattr(self, self.arguments.task):
            print('** Error: invalid task "' + self.arguments.task + '"')
            self.parser.print_usage()
            exit(1)
        getattr(self, self.arguments.task)()

    def hetkmers(self):
        '''
        Calculate unique kmer pairs from a Jellyfish or KMC dump file.
        '''
        parser = argparse.ArgumentParser(
            description='Calculate unique kmer pairs from a Jellyfish or KMC dump file.')
        parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Alphabetically sorted Jellyfish or KMC dump file (stdin).')
        parser.add_argument('-o', help='The pattern used to name the output (kmerpairs).', default='kmerpairs')
        parser.add_argument('-k', help='The length of the kmer.', default=21)
        parser.add_argument('-t', help='Number of processes to use.', default = 4)
        parser.add_argument('--all', dest='all', action='store_const', const = True, default = False,
                          help='Get all kmer pairs one SNP away from each other (default: just the middle one).')
        self.__assignArguments__(parser, sys.argv)

    def plot(self):
        '''
        Generate 2d histogram; infer ploidy and plot a smudgeplot.
        '''
        parser = argparse.ArgumentParser(description='Generate 2d histogram for smudgeplot')
        parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='name of the input tsv file with covarages [default \"coverages_2.tsv\"]"')
        parser.add_argument('-o', help='The pattern used to name the output (smudgeplot).', default='smudgeplot')
        parser.add_argument('-q', help='Remove kmer pairs with coverage over the specified quantile; [default none]', default=1)
        parser.add_argument('-L', help='The lower boundary used when dumping kmers [default min(total_pair_cov) / 2]', type=int, default=0)
        parser.add_argument('-nbins', help='The number of nbins used for smudgeplot matrix (nbins x nbins) [default autodetection]', type=int, default=0)
        # parser.add_argument('-k', help='The length of the kmer.', default=21)
        parser.add_argument('--homozygous', action="store_true", default = False, help="Assume no heterozygosity in the genome - plotting a paralog structure; [default False]")
        self.__assignArguments__(parser, sys.argv)

    def __assignArguments__(self, parser, argv):
        parser.parse_args(argv[2:]) # execution of the local parser
        self.arguments = self.parser.parse_args(argv[1:]) # This is done here to separate "-h" for each parser

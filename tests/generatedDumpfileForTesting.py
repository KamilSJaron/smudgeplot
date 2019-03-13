'''
- load 5000 nt seuqence, where 2500 are unique, 500x2 are duplicated and 500x3 are triplicated
- alternate the sequence as the "second haplotype"
'''

from random import choice
from numpy import linspace
from numpy.random import normal

h1_genome_file = open('tests/data/fake_genome_h1.fa', 'r')

h1_genome_file.readline()
h1_sequence = h1_genome_file.readline().rstrip('\n')

h1_genome_file.close()

k = 21
kmer_database = dict()

for pos in (range(len(h1_sequence) - k + 1)):
    kmer = h1_sequence[pos:pos + k]
    kmer_database[kmer] = kmer_database.get(kmer, 0) + 1

kmer_feqequency_hist = dict()
for value in kmer_database.values():
    kmer_feqequency_hist[value] = kmer_feqequency_hist.get(value, 0) + 1

sum(kmer_feqequency_hist.values())
# 3520 kmers in the dataset now
# 2540 are from unique seq
# 500 are there 2x
# 400 are there 3x

h2_sequence_list = list(h1_sequence)
# set the important one

###### 1 heterozygous ######
# AAAAAAAAAATAAAAAAAAAA
# AAAAAAAAAAGAAAAAAAAAA
# AB
h2_sequence_list[13] = 'G'

# 500 ... 2500
def pick_alternative_nt(nt):
    nts = set(['A', 'C', 'T', 'G'])
    nts.remove(nt)
    return(choice(list(nts)))

###### 50 heterozygous ######
# AB
for position in linspace(530, 1755, 50):
    position = int(position)
    nt = h2_sequence_list[position]
    h2_sequence_list[position] = pick_alternative_nt(nt)

###### 20 linked heterozygous ######
# 10 pairs 10 bases from each other
# AB, but only if all kmers are considered
for position in linspace(2030, 2390, 10):
    position = int(position)
    h2_sequence_list[position] = pick_alternative_nt(h2_sequence_list[position])
    h2_sequence_list[position + 10] = pick_alternative_nt(h2_sequence_list[position + 10])

# 2500 ... 3500 het kmers in duplicated region
# 10 heterozygous in both loci; 10 heteryzygous in one
# AABB type
for position in linspace(2530, 2737, 10):
    position = int(position)
    nt = pick_alternative_nt(h2_sequence_list[position])
    h2_sequence_list[position] = nt
    h2_sequence_list[position + 500] = nt

# AAAB type
for position in linspace(2760, 2967, 10):
    position = int(position)
    nt = pick_alternative_nt(h2_sequence_list[position])
    h2_sequence_list[position] = nt

# 3500 ... 5000 het kmers in triplicated region
# AAAABB type
for position in linspace(3530, 3737, 10):
    position = int(position)
    nt = pick_alternative_nt(h2_sequence_list[position])
    h2_sequence_list[position] = nt
    h2_sequence_list[position + 500] = nt

h2_sequence = "".join(h2_sequence_list)

# add the second haplytype to database
for pos in (range(len(h2_sequence) - k + 1)):
    kmer = h2_sequence[pos:pos + k]
    kmer_database[kmer] = kmer_database.get(kmer, 0) + 1

kmer_feqequency_hist = dict()
for value in kmer_database.values():
    kmer_feqequency_hist[value] = kmer_feqequency_hist.get(value, 0) + 1

with open('tests/data/toy_kmer_k21.dump', 'w') as dump_file:
    sorted_kmers = list(kmer_database.keys())
    sorted_kmers.sort()
    for kmer in sorted_kmers:
        freq = int(normal(kmer_database[kmer] * 100, 10))
        dump_file.write(kmer + '\t' + str(freq) + '\n')



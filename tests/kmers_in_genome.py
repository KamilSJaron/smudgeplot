from collections import defaultdict

kmers_file_name = 'tests/data/toy_middle_kmer_pairs.tsv'
refernce_file_name = 'tests/data/fake_genome_h1.fa'

h1_genome_file = open(refernce_file_name, 'r')

h1_genome_file.readline()
h1_sequence = h1_genome_file.readline().rstrip('\n')

h1_genome_file.close()

k = 21
k_for_db = int((k - 1) / 2)
kmer_database = defaultdict(list)

# make a hashmat kmer -> position
for pos in range(len(h1_sequence) - k + 1):
    kmer = h1_sequence[pos:pos + k_for_db]
    kmer_database[kmer].append(pos)

Ls = []
Rs = []
with open(kmers_file_name) as kmer_file:
    for line in kmer_file:
        L, R = line.rstrip('\n').split('?')
        Ls.append(L)
        Rs.append(R)

heterozygous_positions = []
# assuming that all the kmers are in the original sequence (and they got to be this time)
for i,L in enumerate(Ls):
    het_positions = []
    Rpositions = kmer_database[Rs[i]]
    for Lpos in kmer_database[L]:
        het_pos = Lpos + k_for_db
        if het_pos + 1 in Rpositions:
            het_positions.append(het_pos)
    heterozygous_positions.append(het_positions)

# lambda l: [item for sublist in l for item in sublist]
len(heterozygous_positions)
heterozygous_positions_flat = [item[0] for item in heterozygous_positions]
heterozygous_positions_flat.sort()

def checkPositions(pos_from, pos_to):
    return(str(sum([pos > pos_from and pos < pos_to for pos in heterozygous_positions_flat])))


print("Expected 1, found: " + checkPositions(1,500))
print("Expected 50, found: " + checkPositions(500,2000))
print("Expected 0, found: " + checkPositions(2000,2500))
print("Expected 20, found: " + checkPositions(2500,3500))
print("Expected 10, found: " + checkPositions(3500,5000))

# heterozygous_positions_flat[51:70]
# # 2668 --> positive control
# h1_sequence[2658:2668] in Ls
# h1_sequence[2669:2678] in Rs
# # 2691 --> one missing!
# L = h1_sequence[2681:2691]; L in Ls
# R = h1_sequence[2692:2702]; R in Rs
#
# kmer_database[L]
# kmer_database[R]
# both 2x in the genome
#
# It's
# GTGCTCACGGAACTTACTGTA
# GTGCTCACGGGACTTACTGTA
# and there are no other kmers starting the first 10nt; therefore it SHOULD be included
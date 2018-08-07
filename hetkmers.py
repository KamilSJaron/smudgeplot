from collections import defaultdict, Counter
from itertools import combinations
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pickle

def kmer_to_int(kmer):
  """This function turns a kmer string into an integer."""
  k = len(kmer)
  letter_to_int = {'A':1, 'C':2, 'G':3, 'T':4}
  return sum([4**(k-i-1) * letter_to_int[letter] for i, letter in enumerate(kmer)])

def get_1away_pairs(local_index_to_kmer, k):
  """local_index_to_kmer is a dictionary where the value is a kmer portion, and the key is the index of the original kmer in which the kmer portion is found. get_1away_pairs returns a list of pairs of indices where each pair of indices correspond to a pair of kmers different in exactly one base."""
  
  #This is the base case for the recursion. Return every pair of indices where the kmers at those indices differ at exactly one base.
  if k == 1:
    return [(i,j) for (i,j) in combinations(local_index_to_kmer, 2) if local_index_to_kmer[i][0] != local_index_to_kmer[j][0]]
  
  #Get the two halves of the kmer
  k_L = k//2
  k_R = k-k_L
  
  #initialize dictionaries in which the key is the hash of half of the kmer, and the value is a list of indices of the kmers with that same hash
  kmer_L_hashes = defaultdict(list)
  kmer_R_hashes = defaultdict(list)
  
  #initialize pairs, which will be returned by get_1away_pairs
  pairs = []
  
  #for each kmer, calculate its left hash and right hash, then add its index to the corresponding entries of the dictionary
  for i, kmer in local_index_to_kmer.items():
    kmer_L = kmer[:k_L]
    kmer_R = kmer[k_L:]
    kmer_L_hashes[kmer_to_int(kmer_L)] += [i]
    kmer_R_hashes[kmer_to_int(kmer_R)] += [i]
  
  #for each left hash in which there are multiple kmers with that left hash, find the list of pairs in which the right half differs by 1. (aka, if left half matches, recurse on right half).
  for kmer_L_hash_indices in kmer_L_hashes.values(): #same in first half
    if len(kmer_L_hash_indices) > 1:
      pairs += get_1away_pairs({kmer_L_hash_index:local_index_to_kmer[kmer_L_hash_index][k_L:] for kmer_L_hash_index in kmer_L_hash_indices}, k_R) #differ by 1 in right half
  
  #for each right hash in which there are multiple kmers with that right hash, find the list of pairs in which the left half differs by 1. (aka, if right half matches, recurse on left half).
  for kmer_R_hash_indices in kmer_R_hashes.values(): #same in right half
    if len(kmer_R_hash_indices) > 1:
      pairs += get_1away_pairs({kmer_R_hash_index:local_index_to_kmer[kmer_R_hash_index][:k_L] for kmer_R_hash_index in kmer_R_hash_indices}, k_L) #differ by 1 in left half
  return(pairs)

def get_2away_pairs(local_index_to_kmer, k):
  """local_index_to_kmer is a dictionary where the value is a kmer portion, and the key is the index of the original kmer in which the kmer portion is found. get_2away_pairs returns a list of pairs where each pair of indices corresponds to a pair of kmers different in exactly two bases."""
  
  #These are the base cases for the recursion. If k==1, the kmers obviously can't differ in exactly two bases, so return an empty list. if k==2, return every pair of indices where the kmers at those indices differ at exactly two bases.
  if k == 1:
    return []
  if k == 2:
    return [(i, j) for (i,j) in combinations(local_index_to_kmer, 2) if local_index_to_kmer[i][0] != local_index_to_kmer[j][0] and local_index_to_kmer[i][1] != local_index_to_kmer[j][1]]
  
  #Get the two halves of the kmer
  k_L = k//2
  k_R = k-k_L
  
  #initialize dictionaries in which the key is the hash of half of the kmer, and the value is a list of indices of the kmers with that same hash
  kmer_L_hashes = defaultdict(list)
  kmer_R_hashes = defaultdict(list)
  
  #initialize pairs, which will be returned by get_1away_pairs
  pairs = []
  
  #initialize dictionaries containing the left halves and the right halves (since we will have to check cases where the left half differs by 1 and the right half differs by 1)
  local_index_to_kmer_L = {}
  local_index_to_kmer_R = {}
  
  #for each kmer, calculate its left hash and right hash, then add its index to the corresponding entries of the dictionary
  for i, kmer in local_index_to_kmer.items():
    kmer_L = kmer[:k_L]
    kmer_R = kmer[k_L:]
    local_index_to_kmer_L[i] = kmer_L
    local_index_to_kmer_R[i] = kmer_R
    kmer_L_hashes[kmer_to_int(kmer_L)] += [i]
    kmer_R_hashes[kmer_to_int(kmer_R)] += [i]
  
  #for each left hash in which there are multiple kmers with that left hash, find the list of pairs in which the right half differs by 2. (aka, if left half matches, recurse on right half).
  for kmer_L_hash_indices in kmer_L_hashes.values(): #same in first half
    if len(kmer_L_hash_indices) > 1:
      pairs += get_2away_pairs({kmer_L_hash_index:local_index_to_kmer[kmer_L_hash_index][k_L:] for kmer_L_hash_index in kmer_L_hash_indices}, k_R) #differ by 2 in right half
  
  #for each right hash in which there are multiple kmers with that right hash, find the list of pairs in which the left half differs by 2. (aka, if right half matches, recurse on left half).
  for kmer_R_hash_indices in kmer_R_hashes.values(): #same in second half
    if len(kmer_R_hash_indices) > 1:
      pairs += get_2away_pairs({kmer_R_hash_index:local_index_to_kmer[kmer_R_hash_index][:k_L] for kmer_R_hash_index in kmer_R_hash_indices}, k_L) #differ by 2 in left half
  
  #Find matching pairs where the left half is one away, and the right half is one away
  possible_pairs_L = set(get_1away_pairs(local_index_to_kmer_L,k_L))
  possible_pairs_R = set(get_1away_pairs(local_index_to_kmer_R,k_R))
  pairs += list(possible_pairs_L.intersection(possible_pairs_R))
  return(pairs)


###This code has not been cleaned... needs to be edited!!!
def get_3away_pairs(kmers):
  """kmers is a list of kmers. get_3away_pairs returns a list of pairs where each pair of kmers is different in exactly three bases."""
  k = len(kmers[0])
  if k == 1 or k==2:
    return []
  if k == 3:
    return [pair for pair in combinations(kmers, 2) if pair[0][0] != pair[1][0] and pair[0][1] != pair[1][1] and pair[0][2] != pair[1][2]]
  k_L = k//2
  k_R = k-k_L
  kmer_L_hashes = defaultdict(list)
  kmer_R_hashes = defaultdict(list)
  pairs = []
  kmers_L = []
  kmers_R = []
  for i, kmer in enumerate(kmers):
    kmer_L = kmer[:k_L]
    kmer_R = kmer[k_L:]
    #print(kmer_L)
    #print(kmer_R)
    kmers_L.append(kmer_L)
    kmers_R.append(kmer_R)
    kmer_L_hashes[kmer_to_int(kmer_L)] += [i]
    kmer_R_hashes[kmer_to_int(kmer_R)] += [i]
  for kmer_L_hash in kmer_L_hashes.values(): #same in first half
    if len(kmer_L_hash) > 1:
      kmer_L = kmers[kmer_L_hash[0]][:k_L] #first half
      pairs += [tuple(kmer_L + kmer for kmer in pair) for pair in get_3away_pairs([kmers[i][k_L:] for i in kmer_L_hash])] #differ by 3 in second half
  for kmer_R_hash in kmer_R_hashes.values(): #same in second half
    if len(kmer_R_hash) > 1:
      kmer_R = kmers[kmer_R_hash[0]][k_L:] #second half
      #print(kmer_R)
      pairs += [tuple(kmer + kmer_R for kmer in pair) for pair in get_3away_pairs([kmers[i][:k_L] for i in kmer_R_hash])] #differ by 3 in first half
  possible_pairs = []
  possible_pairs_L = get_1away_pairs(kmers_L)
  possible_pairs_R = get_2away_pairs(kmers_R)
  #print(kmers_L)
  #print(kmers_R)
  #print(possible_pairs_L)
  #print(possible_pairs_R)
  for possible_pair_L in possible_pairs_L:
    for possible_pair_R in possible_pairs_R:
      possible_kmer1 = possible_pair_L[0]+possible_pair_R[0]
      possible_kmer2 = possible_pair_L[1]+possible_pair_R[1]
      if possible_kmer1 in kmers and possible_kmer2 in kmers:
        pairs += [(possible_kmer1, possible_kmer2)]
  possible_pairs = []
  possible_pairs_L = get_2away_pairs(kmers_L)
  possible_pairs_R = get_1away_pairs(kmers_R)
  for possible_pair_L in possible_pairs_L:
    for possible_pair_R in possible_pairs_R:
      possible_kmer1 = possible_pair_L[0]+possible_pair_R[0]
      possible_kmer2 = possible_pair_L[1]+possible_pair_R[1]
      if possible_kmer1 in kmers and possible_kmer2 in kmers:
        pairs += [(possible_kmer1, possible_kmer2)]
  return(pairs)

#Load the particular dumps file you wish
#These were created using jellyfish dump -c -L lower -U upper SRR_k21.jf > SRR_k21.dumps
dumps_file = 'ERR2135445.dumps' #aric1 -L 20 -U 350
dumps_file = 'SRR801084_k21.dumps' #avag1 -L 30 -U 300
dumps_file = 'SRR4242457_k21.dumps' #mare2 -L 13 -U 132
dumps_file = 'SRR4242472_k21.dumps' #ment1 -L 50 -U 350
dumps_file = 'SRR4242474_k21.dumps' #mflo2 -L 60 -U 400
dumps_file = 'SRR4242467_k21.dumps' #minc3 -L 25 -U 210
dumps_file = 'SRR4242462_k21.dumps' #mjav2 -L 80 -U 600
dumps_file = 'ERR2135453_k21.dumps' #rmac1 -L 100 -U 700
dumps_file = 'ERR2135451_k21.dumps' #rmag1 -L 60 -U 500

#create a pandas dataframe of the dumps file
df = pd.read_csv(dumps_file, sep = ' ', header = None)

#get a list of kmers and coverages
kmers = list(df[0])
coverages = list(df[1])

#delete the dataframe
del df

#get the value of k
k = len(kmers[0])

#turn the kmers and coverages into a dictionary, where the key is the index of the original kmer
test_kmers = {i:kmer for (i, kmer) in enumerate(kmers)}
test_coverages = {i:coverage for (i, coverage) in enumerate(coverages)}

#delete the kmers and coverages lists
del kmers
del coverages

#get the one away pairs and save it to a pickle file
one_away_pairs = get_1away_pairs(test_kmers, k)
pickle.dump(one_away_pairs, open('one_away_pairs.p', 'wb'))

#get the two away pairs and save it to a pickle file
#two_away_pairs = get_2away_pairs(test_kmers, k)
#pickle.dump(two_away_pairs, open('two_away_pairs.p', 'wb'))

#To determine the families, use the "connected_components" function of a networkx graph
#add each of the pairs as edges to the graph, and save it to a pickle file
G = nx.Graph()
for one_away_pair in one_away_pairs:
  G.add_edge(*one_away_pair)
#for two_away_pair in two_away_pairs:
#  G.add_edge(*two_away_pair)
pickle.dump(G, open('G.p', 'wb'))

#get the lengths of the families, save it to a pickle file, and then look at the Counter summary
component_lengths = [len(x) for x in nx.connected_components(G)]
pickle.dump(component_lengths, open('component_lengths.p', 'wb'))
Counter(component_lengths)

#get the families with 2-6 members
families_2 = []
families_3 = []
families_4 = []
families_5 = []
families_6 = []
for component in nx.connected_components(G):
  if len(component) == 2:
    families_2.append(list(component))
  elif len(component) == 3:
    families_3.append(list(component))
  elif len(component) == 4:
    families_4.append(list(component))
  elif len(component) == 5:
    families_5.append(list(component))
  elif len(component) == 6:
    families_6.append(list(component))

#save families_2 to a pickle file, then plot a histogram of the coverages
pickle.dump(families_2, open('families_2.p', 'wb'))
coverages_2 = [sum([test_coverages[i] for i in family]) for family in families_2]
pickle.dump(coverages_2, open('coverages_2.p', 'wb'))
plt.hist(coverages_2, bins = 1000)
plt.savefig('coverages_2_hist.png')
plt.close()

#save families_3 to a pickle file, then plot a histogram of the coverages
pickle.dump(families_3, open('families_3.p', 'wb'))
coverages_3 = [sum([test_coverages[i] for i in family]) for family in families_3]
plt.hist(coverages_3, bins = 1000)
plt.savefig('coverages_3_hist.png')
plt.close()

#n, bins, patches = plt.hist(coverages_3, bins = 1000)
#bins[np.argmax(n)]

#save families_4 to a pickle file, then plot a histogram of the coverages
pickle.dump(families_4, open('families_4.p', 'wb'))
coverages_4 = [sum([test_coverages[i] for i in family]) for family in families_4]
plt.hist(coverages_4, bins = 1000)
plt.savefig('coverages_4_hist.png')
plt.close()

#save families_5 to a pickle file, then plot a histogram of the coverages
pickle.dump(families_5, open('families_5.p', 'wb'))
coverages_5 = [sum([test_coverages[i] for i in family]) for family in families_5]
plt.hist(coverages_5, bins = 1000)
plt.savefig('coverages_5_hist.png')
plt.close()

#save families_6 to a pickle file, then plot a histogram of the coverages
pickle.dump(families_6, open('families_6.p', 'wb'))
coverages_6 = [sum([test_coverages[i] for i in family]) for family in families_6]
plt.hist(coverages_6, bins = 1000)
plt.savefig('coverages_6_hist.png')
plt.close()

###some code to load previously saved pickle files
test_kmers = pickle.load(open('test_kmers.p', 'rb'))
test_coverages = pickle.load(open('test_coverages.p', 'rb'))
G = pickle.load(open('G.p', 'rb'))
component_lengths = pickle.load(open('component_lengths.p', 'rb'))
families_2 = pickle.load(open('families_2.p', 'rb'))
coverages_2 = pickle.load(open('coverages_2.p', 'rb'))


###Everything below this is just scratch work
#f = open('ERR2135445_l20_u100.fa', 'r')
#g = open('new.fa', 'w')
#for line in f:
#  if line.startswith('>'):
#    g.write('>' + str(int(line[1:-1])+10000) + '\n')
#  else:
#    g.write(line)
#f.close()
#g.close()

#get_3away_pairs(['AAAAAAAA', 'AACTAAGA', 'AACAATGA', 'AAAAATCG'])


#get_1away_pairs(['AAA', 'AAC'])

#kmers = [''.join([random.choice('ACGT') for _ in range(20)]) for _ in range(10)]

#df2 = df[:1000000]

#for pair in pairs:
#  #f.write(str(df2[df2[0] == pair[0]].iloc[0,1])+'\n')
#  #f.write(str(df2[df2[0] == pair[1]].iloc[0,1])+'\n')
#  [x[1] for x in pairs if x[0] == pair[0]]+[x[0] for x in pairs if x[1] == pair[0]]
#  a = df2[df2[0] == pair[0]].iloc[0,1]/89.2
#  b = df2[df2[0] == pair[1]].iloc[0,1]/89.2
#  f.write(str((a, b, a+b))+'\n')

#Counter([min([Counter(pair[0])[x] for x in ['A', 'C', 'G', 'T']]) for pair in pairs])

#high_complexity_pairs = [pair for pair in pairs if min([Counter(pair[0])[x] for x in ['A', 'C', 'G', 'T']])==5]

#for hcpair in high_complexity_pairs:
#  f.write(str(df2[df2[0] == hcpair[0]]))
#  f.write(str(df2[df2[0] == hcpair[1]]))

#570620  TAAAATAATTTTTTTCTTAAA  115
#878881  TAAAATAATTTTTTTCTAAAA  67
#182

#526664  AATTACCATTCAACCAGTTTC  156
#922303  AATTACCATTCAACCAGATTC  166
#322

#394517  AAGAGAAAAGAAAAAAGTAAT  180
#788086  AAGAGAAAAGAAAAAAGAAAT  180
#360

#420665  AAAAAAAAGTGTTTTACTTTG  119
#946878  AAAAAAAAGTGTTTTACTCTG  95
#214

#594426  ACAAAATATTACCTTTATCTA  117
#768315  ACAAAATATTACCTTTATTTA  152
#269

#536269  ACAGATTGGCTTGTTTGAGCC  103
#711261  ACAGATTGGCTTGTTTGAACC  99

#383862  ATTTCATTTGTTAGAAAAAAA  139
#907248  ATTTCATTTGTTAGAAAAGAA  162

#438051  TCAACAGAAAATAATGGAGCA  152
#962365  TCAACAGAAAATAATGGAACA  143

#425231  AAAAAAAAACGAAAAAATTTT  15
#734086  AAAAAAAAACGAAAAAAATTT  18

#607197  AAAAAAAAACACGACATGTTT  154
#783001  AAAAAAAAACACGACATGCTT  134


#test_kmers = {i:kmer for (i, kmer) in enumerate(kmers[:100000])}

#members = [x[0] for x in one_away_pairs] + [x[1] for x in one_away_pairs] + [x[0] for x in two_away_pairs] + [x[1] for x in two_away_pairs]
#G = nx.Graph()
#for one_away_pair in one_away_pairs:
#  G.add_edge(*one_away_pair)
#for two_away_pair in two_away_pairs:
#  G.add_edge(*two_away_pair)

#component_lengths = [len(x) for x in nx.connected_components(G)]
#Counter(component_lengths)
#families = [list(x) for x in nx.connected_components(G) if len(x) == 2]
#coverages = [df2.iloc[families[i][0], 1]+df2.iloc[families[i][1], 1] for i in range(len(families))]
#plt.hist(coverages, bins = 100)
#plt.savefig('coverages_hist.png')
#plt.close()

#families_3 = [list(x) for x in nx.connected_components(G) if len(x) == 3]
#coverages_3 = [df2.iloc[families_3[i][0], 1]+df2.iloc[families_3[i][1], 1] for i in range(len(families_3))]
#plt.hist(coverages_3, bins = 100)
#plt.savefig('coverages_3_hist.png')
#plt.close()

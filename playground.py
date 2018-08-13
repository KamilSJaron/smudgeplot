#-----
# WHat I tried ot make plots work
# https://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#-------

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

#  dumps_file = 'kmers_dump_L120_U1500.tsv'





plt.hist(coverages_2, bins = 1000)
plt.savefig('coverages_2_hist.png')
plt.close()

#  then plot a histogram of the coverages

plt.hist(coverages_3, bins = 1000)
plt.savefig('coverages_3_hist.png')
plt.close()

#n, bins, patches = plt.hist(coverages_3, bins = 1000)
#bins[np.argmax(n)]

#save families_4 to a pickle file, then plot a histogram of the coverages

plt.hist(coverages_4, bins = 1000)
plt.savefig('coverages_4_hist.png')
plt.close()


plt.hist(coverages_5, bins = 1000)
plt.savefig('coverages_5_hist.png')
plt.close()

#save families_6 to a pickle file, then plot a histogram of the coverages
plt.hist(coverages_6, bins = 1000)
plt.savefig('coverages_6_hist.png')
plt.close()

###some code to load previously saved pickle files
# test_kmers = pickle.load(open('test_kmers.p', 'rb'))
# test_coverages = pickle.load(open('test_coverages.p', 'rb'))
G = pickle.load(open('G.p', 'rb'))
component_lengths = pickle.load(open('component_lengths.p', 'rb'))
families_2 = pickle.load(open('families_2.p', 'rb'))
coverages_2 = pickle.load(open('coverages_2.p', 'rb'))
# one_away_pair = pickle.load(open('one_away_pairs.p', 'rb'))

# perhaps faster way how to calculate coverages_2
# coverages_2 = [test_coverages[cov_i1] + test_coverages[cov_i2] for cov_i1, cov_i2 in families_2]


#-----
#   for coverage in coverages_2:
#

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

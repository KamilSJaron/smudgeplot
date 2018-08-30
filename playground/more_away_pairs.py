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
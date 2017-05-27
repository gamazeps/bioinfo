import random
from collections import defaultdict

def nuc_to_num(n):
    nucleotides = ['A', 'C', 'G' , 'T']
    return nucleotides.index(n)

def laplace_profile_from_motif(motif, k):
    t = len(motif)
    profile = [list(), list(), list(), list()]
    for i in range(0, k):
        scores = [1.0, 1.0, 1.0, 1.0]
        for j in range(0, t):
            scores[nuc_to_num(motif[j][i])] += 1
        for i in range(0, 4):
            profile[i].append(scores[i] / (t+4))
    return profile

def profile_from_motif(motif, k):
    t = len(motif)
    profile = [list(), list(), list(), list()]
    for i in range(0, k):
        scores = [0.0, 0.0, 0.0, 0.0]
        for j in range(0, t):
            scores[nuc_to_num(motif[j][i])] += 1
        for i in range(0, 4):
            profile[i].append(scores[i] / t)
    return profile

def score(motifs, k):
    result = 0
    for i in range(0, k):
        local = defaultdict(int)
        for motif in motifs:
            local[motif[i]] += 1
        m = 0
        local_consensus = ''
        for c in local:
            if local[c] > m:
                m = local[c]
                local_consensus = c
        for (key, v) in local.iteritems():
            if key != local_consensus:
                result += v
    return result

def profile_most_probable(text, k, profile):
    prob = -1.0
    res = ''
    for i in range(0, len(text) - k + 1):
        p = prob_from_profile(text[i: i+k], profile)
        if p > prob:
            prob = p
            res = text[i: i+k]
    return res

def prob_from_profile(pattern, profile):
    k = len(pattern)
    acc = 1
    for i in range(0, k): 
        n = nuc_to_num(pattern[i])
        acc *= profile[n][i]
    return acc

def sample_rand_kmer(text, k):
    i = random.randint(0, len(text) - k)
    return text[i: i+k]

def rand_kmer_with_proba(text, k, profile):
    probas = [prob_from_profile(text[i: i+k], profile)
                  for i in range(0, len(text) - k + 1)]
    total = sum(probas)
    r = random.random() * total
    acc = 0.0 
    for i in range(0, len(text) - k + 1):
        acc += probas[i]
        if acc >= r:
            return text[i: i+k]

def gibbs_sampler(dnas, k, t, N):
    motifs = [sample_rand_kmer(sample, k) for sample in dnas]
    best_motifs = motifs
    for _ in range(0, N):
        i = random.randint(0, t - 1)
        profile = laplace_profile_from_motif(motifs[0: i] + motifs[i+1: t], k)
        motifs[i] = rand_kmer_with_proba(dnas[i], k, profile)
        if score(motifs, k) < score(best_motifs, k):
            best_motifs = motifs
    return best_motifs

def randomized_motif_search(dnas, k, t):
    motifs = [sample_rand_kmer(sample, k) for sample in dnas]
    best_motifs = motifs
    while True:
        profile = laplace_profile_from_motif(motifs, k)
        motifs = [profile_most_probable(dna, k, profile) for dna in dnas]
        if score(motifs, k) < score(best_motifs, k):
            best_motifs = motifs
        else:
            return best_motifs

def iterate_rand_motif(dnas, k, t, N):
    n = 30
    m = 4 ** k
    best_motif = None
    for i in range(0, n):
        motif = gibbs_sampler(dnas, k, t, N)
        if score(motif, k) < m:
            m = score(motif, k)
            #print(m, i)
            best_motifs = motif
    return best_motifs


##############################################################################

def kmer_from_number(n, k):
    nucleotides = ['A', 'C', 'G' , 'T']
    if k <= 1:
        return nucleotides[n]
    return kmer_from_number(n / 4, k - 1) + nucleotides[n%4]


def pprint_list(l):
    res = ''
    for e in l:
        res += str(e) + ' '
    print(res[:-1])

def print_set(s):
    for e in s:
        print(e)

###############################################################################

def parse_file_1():
    workfile = 'test_data.txt'
    workfile = 'dataset_163_4.txt'
    f = open(workfile, 'r')
    lines = f.readlines()
    ints = lines[0].strip().split(' ')
    k = int(ints[0]) 
    t = int(ints[1]) 
    N = int(ints[2]) 

    dnas = [s.strip() for s in lines[1:]]
    return (dnas, k, t, N)

################################################################################

(dnas, k, t, N) = parse_file_1()

print_set(iterate_rand_motif(dnas, k, t, N))

#print(score([
#    "TCTCGGGG",
#    "CCAAGGTG",
#    "TACAGGCG",
#    "TTCAGGTG",
#    "TCCACGTG"
#    ], 8))

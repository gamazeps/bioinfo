from collections import defaultdict

nucleotides = {'A', 'C', 'G' , 'T'}

def hamming_dist(a, b):
    count = 0
    for i in range(0, len(a)):
        if a[i] != b[i]:
            count += 1
    return count

def generate_neighbours(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return nucleotides

    neighbours = set()
    suffix_neighbours = generate_neighbours(pattern[1:], d)

    for s in suffix_neighbours:
        if hamming_dist(s, pattern[1:]) < d:
            for n in nucleotides:
                neighbours.add(n + s)
        else:
            neighbours.add(pattern[0] + s)

    return neighbours

def approximate_pattern_match(pattern, gene, dist):
    res = list()
    for i in range(0, len(gene) - len(pattern) + 1):
        if hamming_dist(gene[i:i+len(pattern)], pattern) <= dist:
            res.append(i)
    return res

def is_in_all_with_mismatch(sample, strings, d):
    for s in strings:
        if (len(approximate_pattern_match(sample, s, d)) < 1):
            return False
    return True


def motif_enumeration(dnas, k, d):
    res = set()

    for dna in dnas:
        for i in range(0, len(dna) - k + 1):
            neighbours = generate_neighbours(dna[i: i+k], d)
            for candidate in neighbours:
                if is_in_all_with_mismatch(candidate, dnas, d):
                    res.add(candidate)

    return res


def median_string(patterns, k):
    dist = 4 ** k
    median = set()
    for i in range(0, 4 ** k):
        pattern = kmer_from_number(i, k)
        if multiple_dists(patterns, pattern) == dist:
            median.add(pattern)
        if multiple_dists(patterns, pattern) < dist:
            dist = multiple_dists(patterns, pattern)
            median = set()
            median.add(pattern)
    return median

def multiple_dists(patterns, pattern):
    acc = 0
    for p in patterns:
        acc += min_hamming_dist(p, pattern)
    return acc

def min_hamming_dist(string, pattern):
    m = len(string)
    k = len(pattern)
    for i in range(0, len(string) - k + 1):
        if hamming_dist(pattern, string[i: i+k]) < m:
            m = hamming_dist(pattern, string[i: i+k])
    return m

def most_probable(text, k, profile):
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

def greedy_motif_search(patterns, k, t):
    best_motifs = [p[0: k] for p in patterns]
    n = len(patterns[0])

    for i in range(0, n - k + 1):
        motifs = [patterns[0][i: i+k]]
        for j in range(1, t):
            profile = laplace_profile_from_motif(motifs, k)
            motifs.append(most_probable(patterns[j], k, profile))
        if score(motifs, k) < score(best_motifs, k):
            #print(i)
            best_motifs = motifs
    return best_motifs

def score(motifs, k):
    result = 0
    for i in range(0, k):
        # establish the most common nucleotide
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
    workfile = 'dataset_156_8.txt'
    workfile = 'dataset_160_9.txt'
    #workfile = 'test_data.txt'
    f = open(workfile, 'r')
    lines = f.readlines()
    ints = lines[0].strip().split(' ')
    k = int(ints[0]) 
    t = int(ints[1]) 

    dnas = [s.strip() for s in lines[1:]]
    return (dnas, k, t)

def parse_file_2():
    workfile = 'dataset_158_9.txt'
    f = open(workfile, 'r')
    lines = f.readlines()
    k = int(lines[0].strip())

    dnas = [s.strip() for s in lines[1:]]
    return (dnas, k)

def parse_file_3():
    workfile = 'dataset_159_3.txt'
    f = open(workfile, 'r')
    lines = f.readlines()
    text = lines[0].strip()
    k = int(lines[1].strip())

    profile = list()
    for s in lines[2:]:
        tmp = list()
        for p in s.strip().split(' '):
            tmp.append(float(p))
        profile.append(tmp)

    return (text, k, profile)

##############################################################################

(k, t, dnas) = parse_file_1()

#print_set(greedy_motif_search(k, t, dnas))


samples =[
        "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
        "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
        "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"
        ]

print(median_string(samples, 7))

#pprint_set(motif_enumeration(samples, k, d))

#print(kmer_from_number(0, 2))
#print(kmer_from_number(1, 2))
#print(kmer_from_number(4, 2))
#print(kmer_from_number(5, 2))


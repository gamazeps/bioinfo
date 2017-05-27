import collections

def hamming_dist(a, b):
    count = 0
    for i in range(0, len(a)):
        if a[i] != b[i]:
            count += 1
    return count

def approximate_pattern_match(pattern, gene, dist):
    res = list()
    for i in range(0, len(gene) - len(pattern) + 1):
        if hamming_dist(gene[i:i+len(pattern)], pattern) <= dist:
            res.append(i)
    return res

def pprint_list(l):
    res = ''
    for e in l:
        res += str(e) + ' '
    print(res[:-1])

def pprint_set(s):
    for e in s:
        print(e)

nucleotides = {'A', 'C', 'G' , 'T'}

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

def my_generate(pattern, d):
    print(d)
    if d <= 0:
        return {pattern}

    neighbours = set()
    for i in range(0, len(pattern)):
        for n in nucleotides:
            tmp = pattern[:i] + n + pattern[i+1:]
            print(tmp)
            neighbours.add(tmp)
            neighbours |= my_generate(tmp, d - 1)

    return neighbours

def frequent_words_with_mismatch(gene, k, d):
    data = collections.defaultdict(int)

    for i in range(0, len(gene) - k + 1):
        ori = gene[i: i + k]
        kmers = generate_neighbours(ori, d)
        for kmer in kmers:
            data[kmer] += 1

    res = list()
    m = 0

    for k in data:
        if data[k] == m:
            res.append(k)
        elif data[k] > m:
            m = data[k]
            res = [k]

    return res


def reverse_cmp(s):
    rev = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
    }
    res = ''
    for e in s[::-1]:
        res += rev[e]
    return res

def frequent_words_with_mismatch_and_rev(gene, k, d):
    data = collections.defaultdict(int)

    for i in range(0, len(gene) - k + 1):
        ori = gene[i: i + k]
        kmers = generate_neighbours(ori, d)
        for kmer in kmers:
            data[kmer] += 1

    for i in range(0, len(gene) - k + 1):
        ori = gene[i: i + k]
        kmers = generate_neighbours(reverse_cmp(ori), d)
        for kmer in kmers:
            data[kmer] += 1

    res = list()
    m = 0

    for k in data:
        if data[k] == m:
            res.append(k)
        elif data[k] > m:
            m = data[k]
            res = [k]

    return res


def skew(s):
    count = 0
    res = list()

    for i in s:
        if i == 'C':
            count += 1
        if i == 'G':
            count -= 1
        res.append(count)

    m = res[0]
    r = list()
    for i in range(0, len(res)):
        if res[i] == m:
            r.append(i + 1)
        elif res[i] > m:
            m = res[i]
            r = list()
            r.append(i + 1)
    return r



def read_stuff():
    workfile = 'dataset_9_6.txt'
    f = open(workfile, 'r')
    pattern = f.readline().strip()
    gene = f.readline().strip()
    d = int(f.readline().strip())
    return (pattern, gene, d)

def read_stuff_2():
    workfile = 'dataset_3014_3.txt'
    f = open(workfile, 'r')
    pattern = f.readline().strip()
    k = int(f.readline().strip())
    return (pattern, k)

def read_stuff_3():
    workfile = 'dataset_9_8.txt'
    f = open(workfile, 'r')
    pattern = f.readline().strip()
    l = f.readline().strip().split(' ')
    return (pattern, int(l[0]), int(l[1]))

#print(hamming_dist('GGGCCGTTGGT', 'GGACCGTTGAC'))
#pprint_list(approximate_pattern_match('ATTCTGGA',
#            'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',
#            3
#            ))
#print(len(approximate_pattern_match('GAGG',
#            'TTTAGAGCCTTCAGAGG',
#            2
#            )))
#
#(p, g, d) = read_stuff()
#print(len(approximate_pattern_match(p, g, d)))

(p, k, d) = read_stuff_3()
#pprint_set(generate_neighbours(p, d))
#print(generate_neighbours('ACG', 1) - my_generate('ACG', 1))
#pprint_list(frequent_words_with_mismatch_and_rev(p, k, d))

#pprint_list(frequent_words_with_mismatch_and_rev('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1))

#print(skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'))

print(hamming_dist('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA', 'CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG'))
print(skew('CATTCCAGTACTTCATGATGGCGTGAAGA'))
print(len(approximate_pattern_match('TGT', 'CGTGACAGTGTATGGGCATCTTT', 1)))
print(len(generate_neighbours('ACGT', 3)))

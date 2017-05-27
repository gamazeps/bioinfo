import collections

def composition(text, k):
    comp = list()
    for i in range(0, len(text) - k + 1):
        comp.append(text[i: i+k])
    return sorted(comp)

def string_spelling(kmers):
    if len(kmers) == 0:
        return ""
    res = kmers[0]
    for i in range(1, len(kmers)):
        res += kmers[i][-1]
    return res

def overlap_graph_problem(kmers_list):
    kmers_set = set(kmers_list)

    # construct the prefix dictionnary:
    with_prefix = collections.defaultdict(list)
    for kmer in kmers_set:
        with_prefix[kmer[:-1]].append(kmer)

    res = dict()
    for kmer in kmers_set:
        res[kmer] = with_prefix[kmer[1:]]

    return res

def de_bruijn_graph(kmers):
    graph = collections.defaultdict(list)

    for kmer in kmers:
        graph[kmer[:-1]].append(kmer[1:])

    return graph

###########################################

def print_ml(c):
    for e in c:
        print(e)

def print_adjency_list(adj):
    for k in adj:
        for e in adj[k]:
            print(k + " -> " + e)

def print_de_bruijn(adj):
    nodes = sorted(adj)

    for k in nodes:
        print(k + " -> " + format_list(adj[k]))

def format_list(l):
    return ",".join(l)

###########################################

def parse_file_1():
    workfile = 'dataset_199_6.txt'
    f = open(workfile, 'r')
    k = int(f.readline().strip())
    text = f.readline().strip()
    return (text, k)

def parse_file_2():
    workfile = 'dataset_200_8.txt'
    f = open(workfile, 'r')
    lines = f.readlines()
    return map(lambda x: x.strip(), lines)

###########################################

(text, k) = parse_file_1()
kmers = parse_file_2()

#print(string_spelling(kmers))

#print_adjency_list(overlap_graph_problem(kmers))

print_de_bruijn(de_bruijn_graph(kmers))

#print_ml(composition(text, k))

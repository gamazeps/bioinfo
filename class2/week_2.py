import collections

########################################################

def eulerian_cycle(graph, root = None):
    unexplored_edges = set()
    unexplored_nodes = set()

    if len(graph) is 0:
        return list()

    for k in graph:
        for elem in graph[k]:
            unexplored_edges.add((k, elem))

    def random_walk(origin):
        path = list()
        path.append(origin)
        curr = origin
        next_node = None

        while len(graph[curr]) is not 0:
            next_node = graph[curr].pop()
            if len(graph[curr]) is not 0:
                unexplored_nodes.add(curr)
            unexplored_edges.remove((curr, next_node))
            path.append(next_node)
            curr = next_node

        return path

    # Gets a random key if needed
    if root is None:
        (first_node, values) = graph.popitem()
        graph[first_node] = values
    else:
        first_node = root

    result = random_walk(first_node)

    while len(unexplored_edges) is not 0:
        new_origin = unexplored_nodes.pop()
        cycle = random_walk(new_origin)
        index = result.index(new_origin)
        result = result[:index] + cycle + result[index + 1:]

    return result

def eulerian_path(graph):
    in_degrees = collections.defaultdict(int)
    out_degrees = collections.defaultdict(int)
    degrees = dict()
    root = None

    for node in graph:
        for out in graph[node]:
            in_degrees[out] +=1

    for node in graph:
        out_degrees[node] = len(graph[node])

    nodes = set(out_degrees.keys()).union(set(in_degrees.keys()))

    for node in nodes:
        if in_degrees[node] < out_degrees[node]:
            root = node

    return eulerian_cycle(graph, root)

###############################################################

def format_path(path):
    return '->'.join(map(str, path))

def format_kmers_path(path):
    res = path[0]
    for elem in path[1:]:
        res += elem[-1:]
    return res

#################################################################

def parse_input_graph(fname):
    workfile = open(fname, 'r')
    lines = workfile.readlines()
    lines = map(lambda x: x.strip(), lines)

    graph = collections.defaultdict(set)
    for line in lines:
        pair = line.split(' -> ')
        for elem in pair[1].split(','):
            graph[int(pair[0])].add(int(elem))

    return graph

def parse_input_kmers(fname):
    workfile = open(fname, 'r')
    _ = int(workfile.readline().strip())
    kmers = workfile.readlines()
    return map(lambda x: x.strip(), kmers)

def parse_input_kdmers(fname):
    workfile = open(fname, 'r')
    [k, gap] = map(int, workfile.readline().strip().split(' '))
    kmers = map(lambda x: x.strip().split('|'), workfile.readlines())
    return (k, gap, kmers)

##############################################################

def de_bruijn_graph(kmers):
    graph = collections.defaultdict(set)

    for kmer in kmers:
        graph[kmer[:-1]].add(kmer[1:])

    return graph

def kd_de_bruijn_graph(kmers):
    graph = collections.defaultdict(set)

    for kmer in kmers:
        graph[(kmer[0][:-1], kmer[1][:-1])].add((kmer[0][1:], [kmer[1][1:]]))

    return graph

##############################################################

def generate_kmers(k):
    result = set(['0', '1'])
    for _ in range(1, k):
        tmp = set()
        for elem in result:
            tmp.add(elem + '0')
            tmp.add(elem + '1')
        result = tmp
    return result

#############################################################

def main():
    (k, gap, graph) = parse_input_kdmers('StringReconstructionFromReadPairs.txt')
    _ = kd_de_bruijn_graph(graph)

if __name__ == '__main__':
    main()

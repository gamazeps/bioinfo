import collections

def eulerian_cycle(graph):
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
            print(curr)
            next_node = graph[curr].pop()
            if len(graph[curr]) is not 0:
                unexplored_nodes.add(curr)
            unexplored_edges.remove((curr, next_node))
            path.append(next_node)
            curr = next_node

        return path

    # Gets a random key
    (first_iter_node, values) = graph.popitem()
    graph[first_iter_node] = values

    result = random_walk(first_iter_node)

    while len(unexplored_edges) is not 0:
        new_origin = unexplored_nodes.pop()
        cycle = random_walk(new_origin)

    return result

    #while len(unexplored_edges) is not 0:
    #    curr_path = random_walk(unexplored_nodes.pop())

def parse_input(fname):
    workfile = open(fname, 'r')
    lines = workfile.readlines()
    lines = map(lambda x: x.strip(), lines)

    graph = collections.defaultdict(set)
    for line in lines:
        pair = line.split(' -> ')
        for elem in pair[1].split(','):
            graph[int(pair[0])].add(int(elem))

    return graph

def main():
    input_graph = parse_input('input.txt')
    print(eulerian_cycle(input_graph))

if __name__ == '__main__':
    main()

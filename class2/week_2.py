import collections

def eulerian_cycle(graph):
    unexplored_edges = set()
    unexplored_nodes = set()

    for k in graph:
        for e in graph[k]:
            unexplored_edges.add((k, e))

    def random_walk(origin):
        path = list()
        path.append(origin)
        curr = origin
        n = None

        while n is not origin:
            curr = n
            n = graph[curr].pop()
            if len(graph[curr]) is not 0:
                unexplored_nodes.add(curr)
            unexplored_nodes.remove((curr, n))
            path.append(n)

        return path

    first_iter_node = graph.keys()[0]
    res = random_walk(first_iter_node)

    while len(unexplored_edges) is not 0:
        curr_path = random_walk(unexplored_nodes.pop())


"""Independent expanded-graph feasibility predicate and corona constructor."""
import networkx as nx


def feasible(G, labels, allow_empty=True):
    zero = [v for v in G if labels[v] == 0]
    if not zero:
        if not allow_empty:
            return False
    elif not nx.is_connected(G.subgraph(zero)):
        return False
    for v in G:
        ns = [labels[u] for u in G[v]]
        if labels[v] == 0 and not (3 in ns or ns.count(2) >= 2):
            return False
        if labels[v] == 1 and not any(x >= 2 for x in ns):
            return False
        if labels[v] > 0:
            near = set(G[v])
            for u in G[v]:
                near.update(G[u])
            near.discard(v)
            if not any(labels[u] > 0 for u in near):
                return False
    return True


def leafy(G, counts):
    H = nx.convert_node_labels_to_integers(G)
    next_id = len(H)
    for v, count in enumerate(counts):
        for _ in range(count):
            H.add_edge(v, next_id)
            next_id += 1
    return H



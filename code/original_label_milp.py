"""Independent MILP with one-commodity flow to enforce zero-set connectivity.

The model retains all four labels. Its constraints do not use the proposed
path, cycle or corona formulas. SciPy/HiGHS optimality is required and every
returned integer labelling is checked by a separate graph predicate.
"""
import networkx as nx
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix
from verify_ocstdr import feasible


def solve(G, no_one=False, ordinary=False, nonempty=False):
    G = nx.convert_node_labels_to_integers(G)
    n = len(G)
    arcs = [(u, v) for u, v in G.edges() for u, v in ((u, v), (v, u))]
    # x(v,label), chosen root r(v), external supply s(v), active flag t,
    # directed flow on original edges. Only x,r,t are integral.
    rid, sid, tid, fid = 4*n, 5*n, 6*n, 6*n+1
    dim = fid + len(arcs)
    objective = np.zeros(dim)
    lower = np.zeros(dim)
    upper = np.ones(dim)
    integral = np.zeros(dim)
    integral[:5*n] = 1
    integral[tid] = 1
    upper[sid:sid+n] = n
    upper[fid:] = n
    for v in G:
        for k in range(4):
            objective[4*v+k] = k
        if no_one:
            upper[4*v+1] = 0
    rows, lbs, ubs = [], [], []

    def add(d, lo=-np.inf, hi=np.inf):
        rows.append(d)
        lbs.append(lo)
        ubs.append(hi)

    for v in G:
        add({4*v+k: 1 for k in range(4)}, 1, 1)
        d = {4*v: -2}
        for u in G[v]:
            d[4*u+2] = 1
            d[4*u+3] = 2
        add(d, 0)
        d = {4*v+1: -1}
        for u in G[v]:
            d[4*u+2] = 1
            d[4*u+3] = 1
        add(d, 0)
        if not ordinary:
            near = set(nx.single_source_shortest_path_length(G, v, cutoff=2)) - {v}
            # sum_{near} (1-x_u0) >= 1-x_v0.
            d = {4*u: -1 for u in near}
            d[4*v] = 1
            add(d, 1-len(near))
            add({rid+v: 1, 4*v: -1}, hi=0)
            add({sid+v: 1, rid+v: -n}, hi=0)
            add({4*v: 1, tid: -1}, hi=0)
            balance = {sid+v: 1, 4*v: -1}
            for j, (a, b) in enumerate(arcs):
                if a == v:
                    balance[fid+j] = balance.get(fid+j, 0)-1
                if b == v:
                    balance[fid+j] = balance.get(fid+j, 0)+1
            add(balance, 0, 0)
    if not ordinary:
        add(dict([(rid+v, 1) for v in G] + [(tid, -1)]), 0, 0)
        for j, (a, b) in enumerate(arcs):
            add({fid+j: 1, 4*a: -n}, hi=0)
            add({fid+j: 1, 4*b: -n}, hi=0)
        if nonempty:
            lower[tid] = 1
    else:
        upper[4*n:] = 0
    mat = lil_matrix((len(rows), dim))
    for i, d in enumerate(rows):
        for j, value in d.items():
            mat[i, j] = value
    result = milp(objective, integrality=integral,
                  bounds=Bounds(lower, upper),
                  constraints=LinearConstraint(mat.tocsr(), lbs, ubs),
                  options={'mip_rel_gap': 0, 'time_limit': 45})
    assert result.success and result.mip_gap == 0, str(result)
    labels = {v: int(np.argmax(result.x[4*v:4*v+4])) for v in G}
    weight = sum(labels.values())
    assert abs(weight-result.fun) < 1e-6
    if not ordinary:
        assert feasible(G, labels, not nonempty)
    else:
        for v in G:
            ns = [labels[u] for u in G[v]]
            assert labels[v] != 0 or 3 in ns or ns.count(2) >= 2
            assert labels[v] != 1 or any(x >= 2 for x in ns)
    if no_one:
        assert 1 not in labels.values()
    return weight, list(labels.values())



"""Reproduce the exact comparisons and tree timing experiments."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import random
import statistics
import sys
import time

import networkx as nx
import numpy as np
import scipy
from core_solvers import adjacency,check_set,rooted_milp,subsets,tree_dp

OUT=Path(__file__).resolve().parents[1]/'outputs'/'results'

from verify_ocstdr import feasible,leafy
from original_label_milp import solve as original_milp

SEED=20261020


def save(name, data):
    (OUT/name).write_text(json.dumps(data,indent=2),encoding='utf-8')


def instance_record(edges,k):
    return dict(m=len(k),edges=[list(e) for e in edges],k=list(k))


def tree_checks():
    total=0
    rerooted=0
    records=[]
    for m in range(1,7):
        trees=[nx.empty_graph(1)] if m==1 else list(nx.nonisomorphic_trees(m))
        start=time.perf_counter()
        for G in trees:
            edges=list(G.edges())
            adj=adjacency(m,edges)
            for k in itertools.product((1,2,3),repeat=m):
                expected=subsets(adj,k)
                actual=tree_dp(adj,k,count=True)
                assert (actual['value'],actual['count'])==(expected['value'],expected['count'])
                if m<=4:
                    for root in range(m):
                        other=tree_dp(adj,k,count=True,root=root)
                        assert (other['value'],other['count'])==(expected['value'],expected['count'])
                        rerooted+=1
                records.append(dict(**instance_record(edges,k),value=actual['value'],
                                    count=actual['count'],selected=sorted(actual['selected'])))
                total+=1
        print('tree exhaustive',m,'cumulative',total,flush=True)
    # Positive costs, huge multiplicities, all-zero costs and ties.
    for m in (1,2,3,8,12):
        for G in (nx.path_graph(m),nx.star_graph(m-1)):
            adj=adjacency(m,list(G.edges()))
            for k in ([2]*m,[3]*m,[10**9 if v%2 else 1 for v in range(m)]):
                a,b=tree_dp(adj,k,True),subsets(adj,k)
                assert (a['value'],a['count'])==(b['value'],b['count'])
    save('tree_exact_instances.json',records)
    return dict(instances=total,reroot_checks=rerooted,extra_edge_cases=30)


def direct_brute(G):
    """Independent direct four-label predicate via precomputed neighbour lists."""
    n=len(G)
    ns=[list(G[v]) for v in range(n)]
    near=[]
    for v in range(n):
        seen=set(ns[v])
        for u in ns[v]:
            seen.update(ns[u])
        seen.discard(v)
        near.append(list(seen))
    best,count=3*n+1,0
    for lab in itertools.product(range(4),repeat=n):
        weight=sum(lab)
        if weight>best:
            continue
        valid=True
        for v in range(n):
            x=lab[v]
            if x==0 and not (any(lab[u]==3 for u in ns[v]) or sum(lab[u]==2 for u in ns[v])>=2):
                valid=False;break
            if x==1 and not any(lab[u]>=2 for u in ns[v]):
                valid=False;break
            if x>0 and not any(lab[u]>0 for u in near[v]):
                valid=False;break
        if not valid:
            continue
        zero={v for v in range(n) if lab[v]==0}
        if zero:
            seen={next(iter(zero))}
            stack=list(seen)
            while stack:
                for u in ns[stack.pop()]:
                    if u in zero and u not in seen:
                        seen.add(u);stack.append(u)
            if seen!=zero:
                continue
        if weight<best:
            best,count=weight,1
        else:
            count+=1
    return best,count


def graph_checks():
    atlas=[g for g in nx.graph_atlas_g() if 1<=len(g)<=5 and nx.is_connected(g)]
    records=[]
    roots=original=direct=0
    rng=random.Random(SEED)
    for index,G in enumerate(atlas):
        m=len(G)
        edges=list(G.edges())
        adj=adjacency(m,edges)
        patterns=list(dict.fromkeys([tuple([1]*m),tuple(1+v%3 for v in range(m)),
                                    tuple(rng.randint(1,4) for _ in range(m))]))
        for k in patterns:
            expected=subsets(adj,k)
            start=time.perf_counter()
            actual=rooted_milp(adj,k)
            elapsed=time.perf_counter()-start
            assert actual['value']==expected['value']
            roots+=m
            record=dict(**instance_record(edges,k),value=actual['value'],
                        core_set_count=expected['count'],milp_seconds=elapsed,
                        rooted_statuses=actual['statuses'])
            expanded=leafy(G,k)
            witness={v:(0 if v in actual['selected'] else 2) for v in G}
            for leaf in set(expanded)-set(G):
                core=next(iter(expanded[leaf]))
                witness[leaf]=2 if core in actual['selected'] else 1
            assert sum(witness.values())==expected['value'] and feasible(expanded,witness)
            if m<=4:
                independent,_=original_milp(expanded)
                assert independent==expected['value']
                original+=1
                record['original_label_milp']=independent
            if len(expanded)<=8:
                value,number=direct_brute(expanded)
                assert value==expected['value']
                direct+=1
                record['direct_label_value']=value
                record['direct_label_optima']=number
            records.append(record)
        print('general atlas',index+1,len(atlas),flush=True)
    save('general_exact_instances.json',records)
    # Explicit counterexample to confusing core-set counts with function counts.
    G=nx.empty_graph(1)
    core=subsets([[]],[1])
    lab=direct_brute(leafy(G,[1]))
    assert core['count']==1 and lab==(3,2)
    return dict(instances=len(records),rooted_models=roots,
                original_label_milp_checks=original,direct_label_checks=direct,
                count_counterexample=dict(core_minimizers=1,function_minimizers=2))


def random_general_checks():
    records=[]
    rng=random.Random(SEED+1)
    for m in (6,8,10,12):
        for density in (.2,.5):
            for replicate in range(3):
                # Connected by construction: recursive tree plus independent missing edges.
                edges={(rng.randrange(v),v) for v in range(1,m)}
                for u in range(m):
                    for v in range(u+1,m):
                        if rng.random()<density:
                            edges.add((u,v))
                edges=sorted(edges)
                k=[rng.randint(1,4) for _ in range(m)]
                adj=adjacency(m,edges)
                start=time.perf_counter();expected=subsets(adj,k);ts=time.perf_counter()-start
                start=time.perf_counter();actual=rooted_milp(adj,k);tm=time.perf_counter()-start
                assert actual['value']==expected['value']
                records.append(dict(**instance_record(edges,k),density_parameter=density,
                                    replicate=replicate,value=actual['value'],
                                    subset_seconds=ts,milp_seconds=tm))
        print('random general',m,flush=True)
    save('random_general_instances.json',records)
    return dict(instances=len(records),rooted_models=sum(r['m'] for r in records))


def family_edges(family,m,rng):
    if family=='path':
        return [(v-1,v) for v in range(1,m)]
    if family=='star':
        return [(0,v) for v in range(1,m)]
    if family=='binary':
        return [((v-1)//2,v) for v in range(1,m)]
    return [(rng.randrange(v),v) for v in range(1,m)]


def benchmarks():
    rng=random.Random(SEED+2)
    raw=[]
    aggregates=[]
    for m in (100,1000,10000,100000):
        for family in ('path','star','binary','recursive'):
            times=[]
            for instance in range(3):
                edges=family_edges(family,m,rng)
                k=[rng.randint(1,4) for _ in range(m)]
                adj=adjacency(m,edges)
                tree_dp(adj,k)  # one untimed warm-up per instance
                signature=hashlib.sha256(json.dumps([edges,k],separators=(',',':')).encode()).hexdigest()
                for repeat in range(5):
                    start=time.perf_counter()
                    answer=tree_dp(adj,k)
                    elapsed=time.perf_counter()-start
                    assert answer['transitions']==m-1
                    times.append(elapsed)
                    raw.append(dict(m=m,family=family,instance=instance,repeat=repeat,
                                    input_sha256=signature,seconds=elapsed,value=answer['value'],
                                    selected_count=len(answer['selected']),transitions=m-1))
            aggregates.append(dict(m=m,family=family,calls=len(times),median_seconds=statistics.median(times),
                                   min_seconds=min(times),max_seconds=max(times)))
            print('benchmark',m,family,'median',statistics.median(times),flush=True)
    save('tree_timing_raw.json',raw)
    save('tree_timing_summary.json',aggregates)
    return aggregates


def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUT,help='Directory for regenerated results')
    OUT=parser.parse_args().output_dir.resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    # Save completed layers immediately so an interrupted run has usable evidence.
    metadata=dict(seed=SEED,python=sys.version,scipy=scipy.__version__,numpy=np.__version__,
                  networkx=nx.__version__,platform=platform.platform(),processor=platform.processor(),
                  logical_cpus=os.cpu_count(),timing='perf_counter wall time; includes rooting, reconstruction and result checks; excludes graph generation; no counting')
    summary=dict(metadata=metadata)
    for name,fn in [('trees',tree_checks),('atlas',graph_checks),('random_general',random_general_checks),('benchmark',benchmarks)]:
        summary[name]=fn()
        save('summary.json',summary)
    summary['wall_seconds']=time.perf_counter()-start
    summary['code_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}
    save('summary.json',summary)
    print('ALL EXACT CHECKS PASSED',summary['wall_seconds'],flush=True)


if __name__=='__main__':
    main()

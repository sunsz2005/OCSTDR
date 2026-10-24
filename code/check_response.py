"""Check the leaf-addition response independently by subset enumeration."""
import argparse
import itertools
import json
from pathlib import Path
import networkx as nx
from core_solvers import adjacency,leaf_response,subsets

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[1]/'outputs'/'results')
out=parser.parse_args().output_dir.resolve()
out.mkdir(parents=True,exist_ok=True)
checked=0
for m in range(1,6):
    graphs=[nx.empty_graph(1)] if m==1 else nx.nonisomorphic_trees(m)
    for G in graphs:
        adj=adjacency(m,list(G.edges()))
        for k in itertools.product((1,2,3),repeat=m):
            for v in range(m):
                if k[v]<2:
                    continue
                r=leaf_response(adj,k,v)
                for t in (0,1,2,5):
                    modified=list(k)
                    modified[v]+=t
                    actual=subsets(adj,modified)['value']
                    prediction=sum(k)+t+2*m+min(r['outside'],r['inside']+t)
                    assert prediction==actual
                    assert t<=actual-r['baseline']<=2*t
                    checked+=1
adj=adjacency(5,[(0,1),(1,2),(2,3),(3,4)])
k=[2,1,2,1,2]
response=leaf_response(adj,k,2)
assert (response['inside'],response['outside'],response['threshold'])==(-2,-1,1)
rows=[]
for t in range(6):
    modified=list(k);modified[2]+=t
    rows.append(dict(t=t,value=subsets(adj,modified)['value']))
result=dict(checked_predictions=checked,path_example=dict(k=k,vertex=2,response=response,values=rows))
(out/'leaf_response.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('PASS:',checked,'leaf-addition predictions; exact path threshold example.')

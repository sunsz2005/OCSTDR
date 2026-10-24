"""Build manuscript tables from executed experiment records only."""
import argparse
from pathlib import Path
import json
import statistics

folder=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--results-dir',type=Path,default=folder/'results',help='Saved or regenerated JSON records')
parser.add_argument('--output-dir',type=Path,default=folder/'outputs'/'tables')
args=parser.parse_args()
results=args.results_dir.resolve()
output=args.output_dir.resolve()
output.mkdir(parents=True,exist_ok=True)
summary=json.loads((results/'summary.json').read_text())
response=json.loads((results/'leaf_response.json').read_text())
atlas=json.loads((results/'general_exact_instances.json').read_text())
general=json.loads((results/'random_general_instances.json').read_text())
tex=[]
tex.append(r'''\begin{table}[H]
\caption{Coverage of the exact cross-checks. Input collections overlap across comparisons.}\label{tab:checks}
\centering\begin{tabular}{lrl}
\hline Check & Number & Unit\\\hline''')
rows=[('Tree DP versus subsets',summary['trees']['instances']),
      ('Alternative tree roots',summary['trees']['reroot_checks']),
      ('Additional edge-case runs',summary['trees']['extra_edge_cases']),
      ('Atlas core MILP versus subsets',summary['atlas']['instances']),
      ('Random core MILP versus subsets',summary['random_general']['instances']),
      ('Original four-label MILP',summary['atlas']['original_label_milp_checks']),
      ('Direct four-label enumeration',summary['atlas']['direct_label_checks']),
      ('Leaf-response predictions',response['checked_predictions'])]
for (label,number),unit in zip(rows,('Inputs','Runs','Runs','Inputs','Inputs','Inputs','Inputs','Predictions')):
    tex.append(label+' & '+format(number,',')+' & '+unit+r'\\')
tex.append(r'\hline\end{tabular}\end{table}')
(output/'table_checks.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')
tex=[r'''\begin{table}[H]
\caption{Tree DP wall times in milliseconds: median of fifteen calls per cell, from three instances and five timed calls each. Counting is disabled.}\label{tab:timing}
\centering\begin{tabular}{r r r r r}
\hline $m$ & Path & Star & Binary & Recursive\\\hline''']
for m in (100,1000,10000,100000):
    row=[format(m,',')]
    for family in ('path','star','binary','recursive'):
        r=next(x for x in summary['benchmark'] if x['m']==m and x['family']==family)
        row.append('%.3f'%(1000*r['median_seconds']))
    tex.append(' & '.join(row)+r'\\')
tex.append(r'\hline\end{tabular}\end{table}')
(output/'table_timing.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')
tex=[r'''\begin{table}[H]
\caption{General-core solution times in milliseconds, medians across six inputs at each order. Rooted MILP includes all root models, construction and checks; enumeration includes all subset tests.}\label{tab:general}
\centering\begin{tabular}{r r r r}
\hline $m$ & Inputs & Subsets & Core MILP\\\hline''']
for m in (6,8,10,12):
    rows=[x for x in general if x['m']==m]
    tex.append('%d & %d & %.3f & %.3f'%(m,len(rows),1000*statistics.median(x['subset_seconds'] for x in rows),1000*statistics.median(x['milp_seconds'] for x in rows))+r'\\')
tex.append(r'\hline\end{tabular}\end{table}')
(output/'table_general.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')
print('Generated three manuscript tables directly from saved results.')

# Generalized coronas: structure and exact algorithms

Reproducibility code for **Structure and exact algorithms for outer-connected semitotal double Roman domination in generalized coronas**, by Feng Sun.

For a connected core graph H with m vertices, positive multiplicities k[v] attach private leaves to each core vertex. The solver minimizes the correction sum of k[v] - 2 over feasible connected core sets, including the empty set, then adds the baseline sum(k) + 2m. Each selected single-leaf core must have an unselected core neighbour.

## Contents

| Path | Purpose |
| --- | --- |
| code/core_solvers.py | Tree DP, optional core-set counting, leaf response, subset enumeration and rooted core MILP |
| code/original_label_milp.py | Independent four-label MILP on the expanded graph |
| code/verify_ocstdr.py | Original graph-based feasibility predicate and private-leaf constructor |
| code/run_experiments.py | Exhaustive tree checks, atlas checks, generated general cores and large tree timings |
| code/check_response.py | Independent leaf-addition comparisons |
| code/make_tables.py | Generate the three manuscript tables from JSON results |
| results/ | Original experiment records, including small inputs and raw timings |
| tables/ | Tables generated from those original records |
| example.py | Small worked example |

## Environment

The recorded experiments used Python **3.8.5**, NumPy **1.24.4**, SciPy **1.10.1** and NetworkX **2.5** on Windows 10, with an Intel Core i7-10510U CPU and 16 GB RAM.

Use Python 3.8 to recreate the recorded dependency environment. The pinned legacy versions below are intended for reproduction of these experiments, rather than an untested current Python installation. SciPy supplies the HiGHS interface; no commercial optimizer is required.

From the repository root:

~~~sh
python -m venv .venv
~~~

Activate the environment:

~~~powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
~~~

~~~sh
# Linux / macOS
source .venv/bin/activate
~~~

Then install and run the example:

~~~sh
python -m pip install -r requirements.txt
python example.py
~~~

Expected: minimum weight 10, correction -1, and two minimizing core sets.

## Recreate tables from saved records

This does not rerun the optimization experiments:

~~~sh
python code/make_tables.py
~~~

The generated LaTeX tables are placed in outputs/tables/. The committed tables/ directory contains the corresponding reference tables. To compile these table snippets in a separate LaTeX document, load the float package for their [H] placement.

## Rerun the complete experiments

~~~sh
python code/run_experiments.py
python code/check_response.py
python code/make_tables.py --results-dir outputs/results --output-dir outputs/tables
~~~

New results go to outputs/results/ and leave the original results/ files intact. Repeating these commands replaces the corresponding files under outputs/. The full experiment includes the 100,000-vertex tree cases; running time depends on the machine.

You can choose different destinations:

~~~sh
python code/run_experiments.py --output-dir my-run/results
python code/check_response.py --output-dir my-run/results
python code/make_tables.py --results-dir my-run/results --output-dir my-run/tables
~~~

## What is checked

- 5,304 tree inputs: DP values and minimizing core-set counts versus subset enumeration.
- 750 additional tree-root checks and 30 edge-case runs.
- 92 atlas inputs and 24 generated general-core inputs: rooted core MILP versus subsets.
- 29 expanded four-label MILP checks and 14 direct label-enumeration checks.
- 11,720 leaf-addition predictions versus fresh enumeration.
- Four tree families at core orders 100, 1,000, 10,000 and 100,000.

These counts refer to different, overlapping input collections; they are not a count of distinct instances obtained by summing all rows.

## Conventions and timing

Core vertices are indexed 0,...,m-1. Inputs are finite, simple, undirected connected graphs with positive integer multiplicities. tree_dp requires a tree. The zero-labelled set may be empty; otherwise it must induce a connected graph.

Optional counting returns the number of minimizing **core sets**, not the number of all minimum four-label functions. The expanded witness is constructed only for small checks.

Large tree timings include rooting, value computation, core-set recovery and post-solution checks. They exclude graph generation and expanded-graph materialization, and disable counting. Each table cell is the median of 15 calls: three weighted inputs, five calls each, after one warm-up per input. Times will vary; saved reference timings are measurements from the recorded run.

Random seeds and draw ordering are in the scripts: master seed 20261020, generated general cores use seed + 1, and tree timing uses seed + 2. Small inputs are saved explicitly. Large input hashes are stored with raw timings. Candidate-root MILPs have a 30-second time limit; the expanded-label MILP has a 45-second limit. Unsuccessful or time-limited feasible solves raise an error rather than being reported as exact.

## Citation

Citation metadata is provided in CITATION.cff. Repository: https://github.com/sunsz2005/OCSTDR. No publication DOI has been assigned. No open-source license has been selected for this repository.

"""Worked path example: multiplicities (1, 3, 1)."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "code"))
from core_solvers import adjacency, tree_dp

answer = tree_dp(adjacency(3, [(0, 1), (1, 2)]), [1, 3, 1], count=True)
assert (answer["value"], answer["correction"], answer["count"]) == (10, -1, 2)
print("Minimum weight:", answer["value"])
print("Correction:", answer["correction"])
print("Minimizing core sets:", answer["count"])
print("Recovered core set:", sorted(answer["selected"]))

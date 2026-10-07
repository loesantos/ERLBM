#!/usr/bin/env python3
"""Plot the double-shear-layer convergence (Fig. 4 of the E-RLBM paper).

Usage:
    python3 scripts/plot_fig4_double_shear_convergence.py [results/fig4_double_shear_Re80000.csv] [--norm vector|magnitude]

Writes results/fig4_double_shear_convergence.pdf (and .png).
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if not a.startswith("--")]
csv_file = args[0] if args else "results/fig4_double_shear_Re80000.csv"
norm = "vector"
if "--norm" in sys.argv:
    norm = sys.argv[sys.argv.index("--norm") + 1]
column = "E_vector" if norm == "vector" else "E_magnitude"

data = defaultdict(list)
with open(csv_file) as f:
    for row in csv.DictReader(f):
        data[row["model"]].append((int(row["L"]), float(row[column])))

style = {
    "ERLBM065": dict(label=r"E-RLBM $\tau_{nh}=0.65$", marker="o", color="tab:blue", ls="-"),
    "ERLBM075": dict(label=r"E-RLBM $\tau_{nh}=0.75$", marker="^", color="tab:orange", ls="--"),
    "REG": dict(label="RLBM", marker="x", color="tab:red", ls="--"),
}

fig, ax = plt.subplots(figsize=(6.4, 4.8))
for model in ["ERLBM065", "ERLBM075", "REG"]:
    if model not in data:
        continue
    pts = sorted(data[model])
    ax.loglog([p[0] for p in pts], [p[1] for p in pts], lw=1, **style[model])

all_L = sorted({p[0] for v in data.values() for p in v})
if all_L:
    L0, L1 = all_L[0], all_L[-1]
    c = 1.5 * max(p[1] * p[0] ** 2 for v in data.values() for p in v if p[0] >= 256)
    ax.loglog([L0, L1], [c / L0 ** 2, c / L1 ** 2], "k--", lw=0.8, label=r"$L^{-2}$")
    ax.set_xticks(all_L)
    ax.set_xticklabels([str(L) for L in all_L])

ax.set_xlabel("Number of lattice nodes per axis L")
ax.set_ylabel(r"Relative $L_2$ error norm")
ax.grid(True, which="both", ls=":", lw=0.5)
ax.legend()
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/fig4_double_shear_convergence.pdf")
fig.savefig("results/fig4_double_shear_convergence.png", dpi=150)
print("results/fig4_double_shear_convergence.pdf written ({} norm)".format(norm))

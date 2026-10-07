#!/usr/bin/env python3
"""Fig. 3 of the E-RLBM paper: decaying Taylor-Green vortex, Re = 50, diffusive scaling.

tau = 0.62 fixed (nu = 0.04), U0 = 2/L, error at t = T = L/U0, D2Q9, vector norm of Eq. (l2).
Initial density from the exact pressure Poisson solution; populations f_eq + f^(1).
Usage: python3 scripts/plot_fig3_taylor_green.py [results/fig3_taylor_green_Re50.csv]
"""
import csv, sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

csv_file = sys.argv[1] if len(sys.argv) > 1 else "results/fig3_taylor_green_Re50.csv"
data = defaultdict(dict)
for r in csv.DictReader(open(csv_file)):
    data[r["model"]][int(r["L"])] = float(r["E_vector"])
L = sorted(data["BGK"])
sty = [("ERLBM065", "o-", "blue", 6, None, r"E-RLBM $\tau_{nh}=0.65$"),
       ("REG", "^--", "red", 6, None, "RLBM"),
       ("BGK", "s:", "gray", 5, "none", "BGK")]
fig, ax = plt.subplots(figsize=(6.5, 4.95))
for m, fmt, col, ms, mfc, lab in sty:
    ax.loglog(L, [data[m][x] for x in L], fmt, color=col, lw=0.8, ms=ms, mfc=mfc, label=lab)
c = 1.6 * max(data[m][x] * x ** 2 for m in data for x in L)
ax.loglog([L[0], L[-1]], [c / L[0] ** 2, c / L[-1] ** 2], "k--", lw=0.8)
ax.set_xticks(L); ax.set_xticklabels([str(x) for x in L])
ax.minorticks_on(); ax.xaxis.set_minor_formatter(NullFormatter())
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.set_xlabel("Number of lattice nodes per axis L", fontsize=12)
ax.set_ylabel(r"Relative $L_2$ - error norm", fontsize=12)
ax.grid(True, which="both", ls="--", lw=0.5, color="gray", alpha=0.6)
ax.legend(fontsize=11)
fig.tight_layout()
import os; os.makedirs("results", exist_ok=True)
fig.savefig("results/fig3_taylor_green_convergence.pdf"); fig.savefig("results/fig3_taylor_green_convergence.png", dpi=150)

#!/usr/bin/env python3
"""Verificacao do coeficiente hiperviscoso diagonal kappa_4(pi/4) = a^2 ( b - 2a ) / 18  ( Eq. 16 do artigo ).

    python3 scripts/plot_kappa4_check.py [results/kappa4_check_L64.csv]

Pontos: ( nu_ef - nu ) / nu medido com populacoes iniciais f_eq + f^(1) ( init = neq ).
Linhas: previsao linear a ( 2a - b ) K^2 / 6,  a = tau - 1/2,  b = tau_nh - 1/2,  K^2 = 8 pi^2 / L^2.
Imprime a razao medido / previsto.   Grava results/kappa4_check_L<L>.pdf
"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

f = sys.argv[1] if len(sys.argv) > 1 else "results/kappa4_check_L64.csv"
rows = [r for r in csv.DictReader(open(f)) if r["init"] == "neq"]
L = int(rows[0]["L"])
K2 = 8 * np.pi ** 2 / L ** 2
bnh = {"BGK": None, "REG": 0.5, "ERLBM065": 0.15}
sty = {"ERLBM065": ("o", "-", "blue", r"E-RLBM $\tau_{nh}=0.65$"),
       "REG": ("^", "--", "red", "RLBM"), "BGK": ("s", ":", "gray", "BGK")}

t = np.linspace(0.5, 1.2, 300)
a = t - 0.5
fig, ax = plt.subplots(figsize=(6.5, 4.95))
for m in ["ERLBM065", "REG", "BGK"]:
    mk, ls, col, lab = sty[m]
    b = a if bnh[m] is None else bnh[m]
    ax.plot(t, 1e4 * a * (2 * a - b) * K2 / 6, ls, color=col, lw=0.9)
    pts = sorted((float(r["tau"]), float(r["rel_visc_error"])) for r in rows if r["model"] == m)
    x, y = np.array(pts).T
    ax.plot(x, 1e4 * y, mk, color=col, ms=6, mfc="none" if m == "BGK" else col, ls="none", label=lab)
    aa = x - 0.5
    bb = aa if bnh[m] is None else bnh[m]
    pred = aa * (2 * aa - bb) * K2 / 6
    print(m, " ".join("%.3f:%.4f" % (xi, yi / pi) if abs(pi) > 1e-12 else "%.3f:--" % xi
                     for xi, yi, pi in zip(x, y, pred)))
ax.axhline(0, color="k", lw=0.5)
ax.axvline(0.65, color="k", lw=0.5, ls="--")
ax.set_xlim(0.5, 1.2)
ax.minorticks_on()
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.set_xlabel(r"$\tau$", fontsize=12)
ax.set_ylabel(r"$10^4\,(\nu_{\mathrm{ef}}-\nu)/\nu$", fontsize=12)
ax.grid(True, which="major", ls="--", lw=0.5, color="gray", alpha=0.6)
ax.legend(fontsize=11, loc="upper left")
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/kappa4_check_L%d.pdf" % L)
fig.savefig("results/kappa4_check_L%d.png" % L, dpi=150)

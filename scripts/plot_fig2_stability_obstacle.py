#!/usr/bin/env python3
"""Fig. 2: mapa de estabilidade do escoamento em torno do obstaculo ( E-RLBM, D3Q19, Ma = 0.1 ).

    python3 scripts/plot_fig2_stability_obstacle.py [data/fig2_map_obstacle.vtk] [results/fig2_stability_obstacle.pdf]

Branco: estavel; preto: instavel.  Linha vermelha: RLBM ( tau_nh = 1 ); curva azul tracejada: BGK ( tau_nh = tau ).
Eixos do mapa: tau_nh = 0.5 + 0.6 i / 400,  log10 nu = -5 + 4.3 j / 400.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["DejaVu Serif"]})

entrada = sys.argv[1] if len(sys.argv) > 1 else "data/fig2_map_obstacle.vtk"
saida = sys.argv[2] if len(sys.argv) > 2 else "results/fig2_stability_obstacle.pdf"
lines = open(entrada).read().split("\n")
nx, ny = [int(v) for v in next(l for l in lines if l.startswith("DIMENSIONS")).split()[1:3]]
k = next(n for n, l in enumerate(lines) if l.startswith("LOOKUP_TABLE"))
m = np.array(" ".join(lines[k + 1:]).split(), dtype=float)[: nx * ny].reshape(ny, nx)

dt, dl = 0.6 / nx, 4.3 / ny
tnh = 0.5 + 0.6 * np.arange(nx) / nx
lognu = -5 + 4.3 * np.arange(ny) / ny

fig, ax = plt.subplots(figsize=(6.5, 4.8))
ax.imshow(m, origin="lower", cmap=ListedColormap(["white", "black"]), vmin=0, vmax=1, interpolation="nearest",
          extent=(0.5 - dt / 2, 1.1 - dt / 2, -5 - dl / 2, -0.7 - dl / 2), aspect="auto", rasterized=True)
t = np.linspace(0.50002, 1.1, 2000)
ax.plot(t, np.log10((t - 0.5) / 3), "--", color="#3466a8", lw=1.8)
ax.axvline(1.0, color="#e0201b", ls="--", lw=1.6)
ax.set_xlim(0.5, 1.1)
ax.set_ylim(-5, -0.7)
ax.set_xticks(np.round(np.arange(0.5, 1.101, 0.1), 1))
ax.set_yticks(np.round(np.arange(-4.8, -0.79, 0.5), 1))
ax.set_xlabel(r"$\tau_{nh}$", fontsize=16)
ax.set_ylabel(r"$\log(\nu)$", fontsize=16)
ax.tick_params(labelsize=11)
fig.tight_layout()
import os
os.makedirs(os.path.dirname(saida) or ".", exist_ok=True)
fig.savefig(saida, dpi=300)
fig.savefig(saida.replace(".pdf", ".png"), dpi=110)

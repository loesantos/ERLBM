#!/usr/bin/env python3
"""Fig. 1: mapa de estabilidade do vortice de Taylor-Green, RLBM x E-RLBM ( tau_nh = 0.65 ), D3Q19, L = 16.

    python3 scripts/plot_fig1_stability_taylor_green.py [data/fig1_map_reg.vtk] [data/fig1_map_erlbm065.vtk] [results/fig1_stability_taylor_green.pdf]

Azul: os dois estaveis;  vermelho: so o E-RLBM estavel;  preto: os dois instaveis;  amarelo: so o RLBM estavel.
O detalhe ampliado mostra o canto de Ma e nu pequenos, onde estao os pontos em que so o RLBM e estavel.
Eixos do mapa: Ma = 0.05 + ( i + 1 ) 0.95 / 400,   log10 nu = -5 + 3 j / 400.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["DejaVu Serif"]})


def le_mapa(path):
    lines = open(path).read().split("\n")
    nx, ny = [int(v) for v in next(l for l in lines if l.startswith("DIMENSIONS")).split()[1:3]]
    k = next(n for n, l in enumerate(lines) if l.startswith("LOOKUP_TABLE"))
    return np.array(" ".join(lines[k + 1:]).split(), dtype=float)[: nx * ny].reshape(ny, nx)


a = sys.argv[1:]
reg = le_mapa(a[0] if len(a) > 0 else "data/fig1_map_reg.vtk")
erl = le_mapa(a[1] if len(a) > 1 else "data/fig1_map_erlbm065.vtk")
saida = a[2] if len(a) > 2 else "results/fig1_stability_taylor_green.pdf"
ny, nx = reg.shape
# 0 ambos estaveis, 1 so E-RLBM, 2 ambos instaveis, 3 so RLBM
cls = np.select([(reg == 0) & (erl == 0), (reg == 1) & (erl == 0), (reg == 1) & (erl == 1)], [0, 1, 2], 3)
print("ambos estaveis %d, so E-RLBM %d, ambos instaveis %d, so RLBM %d" % tuple((cls == c).sum() for c in range(4)))

cmap = ListedColormap(["#3b4cc0", "#b40426", "#1a1a1a", "#ffd21f"])
dx, dl = 0.95 / nx, 3.0 / ny
ext = (0.05 + dx / 2, 1.0 + dx / 2, -5 - dl / 2, -2 - dl / 2)
opts = dict(origin="lower", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest", extent=ext, aspect="auto")

fig, ax = plt.subplots(figsize=(6.5, 4.8))
ax.imshow(cls, rasterized=True, **opts)
ax.set_xlim(0.05, 1.0)
ax.set_ylim(-5, -2)
ax.set_xticks(np.round(np.arange(0.05, 0.96, 0.1), 2))
ax.set_yticks(np.arange(-5, -1.99, 0.5))
ax.grid(True, ls=":", lw=0.5, color="white", alpha=0.5)
ax.set_xlabel(r"$Ma$", fontsize=16)
ax.set_ylabel(r"$\log(\nu)$", fontsize=16)
ax.tick_params(labelsize=11)

# detalhe: canto de Ma e nu pequenos
axi = ax.inset_axes([0.12, 0.62, 0.25, 0.35])
axi.imshow(cls, rasterized=True, **opts)
axi.set_xlim(0.05, 0.20)
axi.set_ylim(-5, -4.3)
axi.set_xticks([0.05, 0.10, 0.15])
axi.set_yticks([-5.0, -4.8, -4.6, -4.4])
axi.tick_params(labelsize=8, colors="white", direction="in")
for s in axi.spines.values():
    s.set_edgecolor("white")
mark_inset(ax, axi, loc1=3, loc2=4, fc="none", ec="white", lw=0.8)

fig.tight_layout()
import os
os.makedirs(os.path.dirname(saida) or ".", exist_ok=True)
fig.savefig(saida, dpi=300)
fig.savefig(saida.replace(".pdf", ".png"), dpi=110)

#!/usr/bin/env python3
"""Fig. 5: evolucao da vorticidade na camada de cisalhamento dupla.

    python3 scripts/plot_fig5_double_shear_evolution.py runs/ERLBM075_L512 512

Le vor_<passo>.vtk em t/T = 0.50, 0.75, 1.00 e 1.25 ( passo = t/T * 32 L ) e grava
results/fig5_double_shear_evolution.pdf ( mapa de cores divergente e isolinhas da vorticidade ).
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtkmap import read_vectors

pasta, L = sys.argv[1], int(sys.argv[2])
tempos = [0.50, 0.75, 1.00, 1.25]
fig, axs = plt.subplots(2, 2, figsize=(9, 9.3))
for ax, t in zip(axs.flat, tempos):
    passo = int(round(t * 32 * L))
    w = read_vectors(os.path.join(pasta, "vor_%d.vtk" % passo))[:, :, 2]
    lim = np.abs(w).max()
    ax.imshow(w, origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="bilinear")
    ax.contour(w, levels=np.linspace(-lim, lim, 16), colors="k", linewidths=0.3)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("t = %.2f T" % t, fontsize=13)
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/fig5_double_shear_evolution.pdf", dpi=200)
fig.savefig("results/fig5_double_shear_evolution.png", dpi=110)

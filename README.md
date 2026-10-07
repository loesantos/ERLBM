# E-RLBM — extended-regularized lattice Boltzmann method

Source code for the simulations of

> L. O. E. dos Santos, D. N. Siebert, R. L. M. Bazarin,
> *Extended-regularized lattice Boltzmann method*, submitted to *Computers & Fluids* (2026).

Every program is built on **SimBoltz**, a C++ library of lattice Boltzmann functions. A snapshot of the
library files used here is included in `SimBoltz_Functions/`.

## The collision operator

The E-RLBM splits the nonequilibrium population into the regularized (second-order Hermite) part `f1`
and its complement `f_gt1 = f - f_eq - f1`. It relaxes the two parts with different times:

    f* = f_eq + (1 - 1/tau) f1 + (1 - 1/tau_nh) f_gt1

- `tau` sets the viscosity, `nu = (tau - 1/2)/3`.
- `tau_nh` relaxes the kinetic (non-hydrodynamic) sector.
- `tau_nh = 1` gives the regularized LBM (RLBM) of Latt & Chopard.
- `tau_nh = tau` gives BGK.

In the library the operator is `coll_Ext_Reg()` in `SimBoltz_Functions/Collision_operators.cpp`. Earlier
versions of the code called it *over-regularized* (`coll_BGK_OverReg`).

## Contents

| path | content |
|---|---|
| `src/main_stability-map_taylor-green.cpp` | Fig. 1: Taylor–Green stability map (D3Q19) |
| `src/main_stability-map_obstacle.cpp` | Fig. 2: stability map of the flow past an obstacle (D3Q19) |
| `src/main_taylor-green.cpp` | Fig. 3: Taylor–Green convergence; check of the diagonal hyperviscous coefficient (D2Q9) |
| `src/main_double-shear-layer.cpp` | Figs. 4 and 5: double shear layer (D2Q9; OpenMP or OpenACC) |
| `Makefile` | build rules |
| `src/erlbm_common.h` | code shared by the programs (library includes, collision choice, geometry, command line) |
| `tools/error_calculator.cpp` | relative L2 error of a velocity field against a finer reference (Fig. 4) |
| `scripts/run_*.sh` | run the cases of each figure and write the data to `results/` |
| `scripts/plot_*.py` | plot the figures from those data |
| `analysis/` | linear and symbolic analysis of D2Q9 (Eqs. 15–17) and the Hermite orthogonality check (Sec. 2.2) |
| `data/` | the data plotted in the paper |
| `SimBoltz_Functions/` | SimBoltz library (snapshot) |

Comments in the source code are in Portuguese.

## Building

Requirements:

- a C++17 compiler with OpenMP (tested with g++ 13);
- Python 3 with NumPy and Matplotlib for the plots;
- SymPy for `analysis/kappa4_symbolic.py`.

    make            # bin/stability_map_taylor_green, bin/stability_map_obstacle, bin/taylor_green,
                    # bin/taylor_green_d3q19, bin/double_shear_layer, bin/error_calculator
    make gpu        # bin/double_shear_layer_gpu (OpenACC, NVIDIA HPC SDK nvc++)

The library is taken from `SimBoltz_Functions/`. To use another copy, run `make SIMBOLTZ=/path/to/SimBoltz_Functions`
(the path must not contain spaces).

Every program prints its options in the header comment of its source file. Options use the form `--name value`.

## Reproducing the figures

All scripts run from any directory and write to `results/`.

| figure | command | cost |
|---|---|---|
| Fig. 1 | `scripts/run_fig1_stability_taylor_green.sh` | 2 × 160 000 runs (cluster); `N=40` for a quick test |
| Fig. 2 | `scripts/run_fig2_stability_obstacle.sh` | 160 000 runs (cluster); `N=40` for a quick test |
| Fig. 3 | `scripts/run_fig3_taylor_green.sh` | about 10 min |
| Sec. 4.2, Eq. (16) check | `scripts/run_kappa4_check.sh` | about 5 min |
| Fig. 4 | `scripts/run_fig4_double_shear_layer.sh reference`, then `scripts/run_fig4_double_shear_layer.sh` | reference L = 2048: hours |
| Fig. 5 | `scripts/run_fig4_double_shear_layer.sh fig5` | L = 512 |

The scripts take `THREADS` (OpenMP threads, 0 = all) and, where relevant, `N` (map points per axis) or `LS`
(grid sizes).

### Fig. 1 — Taylor–Green stability map

- D3Q19 with one lattice layer in z, L = 16, periodic.
- Equilibrium initialization with uniform density.
- The map covers 400 × 400 points:
  - `U0 = 0.05 cs + (i + 1)(cs - 0.05 cs)/400`, which gives 0.05 < Ma ≤ 1;
  - `log10 nu = -5 + 3 j/400`.
- A run is unstable if the total kinetic energy exceeds twice its initial value at a check. Checks are
  made every 100 steps during 100 T.
- Classification of the maps:
  - **Blue:** RLBM and E-RLBM (`tau_nh = 0.65`) both stable.
  - **Red:** only the E-RLBM stable.
  - **Black:** both unstable.
  - **Yellow:** only the RLBM stable. These are 451 points of the paper's map, all in the scattered region
    at Ma < 0.13 and log10 nu < -4.4, which the inset of the figure enlarges.

### Fig. 2 — obstacle stability map

- D3Q19 channel of 20 × 12 nodes with a 2 × 4 block, Ma = 0.1. The block is centred vertically, with 6 fluid
  columns upstream and 12 downstream.
- Boundary conditions:
  - Zou–He velocity inlet;
  - outlet populations copied from the upstream neighbour and rescaled to the initial density;
  - halfway bounce-back on the block;
  - periodic in y.
- The map covers `tau_nh = 0.5 + 0.6 i/400` and `log10 nu = -5 + 4.3 j/400`.
- A run lasts `200 ny/U0` steps. It is unstable if, at any step, the density of any node exceeds ten times the
  initial value or is not finite.

### Fig. 3 — Taylor–Green convergence

- D2Q9, Re = 50, diffusive scaling: `tau = 0.62` fixed and `U0 = 2/L`.
- Error at `t = T = L/U0`.
- Initial density from the exact pressure Poisson solution.
- Initial populations `f_eq + f^(1)`. `INIT=eq` gives equilibrium populations only.

The same program, run with `--tau` instead of `--Re`, measures the error when the amplitude has decayed to `e^-1`
(`2 nu k^2 t = 1`). In that mode it also prints the relative error of the effective viscosity. This is the check of
`kappa_4(pi/4) = a^2 (b - 2a)/18` reported in Sec. 4.2 (`scripts/run_kappa4_check.sh`).

### Figs. 4 and 5 — double shear layer

- D2Q9, Re = 80000, `U0 = 1/32`, `lambda = 80`, `epsilon = 0.05`.
- Uniform density and equilibrium populations.
- The error is computed at `t = T` against an RLBM run with L = 2048. Each grid's nodes coincide with every
  `(2048/L)`-th reference node, so no interpolation is used.
- `tools/error_calculator` prints two norms:
  - the vector norm of the paper, `E_vector = sqrt( sum |u - u_ref|^2 / sum |u_ref|^2 )`;
  - the magnitude norm `E_magnitude = sqrt( sum (|u| - |u_ref|)^2 / sum |u_ref|^2 )`.
- Fig. 5 is the E-RLBM run with `tau_nh = 0.75` at L = 512.

## Data used in the paper (`data/`)

| file | content |
|---|---|
| `fig1_map_reg.vtk`, `fig1_map_erlbm065.vtk` | Fig. 1 maps |
| `fig3_taylor_green_Re50.csv` | Fig. 3 |
| `fig3_taylor_green_Re50_feq_init.csv` | the same runs with equilibrium initialization (for comparison) |
| `fig4_double_shear_Re80000.csv` | Fig. 4 |
| `kappa4_check_L64.csv` | check of Eq. (16) |
| `fig2_map_obstacle.vtk` | Fig. 2 map |

### Previous version of the code

The Fig. 1 maps were computed with an earlier version of the code, which had its own copy of the library
functions. The programs in `src/` reproduce them except for points near the stability boundary:

- **Test:** a 20 × 20 subsample of each map.
- **Result:** 2 differing points out of 400 for the E-RLBM map and 5 out of 400 for the RLBM map. All of them
  are in mixed regions of the original map.
- **Why points differ:** near the boundary, an instability grows so slowly that it reaches the threshold close
  to the end of the 100 T horizon, and whether it does depends on rounding (compiler and flags).
- **NaN check:** the earlier code did not test for non-finite energy. No run of the earlier code became non-finite
  at the 695 stable points surrounded by unstable ones.

The Fig. 2 map in `data/` was computed on a cluster with the earlier code, after two corrections: its instability
test now also counts a non-finite density as unstable, and a race condition in its OpenMP loop was removed. The
corrected map is identical, point by point, to the map computed before the corrections.
`src/main_stability-map_obstacle.cpp` uses the same test. On a 20 × 20 subsample it agrees with the map at 398 of
400 points; the 2 differing points lie on the scattered lower boundary of the stable region.

The convergence data (Figs. 3 and 4) and the Eq. (16) check are reproduced exactly by the current programs.
For Fig. 4, L = 64 gives identical output files; the L = 1024 E-RLBM `tau_nh = 0.75` value was obtained on a
cluster with the earlier code.

## License GPL-3.0

<!-- choose a license, e.g. MIT or GPL-3.0 -->

## Citation

<!-- add the DOI of the paper / Zenodo archive here -->

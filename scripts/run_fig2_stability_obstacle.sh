#!/bin/bash
#=====================================================================================================================#
#   Fig. 2 - mapa de estabilidade do escoamento em torno de um obstaculo ( E-RLBM, tau_nh x nu ), D3Q19, Ma = 0.1.
#
#       scripts/run_fig2_stability_obstacle.sh
#
#   Variaveis:  N ( pontos por eixo, padrao 400 )   THREADS ( 0 = todos )
#   Saida:      results/fig2_map_obstacle.vtk
#=====================================================================================================================#
set -e
cd "$( dirname "$0" )/.."
MK=$( [ -f Makefile ] && echo Makefile || echo erlbm.mk )
N=${N:-400}
THREADS=${THREADS:-0}
make -f "$MK" bin/stability_map_obstacle > /dev/null
mkdir -p results
./bin/stability_map_obstacle --n-tau $N --n-nu $N --threads $THREADS --out results/fig2_map_obstacle.vtk
python3 scripts/plot_fig2_stability_obstacle.py results/fig2_map_obstacle.vtk

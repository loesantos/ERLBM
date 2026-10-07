#!/bin/bash
#=====================================================================================================================#
#   Fig. 1 - mapa de estabilidade do vortice de Taylor-Green ( RLBM e E-RLBM com tau_nh = 0.65 ), D3Q19, L = 16.
#
#       scripts/run_fig1_stability_taylor_green.sh
#
#   Variaveis:  N ( pontos por eixo, padrao 400 )   THREADS ( 0 = todos )
#   Saidas:     results/fig1_map_reg.vtk, results/fig1_map_erlbm065.vtk
#
#   Custo: o mapa 400 x 400 tem 160 000 simulacoes de ate ~55 000 passos cada; no artigo ele foi feito num
#   cluster.  Use N=40 para um teste rapido.
#=====================================================================================================================#
set -e
cd "$( dirname "$0" )/.."
MK=$( [ -f Makefile ] && echo Makefile || echo erlbm.mk )
N=${N:-400}
THREADS=${THREADS:-0}
make -f "$MK" bin/stability_map_taylor_green > /dev/null
mkdir -p results
./bin/stability_map_taylor_green --collision reg   --n-ma $N --n-nu $N --threads $THREADS --out results/fig1_map_reg.vtk
./bin/stability_map_taylor_green --collision erlbm --tau-nh 0.65 --n-ma $N --n-nu $N --threads $THREADS --out results/fig1_map_erlbm065.vtk
python3 scripts/plot_fig1_stability_taylor_green.py results/fig1_map_reg.vtk results/fig1_map_erlbm065.vtk

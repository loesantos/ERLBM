#!/bin/bash
#=====================================================================================================================#
#   Fig. 3 - convergencia no vortice de Taylor-Green: Re = 50, tau = 0.62 fixo, U_0 = 2/L, t = T, D2Q9.
#
#       scripts/run_fig3_taylor_green.sh            populacoes iniciais f_eq + f^(1)  ( artigo )
#       INIT=eq scripts/run_fig3_taylor_green.sh    so f_eq ( para comparacao )
#
#   Saida: results/fig3_taylor_green_Re50.csv  ( ou results/fig3_taylor_green_Re50_feq_init.csv )
#=====================================================================================================================#
set -e
cd "$( dirname "$0" )/.."
MK=$( [ -f Makefile ] && echo Makefile || echo erlbm.mk )
INIT=${INIT:-neq}
make -f "$MK" bin/taylor_green > /dev/null
mkdir -p results
CSV=results/fig3_taylor_green_Re50.csv
[ "$INIT" = "eq" ] && CSV=results/fig3_taylor_green_Re50_feq_init.csv
echo "model,L,tau,E_vector,E_magnitude" > "$CSV"
for L in 32 64 128 256; do
	for m in "BGK bgk 1.0" "REG reg 1.0" "ERLBM065 erlbm 0.65"; do
		set -- $m
		read -r l tau tnh u0 steps ev em dnu < <( ./bin/taylor_green --L $L --collision $2 --tau-nh $3 --init $INIT )
		echo "$1,$l,$tau,$ev,$em" | tee -a "$CSV"
	done
done
python3 scripts/plot_fig3_taylor_green.py "$CSV"

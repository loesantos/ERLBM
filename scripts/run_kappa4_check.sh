#!/bin/bash
#=====================================================================================================================#
#   Verificacao do coeficiente hiperviscoso diagonal, Eq. (16):  ( nu_ef - nu ) / nu = a ( 2a - b ) K^2 / 6
#   Taylor-Green, L = 64, U_0 = 1e-3, erro medido em 2 nu k^2 t = 1, varredura em tau ( Sec. 4.2 do artigo ).
#
#       scripts/run_kappa4_check.sh
#
#   Saida: results/kappa4_check_L64.csv   ( colunas: model, init, L, tau, tau_nh, U0, steps, E_vector,
#          E_magnitude, rel_visc_error )  e  results/kappa4_check_L64.pdf
#=====================================================================================================================#
set -e
cd "$( dirname "$0" )/.."
MK=$( [ -f Makefile ] && echo Makefile || echo erlbm.mk )
L=${L:-64}
make -f "$MK" bin/taylor_green > /dev/null
mkdir -p results
CSV=results/kappa4_check_L$L.csv
echo "model,init,L,tau,tau_nh,U0,steps,E_vector,E_magnitude,rel_visc_error" > "$CSV"
for tau in 0.505 0.51 0.52 0.535 0.55 0.575 0.6 0.625 0.65 0.7 0.75 0.8 0.9 1.0 1.1 1.2; do
	for m in "BGK bgk $tau" "REG reg 1.0" "ERLBM065 erlbm 0.65"; do
		set -- $m
		for init in eq neq; do
			echo "$1,$init,$( ./bin/taylor_green --L $L --tau $tau --u0 0.001 --collision $2 --tau-nh $3 --init $init | tr ' ' ',' )" >> "$CSV"
		done
	done
	echo "   tau = $tau"
done
python3 scripts/plot_kappa4_check.py "$CSV"

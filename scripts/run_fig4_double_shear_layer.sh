#!/bin/bash
#=====================================================================================================================#
#   Figs. 4 e 5 - camada de cisalhamento dupla, Re = 80000, U_0 = 1/32, D2Q9.
#
#       scripts/run_fig4_double_shear_layer.sh reference     referencia RLBM L = 2048 ( varias horas )
#       scripts/run_fig4_double_shear_layer.sh               L = 64 ... 1024: RLBM e E-RLBM ( tau_nh = 0.65, 0.75 )
#       scripts/run_fig4_double_shear_layer.sh fig5          so o caso da Fig. 5 ( E-RLBM tau_nh = 0.75, L = 512 )
#
#   Variaveis:  LS ( padrao "64 128 256 512 1024" )   REF ( padrao runs/REF_L2048/vel_1.000000.vtk )
#               THREADS ( 0 = todos )
#
#   Cada caso roda em runs/<modelo>_L<L>/.  Os erros em t = T vao para results/fig4_double_shear_Re80000.csv
#=====================================================================================================================#
set -e
cd "$( dirname "$0" )/.."
MK=$( [ -f Makefile ] && echo Makefile || echo erlbm.mk )
RE=80000
LS=${LS:-"64 128 256 512 1024"}
THREADS=${THREADS:-0}
REF=${REF:-runs/REF_L2048/vel_1.000000.vtk}
make -f "$MK" bin/double_shear_layer bin/error_calculator > /dev/null

run_case ()     # nome colisao L tau_nh
{
	local name=$1 col=$2 L=$3 tnh=$4
	local dir=runs/${name}_L${L}
	#  tau = 1/2 + 3 U_0 L / Re  ( mesmo arredondamento usado nas rodadas do artigo )
	local tau
	tau=$( awk -v L="$L" -v Re="$RE" 'BEGIN { printf "%.12g", 0.5 + 3.0 * ( 1.0 / 32.0 ) * L / Re }' )
	mkdir -p "$dir"
	( cd "$dir" && ../../bin/double_shear_layer --L "$L" --collision "$col" --tau-nh "$tnh" --tau "$tau" \
	                                           --threads "$THREADS" > log.txt 2>&1 )
	echo "   $dir  ( tau = $tau, tau_nh = $tnh )"
}

if [ "$1" = "reference" ]; then
	run_case REF reg 2048 1.0
	exit 0
fi

if [ "$1" = "fig5" ]; then
	run_case ERLBM075 erlbm 512 0.75
	python3 scripts/plot_fig5_double_shear_evolution.py runs/ERLBM075_L512 512
	exit 0
fi

[ -f "$REF" ] && REF=$( cd "$( dirname "$REF" )" && pwd )/$( basename "$REF" )
mkdir -p results
CSV=results/fig4_double_shear_Re80000.csv
echo "model,L,tau,tau_nh,E_vector,E_magnitude" > "$CSV"

for L in $LS; do
	for c in "REG reg 1.0" "ERLBM065 erlbm 0.65" "ERLBM075 erlbm 0.75"; do
		set -- $c
		run_case "$1" "$2" "$L" "$3"
		if [ -f "$REF" ]; then
			( cd "runs/$1_L$L" && ../../bin/error_calculator vel_1.000000.vtk "$REF" > erro.txt )
			read -r _ ev em < "runs/$1_L$L/erro_L2.dat"
			tau=$( awk -v L="$L" -v Re="$RE" 'BEGIN { printf "%.12g", 0.5 + 3.0 * ( 1.0 / 32.0 ) * L / Re }' )
			echo "$1,$L,$tau,$3,$ev,$em" >> "$CSV"
		else
			echo "      ( referencia $REF nao encontrada: erro nao calculado )"
		fi
	done
done
[ -f "$REF" ] && python3 scripts/plot_fig4_double_shear_convergence.py "$CSV"

#=====================================================================================================================#
#
#   E-RLBM - programas do artigo
#
#
#       make            compila tudo em bin/  ( CPU / OpenMP )
#       make gpu        double_shear_layer com OpenACC ( nvc++ )
#       make clean
#
#   Executaveis:
#
#       bin/stability_map_taylor_green    Fig. 1   ( D3Q19 )
#       bin/stability_map_obstacle        Fig. 2   ( D3Q19 )
#       bin/taylor_green                  Fig. 3 e verificacao de kappa_4   ( D2Q9 )
#       bin/taylor_green_d3q19            o mesmo em D3Q19 ( nz = 1 )
#       bin/double_shear_layer            Figs. 4 e 5   ( D2Q9 )
#       bin/error_calculator              erro L2 relativo contra o campo de referencia ( Fig. 4 )
#
#   SIMBOLTZ aponta para a biblioteca SimBoltz ( uma copia acompanha o repositorio ).
#
#=====================================================================================================================#

SIMBOLTZ ?= SimBoltz_Functions

CXX      ?= g++
NVCXX    ?= nvc++
GPUARCH  ?= ccnative
EXTRA    ?=

COMMON   := -std=c++17 -O3 -I"$(SIMBOLTZ)" -Isrc $(EXTRA)
CPUFLAGS := -march=native -fopenmp

LIB      := $(wildcard $(SIMBOLTZ)/*.cpp) src/erlbm_common.h

BINS     := bin/stability_map_taylor_green bin/stability_map_obstacle bin/taylor_green bin/taylor_green_d3q19 \
            bin/double_shear_layer bin/error_calculator

.PHONY: all gpu clean

all: $(BINS)

bin:
	mkdir -p bin

bin/stability_map_taylor_green: src/main_stability-map_taylor-green.cpp $(LIB) | bin
	$(CXX) $(COMMON) $(CPUFLAGS) $< -o $@

bin/stability_map_obstacle: src/main_stability-map_obstacle.cpp $(LIB) | bin
	$(CXX) $(COMMON) $(CPUFLAGS) $< -o $@

bin/taylor_green: src/main_taylor-green.cpp $(LIB) | bin
	$(CXX) $(COMMON) $(CPUFLAGS) -DLATT=9 $< -o $@

bin/taylor_green_d3q19: src/main_taylor-green.cpp $(LIB) | bin
	$(CXX) $(COMMON) $(CPUFLAGS) -DLATT=19 $< -o $@

bin/double_shear_layer: src/main_double-shear-layer.cpp $(LIB) | bin
	$(CXX) $(COMMON) $(CPUFLAGS) $< -o $@

bin/error_calculator: tools/error_calculator.cpp | bin
	$(CXX) -std=c++17 -O2 $< -o $@

gpu: src/main_double-shear-layer.cpp $(LIB) | bin
	$(NVCXX) $(COMMON) -acc=gpu -gpu=$(GPUARCH) -Minfo=accel $< -o bin/double_shear_layer_gpu

clean:
	rm -f $(BINS) bin/double_shear_layer_gpu

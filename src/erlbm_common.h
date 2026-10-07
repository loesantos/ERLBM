//=============================== E-RLBM - funcoes comuns aos programas do artigo ====================================//
//
//  Incluido por todos os main_*.cpp.  Antes de incluir este arquivo o programa define a rede:
//
//      #define nvel 9      #define dim 2       ( D2Q9 )
//      #define nvel 19     #define dim 3       ( D3Q19 )
//
//  Aqui ficam: a inclusao da biblioteca SimBoltz, a escolha do operador de colisao em tempo de
//  execucao, a montagem de geometrias simples na memoria ( sem arquivo .vtk de entrada ), a
//  alocacao da rede e uma leitura minima de argumentos da linha de comando.
//
//====================================================================================================================//

#ifndef ERLBM_COMMON_H
#define ERLBM_COMMON_H

#include "Definitions.cpp"

#include "LBM_functions.cpp"

#include "Boundary_conditions.cpp"

#include "Collision_operators.cpp"

#include "Initial_conditions.cpp"

#include "Other_functions.cpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>


//------ Operador de colisao ----------------------------------------------------------------------------------------//
//
//      bgk    ->  BGK                                     coll_BGK
//      mrt    ->  MRT ( d'Humieres )                      coll_MRT      ( so D2Q9 )
//      reg    ->  regularizado ( Latt & Chopard 2006 )    coll_Reg      ( = E-RLBM com tau_nh = 1 )
//      erlbm  ->  regularizado estendido ( E-RLBM )       coll_Ext_Reg  ( antigo "over-regularized" )

enum COLISAO { COL_BGK = 1, COL_MRT = 2, COL_REG = 3, COL_EXT_REG = 4 };

inline int le_colisao ( const string& nome )
{
	if ( nome == "bgk" )   return COL_BGK;
	if ( nome == "mrt" )   return COL_MRT;
	if ( nome == "reg" )   return COL_REG;
	if ( nome == "erlbm" ) return COL_EXT_REG;

	cerr << "\nOperador de colisao desconhecido: " << nome << "  ( use bgk, mrt, reg ou erlbm )" << endl;

	exit ( 1 );
}

inline const char* nome_colisao ( int col )
{
	switch ( col )
	{
		case COL_BGK:     return "BGK";
		case COL_MRT:     return "MRT";
		case COL_REG:     return "RLBM";
		case COL_EXT_REG: return "E-RLBM";
	}
	return "?";
}

#pragma acc routine seq
inline void colisao_sitio ( double *f, double acc_x, double acc_y, double acc_z,
                            LATTICE lattice, PARAMETERS parameters, int col )
{
	switch ( col )
	{
		case COL_BGK:     coll_BGK     ( f, acc_x, acc_y, acc_z, lattice, parameters ); break;
#if nvel == 9
		case COL_MRT:     coll_MRT     ( f, acc_x, acc_y, acc_z, lattice, parameters ); break;
#endif
		case COL_REG:     coll_Reg     ( f, acc_x, acc_y, acc_z, lattice, parameters ); break;
		case COL_EXT_REG: coll_Ext_Reg ( f, acc_x, acc_y, acc_z, lattice, parameters ); break;
	}
}


//------ Geometria na memoria ---------------------------------------------------------------------------------------//
//
//  geometry.ini[ x + y nx + z nx ny ] = 0 para solido e = indice ( 1, 2, ... ) do sitio fluido, na mesma ordem
//  ( x mais rapido ) usada por read_geo().  As fronteiras sao periodicas em def_dir_prop() .

inline void numera_fluido ( GEOMETRY& geometry )
{
	int n = 0;

	for ( int i = 0; i < geometry.nx * geometry.ny * geometry.nz; i++ )
		if ( geometry.ini[i] ) geometry.ini[i] = ++n;

	geometry.fluid = n;

	geometry.phi = ( double ) n / ( double ) ( geometry.nx * geometry.ny * geometry.nz );
}

inline void geometria_periodica ( GEOMETRY& geometry, int nx, int ny, int nz = 1 )
{
	geometry.file = "";
	geometry.ftesc = 1.0;
	geometry.nx = nx;
	geometry.ny = ny;
	geometry.nz = nz;

	geometry.ini = new int[ nx * ny * nz ];

	for ( int i = 0; i < nx * ny * nz; i++ ) geometry.ini[i] = 1;

	numera_fluido ( geometry );
}


//------ Rede -------------------------------------------------------------------------------------------------------//
//
//  Define as velocidades, os pesos e o mapa de propagacao ( ini_stream ).  As populacoes ( inif, inif_new )
//  sao alocadas a parte, porque os mapas de estabilidade rodam varias simulacoes em paralelo sobre a mesma
//  geometria: c_i, Q_i e ini_stream sao so lidos e podem ser compartilhados entre as threads.

inline void define_rede ( GEOMETRY& geometry, LATTICE& lattice )
{
	lattice.c_i = new double[ nvel * dim ];

#if nvel == 9
	def_lattice_d2q9 ( lattice );
#elif nvel == 19
	def_lattice_d3q19 ( lattice );
#else
	#error "rede nao prevista: use D2Q9 ( nvel 9 ) ou D3Q19 ( nvel 19 )"
#endif

	lattice.Q_i = new double[ nvel * dim * dim ];

	calc_Q ( lattice );

	lattice.ini_stream = new int[ geometry.fluid * nvel ];

	def_dir_prop ( geometry, lattice );
}

inline void aloca_populacoes ( const GEOMETRY& geometry, LATTICE& lattice )
{
	lattice.inif     = new double[ geometry.fluid * nvel ];
	lattice.inif_new = new double[ geometry.fluid * nvel ];
}

inline void libera_populacoes ( LATTICE& lattice )
{
	delete[] lattice.inif;
	delete[] lattice.inif_new;

	lattice.inif = nullptr;
	lattice.inif_new = nullptr;
}


//------ Linha de comando -------------------------------------------------------------------------------------------//
//
//  Argumentos na forma  --nome valor .  arg_str / arg_double / arg_int devolvem o valor ou o padrao.

struct ARGUMENTOS
{
	int argc;
	char** argv;

	const char* procura ( const char* nome ) const
	{
		for ( int i = 1; i < argc - 1; i++ )
			if ( strcmp ( argv[i], nome ) == 0 ) return argv[ i + 1 ];
		return nullptr;
	}

	bool tem ( const char* nome ) const
	{
		for ( int i = 1; i < argc; i++ )
			if ( strcmp ( argv[i], nome ) == 0 ) return true;
		return false;
	}

	string arg_str ( const char* nome, const string& padrao ) const
	{
		const char* v = procura ( nome );
		return v ? string ( v ) : padrao;
	}

	double arg_double ( const char* nome, double padrao ) const
	{
		const char* v = procura ( nome );
		return v ? atof ( v ) : padrao;
	}

	int arg_int ( const char* nome, int padrao ) const
	{
		const char* v = procura ( nome );
		return v ? atoi ( v ) : padrao;
	}
};

inline void define_threads ( int n )
{
	if ( n > 0 ) omp_set_num_threads ( n );
}

#endif

//=============================== Mapa de estabilidade - vortice de Taylor-Green ( Fig. 1 ) =========================//
//
//  Vortice de Taylor-Green decaindo, dominio periodico L x L x 1, rede D3Q19 ( uma camada em z ).
//
//      u_x =  U_0 sin( k x ) cos( k y ),    u_y = - U_0 cos( k x ) sin( k y ),    k = 2 pi / L
//
//  Populacoes iniciais no equilibrio, densidade uniforme.  Para cada par ( Ma, nu ) do mapa a simulacao corre
//  por 100 T, T = L / U_0, e e classificada como instavel se a energia cinetica total ultrapassar 2 vezes o
//  valor inicial numa das verificacoes, feitas a cada 100 passos.  Energia nao finita tambem conta como
//  instabilidade ( a versao usada no artigo nao fazia esse teste; ver README ).
//
//  Eixos do mapa ( os mesmos da Fig. 1 ):
//
//      i = 0 ... n_ma - 1 :   U_0 = U_min + ( i + 1 ) ( c_s - U_min ) / n_ma ,   U_min = 0.05 c_s
//      j = 0 ... n_nu - 1 :   log10 nu = -5 + 3 j / n_nu ,                        tau = 3 nu + 1/2
//
//  Uso:
//
//      ./stability_map_taylor_green --collision erlbm --tau-nh 0.65 --out map_erlbm.vtk
//      ./stability_map_taylor_green --collision reg                 --out map_reg.vtk
//
//      opcoes:  --L 16   --n-ma 400   --n-nu 400   --threads 0 ( todos )
//
//  Saida:  mapa VTK ( STRUCTURED_POINTS, n_ma x n_nu, i mais rapido ):  0 = estavel,  1 = instavel.
//
//  Paralelismo: cada ponto do mapa e uma simulacao independente, com populacoes proprias; as threads
//  dividem os pontos do mapa.  As funcoes da biblioteca chamadas dentro de cada simulacao rodam em serie.
//
//====================================================================================================================//

#define nvel 	19
#define dim 	3

#include "erlbm_common.h"

#include <vector>


//------ Uma simulacao: devolve 1 se instavel, 0 se estavel ----------------------------------------------------------//

int roda_ponto ( const GEOMETRY& geometry, LATTICE lattice, const PARAMETERS& parameters, int col, double u_0 )
{
	aloca_populacoes ( geometry, lattice );

	const double L = ( double ) geometry.nx;

	const double rho_ini = parameters.rho_ini;

	for ( int y = 0; y < geometry.ny; y++ )
	{
		for ( int x = 0; x < geometry.nx; x++ )
		{
			const int *meio = geometry.ini + x + y * geometry.nx;

			const double u_x =  u_0 * sin( 2. * M_PI * ( double ) x / L ) * cos( 2. * M_PI * ( double ) y / L );

			const double u_y = -u_0 * cos( 2. * M_PI * ( double ) x / L ) * sin( 2. * M_PI * ( double ) y / L );

			dist_eq ( lattice.inif + ( *meio - 1 ) * nvel, u_x, u_y, 0.0, rho_ini, lattice );
		}
	}

	const int n_steps = ( int ) ( ( double ) ( 100 * geometry.ny ) / u_0 );		// 100 T

	const double energy_ini = kinectic_energy ( geometry, lattice );

	int instavel = 0;

	for ( int step = 0; step < n_steps; step++ )
	{
		if ( step % 100 == 0 && step > 0 )
		{
			const double energy = kinectic_energy ( geometry, lattice );

			if ( energy > 2. * energy_ini || ! std::isfinite( energy ) )
			{
				instavel = 1;

				break;
			}
		}

		for ( int pto = 0; pto < geometry.fluid; pto++ )
		{
			double *f = lattice.inif + pto * nvel;

			colisao_sitio ( f, 0., 0., 0., lattice, parameters, col );

			double mx, my, mz;

			propag_site ( lattice, pto, mx, my, mz );
		}

		double *temp = lattice.inif;

		lattice.inif = lattice.inif_new;

		lattice.inif_new = temp;
	}

	libera_populacoes ( lattice );

	return instavel;
}


int main ( int argc, char* argv[] )
{
	ARGUMENTOS args { argc, argv };

	const int    col    = le_colisao ( args.arg_str ( "--collision", "erlbm" ) );
	const double tau_nh = args.arg_double ( "--tau-nh", 0.65 );
	const int    Li     = args.arg_int ( "--L", 16 );
	const int    n_ma   = args.arg_int ( "--n-ma", 400 );
	const int    n_nu   = args.arg_int ( "--n-nu", 400 );
	const string saida  = args.arg_str ( "--out", "stability_map_taylor_green.vtk" );

	define_threads ( args.arg_int ( "--threads", 0 ) );

	omp_set_max_active_levels ( 1 );

	//------ Geometria e rede ( compartilhadas ) --------------------------------------------------------------------//

	GEOMETRY geometry;

	geometria_periodica ( geometry, Li, Li, 1 );

	LATTICE lattice;

	define_rede ( geometry, lattice );

	//------ Eixos do mapa ------------------------------------------------------------------------------------------//

	const double c_s   = 1. / sqrt( 3. );
	const double u_min = 0.05 * c_s;

	const double log_nu_min = -5.;
	const double log_nu_max = -2.;

	cout << "\nMapa de estabilidade - Taylor-Green ( D3Q19, L = " << Li << " )" << endl;
	cout << "   colisao = " << nome_colisao ( col );
	if ( col == COL_EXT_REG ) cout << "  ( tau_nh = " << tau_nh << " )";
	cout << "\n   mapa    = " << n_ma << " x " << n_nu << "   ( Ma x log10 nu )" << endl;
	cout << "   threads = " << omp_get_max_threads() << endl;

	vector<int> mapa ( ( size_t ) n_ma * n_nu, 0 );

	int feitos = 0;

	#pragma omp parallel for collapse( 2 ) schedule( dynamic )
	for ( int j = 0; j < n_nu; j++ )
	{
		for ( int i = 0; i < n_ma; i++ )
		{
			const double u_0 = u_min + ( i + 1 ) * ( c_s - u_min ) / n_ma;

			const double log_nu = log_nu_min + j * ( log_nu_max - log_nu_min ) / n_nu;

			PARAMETERS parameters;

			parameters.tau     = 3. * pow( 10., log_nu ) + 0.5;
			parameters.visc    = pow( 10., log_nu );
			parameters.tau_nh  = ( col == COL_REG ) ? 1.0 : tau_nh;
			parameters.rho_ini = 1.0;

			mapa[ i + ( size_t ) j * n_ma ] = roda_ponto ( geometry, lattice, parameters, col, u_0 );

			#pragma omp atomic
			feitos++;

			if ( feitos % n_ma == 0 )
			{
				#pragma omp critical
				cout << "   " << feitos << " / " << n_ma * n_nu << " pontos" << endl;
			}
		}
	}

	//------ Gravacao -----------------------------------------------------------------------------------------------//

	ofstream fmap ( saida );

	fmap << "# vtk DataFile Version 2.0" << endl;
	fmap << "Stability map Taylor-Green: " << nome_colisao ( col ) << " tau_nh = " << tau_nh
	     << " ; x: U0 = Umin + (i+1)(cs-Umin)/" << n_ma << ", Umin = 0.05 cs ; y: log10 nu = -5 + 3 j/" << n_nu << endl;
	fmap << "ASCII" << endl;
	fmap << "DATASET STRUCTURED_POINTS" << endl;
	fmap << "DIMENSIONS " << n_ma << " " << n_nu << " 1" << endl;
	fmap << "ASPECT_RATIO 1 1 1" << endl;
	fmap << "ORIGIN 0 0 0" << endl;
	fmap << "POINT_DATA " << n_ma * n_nu << endl;
	fmap << "SCALARS Stability int" << endl;
	fmap << "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < n_nu; j++ )
	{
		for ( int i = 0; i < n_ma; i++ ) fmap << mapa[ i + ( size_t ) j * n_ma ] << " ";

		fmap << endl;
	}

	int n_inst = 0;

	for ( int v : mapa ) n_inst += v;

	cout << "\n" << saida << " gravado: " << n_inst << " de " << n_ma * n_nu << " pontos instaveis." << endl;

	return 0;
}

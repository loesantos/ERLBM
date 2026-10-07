//=============================== Mapa de estabilidade - escoamento em torno de um obstaculo ( Fig. 2 ) ==============//
//
//  Canal 20 x 12 ( x 1 ) com um bloco 2 x 4, rede D3Q19, Ma = 0.1.  Varia tau_nh ( eixo x do mapa ) e a
//  viscosidade ( eixo y ), independentemente.
//
//  Geometria ( nx = 20 + 2 colunas de fronteira, ny = 12, nz = 1 ):
//
//      x = 0 e x = nx - 1     colunas solidas ( fecham o dominio em x; as condicoes abaixo sobrescrevem
//                             as colunas fluidas adjacentes )
//      x = 1                  entrada: Zou-He com velocidade imposta  u = ( U_0, 0, 0 )
//      x = nx - 2             saida: populacoes copiadas do vizinho a montante e reescaladas para rho_ini
//      bloco                  x = 7 e 8,  y = 4 ... 7  ( centrado em y ), bounce-back a meio caminho
//      y, z                   periodicas
//
//  Fluido inicialmente em repouso, populacoes no equilibrio com rho_ini = 1.  Cada simulacao corre por
//  200 ny / U_0 passos e e classificada como instavel se, em algum passo, a densidade de algum sitio
//  ultrapassar 10 rho_ini ou nao for finita.
//
//  Eixos do mapa ( os mesmos da Fig. 2 ):
//
//      i = 0 ... n_tau - 1 :   tau_nh   = 0.5 + 0.6 i / n_tau
//      j = 0 ... n_nu  - 1 :   log10 nu = -5 + 4.3 j / n_nu ,       tau = 3 nu + 1/2
//
//  Uso:
//
//      ./stability_map_obstacle --out map_obstacle.vtk        opcoes:  --n-tau 400  --n-nu 400  --ma 0.1  --threads 0
//
//  Saida:  mapa VTK ( n_tau x n_nu, i mais rapido ):  0 = estavel,  1 = instavel.
//
//====================================================================================================================//

#define nvel 	19
#define dim 	3

#include "erlbm_common.h"

#include <vector>


//------ Geometria do canal com o bloco -----------------------------------------------------------------------------//

void geometria_obstaculo ( GEOMETRY& geometry )
{
	const int nx_canal = 20;
	const int ny       = 12;

	geometry.file  = "";
	geometry.ftesc = 1.0;
	geometry.nx    = nx_canal + 2;
	geometry.ny    = ny;
	geometry.nz    = 1;

	geometry.ini = new int[ geometry.nx * geometry.ny ];

	for ( int y = 0; y < geometry.ny; y++ )
	{
		for ( int x = 0; x < geometry.nx; x++ )
		{
			int s = 1;

			if ( x == 0 || x == geometry.nx - 1 ) s = 0;					// colunas de fechamento

			if ( ( x == 7 || x == 8 ) && y >= 4 && y <= 7 ) s = 0;			// bloco 2 x 4

			geometry.ini[ x + y * geometry.nx ] = s;
		}
	}

	numera_fluido ( geometry );
}


//------ Uma simulacao: devolve 1 se instavel, 0 se estavel ----------------------------------------------------------//

int roda_ponto ( const GEOMETRY& geometry, LATTICE lattice, const PARAMETERS& parameters, double u_imp )
{
	aloca_populacoes ( geometry, lattice );

	const double rho_ini = parameters.rho_ini;

	for ( int pto = 0; pto < geometry.fluid; pto++ )
		dist_eq ( lattice.inif + pto * nvel, 0., 0., 0., rho_ini, lattice );

	const int n_steps = ( int ) ( ( double ) ( 200 * geometry.ny ) / u_imp );

	int instavel = 0;

	for ( int step = 0; step < n_steps && ! instavel; step++ )
	{
		//------ Condicoes de contorno de entrada e saida ------------------------------------------------//

		ZouHeD3Q19_vx_x ( 1, u_imp, geometry, lattice );

		devnull_rho_x ( geometry.nx - 2, rho_ini, geometry, lattice );

		//------ Colisao, verificacao e propagacao -------------------------------------------------------//

		for ( int pto = 0; pto < geometry.fluid; pto++ )
		{
			double *f = lattice.inif + pto * nvel;

			coll_Ext_Reg ( f, 0., 0., 0., lattice, parameters );

			const double dens = density ( f );

			if ( fabs ( dens ) > 10. * rho_ini || ! std::isfinite( dens ) ) instavel = 1;

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

	const int    n_tau = args.arg_int ( "--n-tau", 400 );
	const int    n_nu  = args.arg_int ( "--n-nu", 400 );
	const double Ma    = args.arg_double ( "--ma", 0.1 );
	const string saida = args.arg_str ( "--out", "stability_map_obstacle.vtk" );

	define_threads ( args.arg_int ( "--threads", 0 ) );

	omp_set_max_active_levels ( 1 );		// lacos paralelos internos da biblioteca ( contorno ) em serie

	GEOMETRY geometry;

	geometria_obstaculo ( geometry );

	LATTICE lattice;

	define_rede ( geometry, lattice );

	const double u_imp = Ma / sqrt( 3. );

	const double tau_nh_min = 0.5, tau_nh_max = 1.1;
	const double log_nu_min = -5., log_nu_max = -0.7;

	cout << "\nMapa de estabilidade - obstaculo ( D3Q19, " << geometry.nx - 2 << " x " << geometry.ny
	     << ", Ma = " << Ma << " ), E-RLBM" << endl;
	cout << "   mapa    = " << n_tau << " x " << n_nu << "   ( tau_nh x log10 nu )" << endl;
	cout << "   threads = " << omp_get_max_threads() << endl;

	vector<int> mapa ( ( size_t ) n_tau * n_nu, 0 );

	int feitos = 0;

	#pragma omp parallel for collapse( 2 ) schedule( dynamic )
	for ( int j = 0; j < n_nu; j++ )
	{
		for ( int i = 0; i < n_tau; i++ )
		{
			const double log_nu = log_nu_min + j * ( log_nu_max - log_nu_min ) / n_nu;

			PARAMETERS parameters;

			parameters.tau     = 3. * pow( 10., log_nu ) + 0.5;
			parameters.visc    = pow( 10., log_nu );
			parameters.tau_nh  = tau_nh_min + ( double ) i * ( tau_nh_max - tau_nh_min ) / n_tau;
			parameters.rho_ini = 1.0;

			mapa[ i + ( size_t ) j * n_tau ] = roda_ponto ( geometry, lattice, parameters, u_imp );

			#pragma omp atomic
			feitos++;

			if ( feitos % n_tau == 0 )
			{
				#pragma omp critical
				cout << "   " << feitos << " / " << n_tau * n_nu << " pontos" << endl;
			}
		}
	}

	ofstream fmap ( saida );

	fmap << "# vtk DataFile Version 2.0" << endl;
	fmap << "Stability map obstacle: E-RLBM, Ma = " << Ma << " ; x: tau_nh = 0.5 + 0.6 i/" << n_tau
	     << " ; y: log10 nu = -5 + 4.3 j/" << n_nu << endl;
	fmap << "ASCII" << endl;
	fmap << "DATASET STRUCTURED_POINTS" << endl;
	fmap << "DIMENSIONS " << n_tau << " " << n_nu << " 1" << endl;
	fmap << "ASPECT_RATIO 1 1 1" << endl;
	fmap << "ORIGIN 0 0 0" << endl;
	fmap << "POINT_DATA " << n_tau * n_nu << endl;
	fmap << "SCALARS Stability int" << endl;
	fmap << "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < n_nu; j++ )
	{
		for ( int i = 0; i < n_tau; i++ ) fmap << mapa[ i + ( size_t ) j * n_tau ] << " ";

		fmap << endl;
	}

	int n_inst = 0;

	for ( int v : mapa ) n_inst += v;

	cout << "\n" << saida << " gravado: " << n_inst << " de " << n_tau * n_nu << " pontos instaveis." << endl;

	return 0;
}

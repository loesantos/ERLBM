//=============================== Vortice de Taylor-Green - precisao ( Fig. 3 e verificacao de kappa_4 ) =============//
//
//  Dominio periodico L x L, rede D2Q9 ( ou D3Q19 com nz = 1, compilando com -DLATT=19 ).
//
//      u_x =  U_0 sin( k x ) cos( k y ) e^{-2 nu k^2 t},   u_y = - U_0 cos( k x ) sin( k y ) e^{-2 nu k^2 t},   k = 2 pi / L
//
//  Densidade inicial ( --density poisson, padrao ):  solucao exata da equacao de Poisson da pressao,
//
//      rho = 1 + ( 3 U_0^2 / 4 ) ( cos 2kx + cos 2ky )
//
//  Populacoes iniciais ( --init neq, padrao ):  f = f_eq + f^(1),   f^(1) = w / ( 2 cs^4 ) H : Pi^(1),
//  Pi^(1) = - rho cs^2 tau ( grad u + grad u^T ).  Para o Taylor-Green Sxy = 0 e Sxx = - Syy = 2 U_0 k cos kx cos ky.
//  Com --init eq as populacoes sao so f_eq ( deixa um erro de amplitude ~ tau ( 1 - tau ) K^2 / 3, K^2 = 2 k^2 ).
//
//  Dois modos de uso:
//
//  1) Fig. 3 ( escala difusiva, padrao ):  Re fixo, U_0 = 2 / L ( ou --u0 ), tau = 1/2 + 3 U_0 L / Re ,
//     erro medido em t = T = L / U_0 .
//
//          ./taylor_green --L 64 --collision erlbm --tau-nh 0.65              ( Re = 50 )
//
//  2) Verificacao do coeficiente hiperviscoso ( --tau dado ):  erro medido quando a amplitude cai a e^{-1},
//     isto e, em 2 nu k^2 t = 1 .
//
//          ./taylor_green --L 64 --tau 0.8 --u0 0.001 --collision erlbm --tau-nh 0.65
//
//  Opcoes:  --collision bgk | reg | erlbm   --tau-nh   --Re ( 50 )   --u0   --tau   --init neq | eq
//           --density poisson | uniform      --threads
//
//  Saida ( uma linha ):   L  tau  tau_nh  U_0  passos  E_vetorial  E_modulo  ( nu_ef - nu ) / nu
//
//      E_vetorial = sqrt( sum |u - u_ex|^2 / sum |u_ex|^2 )       ( Eq. (l2) do artigo )
//      E_modulo   = sqrt( sum ( |u| - |u_ex| )^2 / sum |u_ex|^2 )
//      nu_ef      : viscosidade que reproduz a amplitude do modo de Taylor-Green no instante medido
//
//====================================================================================================================//

#ifndef LATT
#define LATT 9
#endif

#if LATT == 9
#define nvel 	9
#define dim 	2
#else
#define nvel 	19
#define dim 	3
#endif

#include "erlbm_common.h"

#include <cstdio>


//------ Populacoes iniciais ----------------------------------------------------------------------------------------//

void inicializa_taylor_green ( const GEOMETRY& geometry, LATTICE lattice, double u_0, double tau,
                               bool poisson, bool neq )
{
	const int Li = geometry.nx;

	const double L = ( double ) Li;

	const double k = 2. * M_PI / L;

	const double cs2 = 1. / 3.;

	for ( int y = 0; y < Li; y++ )
	{
		for ( int x = 0; x < Li; x++ )
		{
			double rho = 1.0;

			if ( poisson ) rho = 1.0 + 0.75 * u_0 * u_0 * ( cos( 4. * M_PI * x / L ) + cos( 4. * M_PI * y / L ) );

			const double u_x =  u_0 * sin( 2. * M_PI * x / L ) * cos( 2. * M_PI * y / L );
			const double u_y = -u_0 * cos( 2. * M_PI * x / L ) * sin( 2. * M_PI * y / L );

			double *f = lattice.inif + ( geometry.ini[ x + y * Li ] - 1 ) * nvel;

			dist_eq ( f, u_x, u_y, 0.0, rho, lattice );

			if ( ! neq ) continue;

			const double Sxx = 2. * u_0 * k * cos( k * x ) * cos( k * y );

			const double Pxx = - rho * cs2 * tau * Sxx;

			const double Pyy = - Pxx;

			for ( int i = 0; i < nvel; i++ )
			{
				const double cx = lattice.c_i[ i * dim ];
				const double cy = lattice.c_i[ i * dim + 1 ];

				f[i] += lattice.w[i] / ( 2. * cs2 * cs2 ) * ( ( cx * cx - cs2 ) * Pxx + ( cy * cy - cs2 ) * Pyy );
			}
		}
	}
}


int main ( int argc, char* argv[] )
{
	ARGUMENTOS args { argc, argv };

	const int    Li       = args.arg_int ( "--L", 64 );
	const int    col      = le_colisao ( args.arg_str ( "--collision", "erlbm" ) );
	const double tau_nh   = args.arg_double ( "--tau-nh", 0.65 );
	const double Re       = args.arg_double ( "--Re", 50. );
	const double u_0      = args.arg_double ( "--u0", 2. / Li );
	const bool   modo_tau = args.tem ( "--tau" );
	const bool   neq      = args.arg_str ( "--init", "neq" ) == "neq";
	const bool   poisson  = args.arg_str ( "--density", "poisson" ) == "poisson";

	define_threads ( args.arg_int ( "--threads", 0 ) );

	//------ Parametros ---------------------------------------------------------------------------------------------//

	PARAMETERS p;

	if ( modo_tau )
	{
		p.tau  = args.arg_double ( "--tau", 0.6 );
		p.visc = ( p.tau - 0.5 ) / 3.;
	}
	else
	{
		p.visc = u_0 * Li / Re;
		p.tau  = 0.5 + 3. * p.visc;
	}

	p.tau_nh  = ( col == COL_REG ) ? 1.0 : ( col == COL_BGK ? p.tau : tau_nh );
	p.rho_ini = 1.0;

	//------ Geometria, rede e condicao inicial ---------------------------------------------------------------------//

	GEOMETRY geometry;

	geometria_periodica ( geometry, Li, Li, 1 );

	LATTICE lattice;

	define_rede ( geometry, lattice );

	aloca_populacoes ( geometry, lattice );

	inicializa_taylor_green ( geometry, lattice, u_0, p.tau, poisson, neq );

	//------ Numero de passos ---------------------------------------------------------------------------------------//

	const double L = Li;

	const double k = 2. * M_PI / L;

	const int nT = modo_tau ? ( int ) round( 1. / ( 2. * p.visc * k * k ) )		// 2 nu k^2 t = 1
	                        : ( int ) round( L / u_0 );								// t = T

	//------ Laco principal -----------------------------------------------------------------------------------------//
	//
	//  A medida e feita depois da colisao do passo nT e antes da troca de ponteiros: lattice.inif guarda as
	//  populacoes pos-colisao, cujos momentos sao os do instante nT ( a colisao os conserva ).

	for ( int step = 0; step <= nT; step++ )
	{
		#pragma omp parallel for
		for ( int pto = 0; pto < geometry.fluid; pto++ )
		{
			double *f = lattice.inif + pto * nvel;

			colisao_sitio ( f, 0., 0., 0., lattice, p, col );

			double mx, my, mz;

			propag_site ( lattice, pto, mx, my, mz );
		}

		if ( step == nT )
		{
			double s_mag = 0, s_vec = 0, s_r = 0, s_proj = 0, s_phi = 0;

			const double e = exp( -2. * p.visc * k * k * step );

			for ( int y = 0; y < Li; y++ )
			{
				for ( int x = 0; x < Li; x++ )
				{
					double vx, vy, vz, rho;

					calcula ( lattice.inif + ( x + y * Li ) * nvel, vx, vy, vz, rho, lattice );

					const double rx =  u_0 * e * sin( k * x ) * cos( k * y );
					const double ry = -u_0 * e * cos( k * x ) * sin( k * y );

					const double us = sqrt( vx * vx + vy * vy + vz * vz ), ur = sqrt( rx * rx + ry * ry );

					s_mag += ( us - ur ) * ( us - ur );
					s_vec += ( vx - rx ) * ( vx - rx ) + ( vy - ry ) * ( vy - ry ) + vz * vz;
					s_r   += ur * ur;

					const double px = sin( k * x ) * cos( k * y ), py = -cos( k * x ) * sin( k * y );

					s_proj += vx * px + vy * py;
					s_phi  += px * px + py * py;
				}
			}

			const double A = s_proj / s_phi / u_0;

			const double nu_ef = - log( A ) / ( 2. * k * k * step );

			printf ( "%d %.6g %.6g %.6g %d %.6g %.6g %.6g\n", Li, p.tau, p.tau_nh, u_0, step,
			         sqrt( s_vec / s_r ), sqrt( s_mag / s_r ), nu_ef / p.visc - 1. );
		}

		double *t = lattice.inif; lattice.inif = lattice.inif_new; lattice.inif_new = t;
	}

	return 0;
}

//=============================== Camada de cisalhamento dupla ( Figs. 4 e 5 ) =======================================//
//
//  Camada de cisalhamento dupla periodica ( Minion & Brown 1997 ), rede D2Q9, dominio L x L periodico.
//
//          u_x = U_0 tanh( lambda ( y/L - 1/4 ) )        se  y <= L / 2
//          u_x = U_0 tanh( lambda ( 3/4 - y/L ) )        se  y >  L / 2
//          u_y = epsilon U_0 sin( 2 pi ( x/L + 1/4 ) )
//
//  com U_0 = 1/32, lambda = 80, epsilon = 0.05.  Densidade inicial uniforme, populacoes no equilibrio.
//  Tempo adimensional t* = passo U_0 / L ( t* = 1 e T = L / U_0 ).  Re = U_0 L / nu,  tau = 1/2 + 3 nu.
//
//  Uso:
//
//      ./double_shear_layer --L 256 --collision erlbm --tau-nh 0.65          ( Re = 80000 )
//
//  Opcoes:  --collision bgk | mrt | reg | erlbm   --tau-nh   --Re ( 80000 )   --tau ( em vez de Re )
//           --t-final ( 1.5 )   --dt-out ( 1/16 )   --vorticity 0 | 1 ( 1 )   --threads ( 0 = todos )
//
//  Saidas ( na pasta corrente ):
//
//      k_energy_lbm.dat      t*  e energia cinetica total ( sum 0.5 |u|^2 ), a cada dt-out
//      vel_<t*>.vtk          campo de velocidade a cada dt-out ( vel_1.000000.vtk e o usado no erro da Fig. 4 )
//      vor_<passo>.vtk       vorticidade a cada dt-out ( Fig. 5 )
//
//  O erro em relacao a uma referencia e calculado por tools/error_calculator.
//
//  Compilacao com OpenACC ( GPU ):  nvc++ -std=c++17 -O3 -acc=gpu -gpu=ccnative -Minfo=accel ...
//
//====================================================================================================================//

#define nvel 	9
#define dim 	2

#include "erlbm_common.h"

#include <chrono>


const double U_0     = 1. / 32.;		// velocidade de referencia ( unidades de rede )
const double LAMBDA  = 80.;				// espessura das camadas de cisalhamento ( ~ 1 / lambda )
const double EPSILON = 0.05;			// amplitude da perturbacao em u_y

const double U_LIMITE = 0.4;			// acima disto ( velocidade rms ) a simulacao e dada por instavel e para


//------ Campo inicial ----------------------------------------------------------------------------------------------//

void initial_double_shear_layer ( GEOMETRY geometry, LATTICE lattice, double L, double rho_ini )
{
	for ( int y = 0; y < geometry.ny; y++ )
	{
		for ( int x = 0; x < geometry.nx; x++ )
		{
			const int *meio = geometry.ini + x + y * geometry.nx;

			double u_x = U_0 * tanh( LAMBDA * ( ( double ) y / L - 1. / 4. ) );

			if ( y > geometry.ny / 2 )
				u_x = U_0 * tanh( LAMBDA * ( 3. / 4. - ( double ) y / L ) );

			const double u_y = EPSILON * U_0 * sin( 2. * M_PI * ( ( double ) x / L + 1. / 4. ) );

			dist_eq ( lattice.inif + ( *meio - 1 ) * nvel, u_x, u_y, 0.0, rho_ini, lattice );
		}
	}
}


int main ( int argc, char* argv[] )
{
	ARGUMENTOS args { argc, argv };

	const int    Li       = args.arg_int ( "--L", 256 );
	const int    col      = le_colisao ( args.arg_str ( "--collision", "erlbm" ) );
	const double tau_nh   = args.arg_double ( "--tau-nh", 0.65 );
	const double Re_in    = args.arg_double ( "--Re", 80000. );
	const double T_FINAL  = args.arg_double ( "--t-final", 1.5 );
	const double DT_SAIDA = args.arg_double ( "--dt-out", 1. / 16. );
	const bool   grava_vor = args.arg_int ( "--vorticity", 1 ) != 0;

	define_threads ( args.arg_int ( "--threads", 0 ) );

	const double L = ( double ) Li;

	PARAMETERS parameters;

	parameters.tau     = args.tem ( "--tau" ) ? args.arg_double ( "--tau", 0.5 ) : 0.5 + 3. * U_0 * L / Re_in;
	parameters.visc    = ( parameters.tau - 0.5 ) / 3.;
	parameters.tau_nh  = ( col == COL_REG ) ? 1.0 : tau_nh;
	parameters.rho_ini = 1.0;

	const double Re = U_0 * L / parameters.visc;

#ifdef _OPENACC
	const bool on_device = ( acc_get_device_type() != acc_device_host );

	cout << "\nOpenACC: " << ( on_device ? "executando no dispositivo" : "executando no hospedeiro" ) << endl;
#else
	cout << "\nBackend: CPU / OpenMP ( " << omp_get_max_threads() << " threads )" << endl;
#endif

	//------ Geometria e rede ---------------------------------------------------------------------------------------//

	GEOMETRY geometry;

	geometria_periodica ( geometry, Li, Li, 1 );

	LATTICE lattice;

	define_rede ( geometry, lattice );

	aloca_populacoes ( geometry, lattice );

	const int n_fluid = geometry.fluid;

	const int size_to_alloc = n_fluid * nvel;

	//------ Passos e gravacao ---------------------------------------------------------------------------------------//

	const int passos_por_t = ( int ) round( L / U_0 );						// passos por unidade de t*

	const int n_steps = ( int ) round( T_FINAL * passos_por_t ) + 1;		// inclui t* = T_FINAL

	const int interval = max( 1, ( int ) round( DT_SAIDA * passos_por_t ) );

	cout << "\nCaso: double shear layer ( D2Q9 )" << endl;
	cout << "   colisao              = " << nome_colisao ( col ) << endl;
	cout << "   L                    = " << L << endl;
	cout << "   tau / tau_nh         = " << parameters.tau << " / " << parameters.tau_nh << endl;
	cout << "   Re = U_0 L / nu      = " << Re << endl;
	cout << "   passos               = " << n_steps << "  ( t* final = " << T_FINAL << " )" << endl;
	cout << "   gravacao a cada      = " << interval << " passos" << endl;

	//------ Condicao inicial ---------------------------------------------------------------------------------------//

	initial_double_shear_layer ( geometry, lattice, L, parameters.rho_ini );

#ifdef _OPENACC

	#pragma acc enter data copyin( lattice )

	#pragma acc enter data copyin( lattice.c_i[ 0 : nvel * dim ] )
	#pragma acc enter data copyin( lattice.Q_i[ 0 : nvel * dim * dim ] )

	#pragma acc enter data copyin( lattice.ini_stream[ 0 : size_to_alloc ] )

	#pragma acc enter data copyin( lattice.inif    [ 0 : size_to_alloc ] )
	#pragma acc enter data create( lattice.inif_new[ 0 : size_to_alloc ] )

#endif

	ofstream file_energy ( "k_energy_lbm.dat" );

	file_energy << setprecision ( 10 );

	const auto time_start = chrono::steady_clock::now();

	int step = 0;

	//============================== Laco principal ==================================================================//

	for ( step = 0; step < n_steps; step++ )
	{
		const double t_adm = ( double ) step * U_0 / L;

#ifdef _OPENACC
		#pragma acc parallel loop gang vector present( lattice )
#else
		#pragma omp parallel for
#endif
		for ( int pto = 0; pto < n_fluid; pto++ )
		{
			double *f = lattice.inif + pto * nvel;

			colisao_sitio ( f, 0., 0., 0., lattice, parameters, col );

			double mx, my, mz;

			propag_site ( lattice, pto, mx, my, mz );
		}

		//------ Energia e campos -------------------------------------------------------------------------------//
		//
		//  Gravados depois da colisao e antes da troca de ponteiros: lattice.inif guarda as populacoes
		//  pos-colisao do passo "step", cujos momentos ( rho, rho u ) sao os do proprio passo.

		if ( step % interval == 0 )
		{
#ifdef _OPENACC
			#pragma acc update self( lattice.inif[ 0 : size_to_alloc ] )
#endif
			const double k_energy = kinectic_energy ( geometry, lattice );

			cout << "\nt* = " << t_adm << "	energia cinetica = " << k_energy << endl;

			file_energy << t_adm << " " << k_energy << endl;

			if ( ! std::isfinite( k_energy ) || sqrt ( 2. * k_energy / n_fluid ) > U_LIMITE )
			{
				cout << "\nSimulacao instavel no passo " << step << " ( t* = " << t_adm << " )." << endl;

				break;
			}

			rec_velocity ( geometry, lattice, t_adm );

			if ( grava_vor ) rec_vorticity ( geometry, lattice, ( unsigned int ) step );
		}

		double *temp = lattice.inif;

		lattice.inif = lattice.inif_new;

		lattice.inif_new = temp;

#ifdef _OPENACC
		if ( on_device )
		{
			#pragma acc serial present( lattice )
			{
				double *t = lattice.inif;

				lattice.inif = lattice.inif_new;

				lattice.inif_new = t;
			}
		}
#endif
	}

	file_energy.close();

	const double elapsed = chrono::duration<double> ( chrono::steady_clock::now() - time_start ).count();

	cout << "\nPassos executados = " << step << endl;
	cout << "Tempo do laco principal = " << elapsed << " s" << endl;
	cout << "Desempenho medio = " << ( elapsed > 0 ? ( double ) n_fluid * step / ( elapsed * 1.0e6 ) : 0. ) << " MLUPS" << endl;

#ifdef _OPENACC
	#pragma acc exit data delete( lattice.inif    [ 0 : size_to_alloc ] )
	#pragma acc exit data delete( lattice.inif_new[ 0 : size_to_alloc ] )
	#pragma acc exit data delete( lattice.ini_stream[ 0 : size_to_alloc ] )
	#pragma acc exit data delete( lattice.Q_i[ 0 : nvel * dim * dim ] )
	#pragma acc exit data delete( lattice.c_i[ 0 : nvel * dim ] )
	#pragma acc exit data delete( lattice )
#endif

	return 0;
}

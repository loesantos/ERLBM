//=============================== Erro L2 relativo do campo de velocidade ============================================//
//
//  Compara o campo de velocidade de uma simulacao ( malha L x L ) com o de referencia ( malha L_ref x L_ref ),
//  amostrando a referencia nos mesmos pontos fisicos:  ( x, y ) na malha L  <->  ( s x, s y ) na de referencia,
//  com s = L_ref / L ( inteiro ).  Sao calculadas duas normas:
//
//      vetorial ( Eq. do artigo ):   E_vec = sqrt(  sum | u_s - u_r |^2      /  sum |u_r|^2  )
//      do modulo ( versao 2025 ):    E_mod = sqrt(  sum ( |u_s| - |u_r| )^2  /  sum |u_r|^2  )
//
//  A versao anterior deste programa calculava apenas E_mod, que e a usada na Fig. 4 original.
//  Pela desigualdade triangular E_mod <= E_vec.
//
//  Uso:
//          ./bin/error_calculator  [ vel_simulacao.vtk ]  [ vel_referencia.vtk ]
//
//  Sem argumentos: vel_1.000000.vtk e vel_ref.vtk ( como na versao anterior ).
//  Os dois arquivos devem estar em VTK legacy ASCII ( rec_velocity da SimBoltz ).
//
//====================================================================================================================//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

bool read_velocity ( const string&, vector<double>&, int&, int&, int& );


int main ( int argc, char* argv[] )
{
	string name_vel_cal = "vel_1.000000.vtk";

	string name_vel_ref = "vel_ref.vtk";

	if ( argc > 1 ) name_vel_cal = argv[1];
	if ( argc > 2 ) name_vel_ref = argv[2];

	//----------------------------------------------------------------------------------------------------------------//

	vector<double> vel_cal, vel_ref;

	int nx, ny, nz, nx_ref, ny_ref, nz_ref;

	if ( ! read_velocity ( name_vel_cal, vel_cal, nx, ny, nz ) ) return 1;

	if ( ! read_velocity ( name_vel_ref, vel_ref, nx_ref, ny_ref, nz_ref ) ) return 1;

	if ( nz != 1 || nz_ref != 1 || nx_ref % nx != 0 || ny_ref % ny != 0 || nx_ref / nx != ny_ref / ny )
	{
		cout << "\nAs malhas nao sao compativeis: a referencia tem de ser um refinamento inteiro da simulacao."
		     << endl;

		return 1;
	}

	const int amost = nx_ref / nx;

	//----------------------------------------------------------------------------------------------------------------//

	double sum_dif2 = 0.0;		// ( |u_s| - |u_r| )^2

	double sum_vec2 = 0.0;		// | u_s - u_r |^2

	double sum_ur2 = 0.0;

	for ( int y = 0; y < ny; y++ )
	{
		for ( int x = 0; x < nx; x++ )
		{
			const size_t pos_cal = 3 * ( ( size_t ) x + ( size_t ) y * nx );

			const size_t pos_ref = 3 * ( ( size_t ) x * amost + ( size_t ) y * amost * nx_ref );

			const double vx_ref = vel_ref[ pos_ref     ];
			const double vy_ref = vel_ref[ pos_ref + 1 ];

			const double u_r = sqrt( vx_ref * vx_ref + vy_ref * vy_ref );

			const double vx_cal = vel_cal[ pos_cal     ];
			const double vy_cal = vel_cal[ pos_cal + 1 ];

			const double u_s = sqrt( vx_cal * vx_cal + vy_cal * vy_cal );

			sum_dif2 = sum_dif2 + ( u_s - u_r ) * ( u_s - u_r );

			sum_vec2 = sum_vec2 + ( vx_cal - vx_ref ) * ( vx_cal - vx_ref ) + ( vy_cal - vy_ref ) * ( vy_cal - vy_ref );

			sum_ur2 = sum_ur2 + u_r * u_r;
		}
	}

	const double erro_mod = sqrt( sum_dif2 / sum_ur2 );

	const double erro_vec = sqrt( sum_vec2 / sum_ur2 );

	cout << "\nL = " << nx << "   L_ref = " << nx_ref << endl;

	cout << setprecision( 6 );

	cout << "\nErro (vetorial, Eq. do artigo)  = " << erro_vec << endl;

	cout <<   "Erro (modulo, versao 2025)      = " << erro_mod << endl;

	//  Uma linha para juntar numa tabela:  L  E_vec  E_mod

	ofstream ftab ( "erro_L2.dat" );

	ftab << nx << " " << setprecision( 10 ) << erro_vec << " " << erro_mod << endl;

	return 0;
}


//====================================================================================================================//

bool read_velocity ( const string& nome_vel, vector<double>& vel, int& nx, int& ny, int& nz )
{
	ifstream fmatriz( nome_vel );

	if ( ! fmatriz )
	{
		cout << "\nNao foi possivel abrir " << nome_vel << endl;

		return false;
	}

	string line, dump;

	stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );

	getline( fmatriz, line );

	dados << line;

	dados >> dump >> nx >> ny >> nz;

	cout << "\n" << nome_vel << ":  x = " << nx << ";  y = " << ny << ";  z = " << nz << endl;

	//  Pula ASPECT_RATIO, ORIGIN, POINT_DATA e VECTORS.

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );

	if ( dump.find( "VECTORS" ) == string::npos )
	{
		cout << "\nCabecalho inesperado em " << nome_vel << " (esperado VECTORS, em ASCII)." << endl;

		return false;
	}

	vel.assign( 3 * ( size_t ) nx * ny * nz, 0.0 );

	for ( size_t k = 0; k < vel.size(); k++ ) fmatriz >> vel[k];

	if ( ! fmatriz )
	{
		cout << "\n" << nome_vel << " terminou antes do esperado." << endl;

		return false;
	}

	return true;
}

//====================================================================================================================//

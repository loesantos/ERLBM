


//================================ Etapa de colisão usando BGK =======================================================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_BGK ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{
	double tau = parameters.tau;

	double f_eq[nvel];

	double op_col[nvel];
	
	double S[nvel];

	double vx, vy, vz, rho;

	double one_over_tau = 1.0 / tau;

	calcula ( f, vx, vy, vz, rho, lattice );

	double vx_alt = vx + 0.5 * acc_x;
	double vy_alt = vy + 0.5 * acc_y;
	double vz_alt = vz + 0.5 * acc_z;
	
	double Fx = acc_x * rho;    
	double Fy = acc_y * rho;
	double Fz = acc_z * rho;
	
	source ( Fx, Fy, Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );
	
	if ( nvel == 77 ) dist_eq_sixth ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
	
	else dist_eq ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
	
	for ( int i = 0; i < nvel; i++ )
	{
		op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;

		f[i] = f[i] + op_col[i] + S[i];
	}
	
}

//====================================================================================================================//




//================================ Etapa de colisão usando BGK for heat conduction ===================================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_BGK ( double *f,  LATTICE lattice, PARAMETERS parameters )
{
	double tau = parameters.tau;
	
	double one_over_tau = 1.0 / tau;

	double f_eq[nvel];

	double T = density ( f );
	
	dist_eq ( f_eq, T, lattice );	
	
	for ( int i = 0; i < nvel; i++ ) f[i] = f[i] + ( f_eq[i] - f[i] ) * one_over_tau;
}

//====================================================================================================================//



//================================ Etapa de colisão usando BGK =======================================================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_BGK_thermal ( double *f, double acc_x, double acc_y, double acc_z, double T_ref, LATTICE lattice, 
						PARAMETERS parameters )
{
	double tau = parameters.tau;
	
	double op_col[nvel];
	
	double S[nvel];

	double vx, vy, vz, rho, Tmp;
	
	double one_over_tau = 1.0 / tau;

	calcula_th ( f, vx, vy, vz, rho, Tmp, lattice );
	
	double theta = Tmp / T_ref - 1.0;
			
	double vx_alt = vx + 0.5 * acc_x;
	double vy_alt = vy + 0.5 * acc_y;
	double vz_alt = vz + 0.5 * acc_z;
	
	double Fx = acc_x * rho;    
	double Fy = acc_y * rho;
	double Fz = acc_z * rho;
	
	source ( Fx, Fy, Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );
	
	theta = 0; // ISOTÉRMICO

	double f_eq[nvel];
	
	dist_eq_thermal ( f_eq, vx_alt, vy_alt, vz_alt, rho, theta, lattice ); 
	
	for ( int i = 0; i < nvel; i++ )
	{
		op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;

		f[i] = f[i] + op_col[i] + S[i];
	}
}

//====================================================================================================================//




//================================ Etapa de colisão usando BGK =======================================================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_BGK_fourth ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{
	double tau = parameters.tau;
	
    double f_eq[nvel];

    double op_col[nvel];
    
    double S[nvel];

    double vx, vy, vz, rho;

    double one_over_tau = 1.0 / tau;

    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    dist_eq_fourth ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    for ( int i = 0; i < nvel; i++ )
    {
        op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;

        f[i] = f[i] + op_col[i] + S[i];
    }

}

//====================================================================================================================//



//=============================== Collision step using the TRT model =================================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_TRT ( double* f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{
	double tau_sim = parameters.tau;
	
	const double inv_tau = 1.0 / tau_sim;

    const double tau_ant = ( ( 8.0 - inv_tau ) / ( 8.0 * (  2.0 - inv_tau ) ) ); 
	
    double vx, vy, vz, rho;

    calcula ( f, vx, vy, vz, rho, lattice );

    double op_col_sim[nvel];

    op_bgk_even ( f, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau_sim, op_col_sim, lattice );

    double op_col_ant[nvel];

    op_bgk_odd ( f, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau_ant, op_col_ant, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f[i] + op_col_sim[i] + op_col_ant[i];

    }
}

//====================================================================================================================//



//=============================== Return the BGK even collision operator (TRT model) ================================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//====================================================================================================================//

#pragma acc routine seq
void op_bgk_even ( double f[nvel], double vx, double vy, double vz, double rho, double acc_x, double acc_y, double acc_z, 
				double tau, double op_col[nvel], LATTICE lattice )
{
    double f_eq[nvel];
    
    double S[nvel]{};

    double one_over_tau = 1.0 / tau;

    double vx_acc = vx + acc_x * 0.5;
    double vy_acc = vy + acc_y * 0.5;
    double vz_acc = vz + acc_z * 0.5;
    
    double Fx = acc_x * rho;
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;  

	if ( acc_x != 0 || acc_y != 0 || acc_z != 0 )
	{
		source ( Fx, Fy, Fz, vx_acc, vy_acc, vz_acc, rho, tau, S, lattice );  
	}

    dist_eq ( f_eq, vx_acc, vy_acc, vz_acc, rho, lattice );

    op_col[0] = ( f_eq[0] - f[0] ) * one_over_tau + S[0];

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_plus = ( f[i] + f[i+1] ) * 0.5;

        double f_eq_plus = ( f_eq[i] + f_eq[i+1] ) * 0.5;
        
        double S_plus = ( S[i] + S[i+1] ) * 0.5;

        op_col[i] = ( f_eq_plus - f_plus ) * one_over_tau + S_plus;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_plus = ( f[i] + f[i-1] ) * 0.5;

        double f_eq_plus = ( f_eq[i] + f_eq[i-1] ) * 0.5;
        
        double S_plus = ( S[i] + S[i-1] ) * 0.5;

        op_col[i] = ( f_eq_plus - f_plus ) * one_over_tau + S_plus;
    }

}

//====================================================================================================================//



//=============================== Calcula o operador BGK antisimétrico (modelo TRT) ==================================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//====================================================================================================================//

#pragma acc routine seq
void op_bgk_odd ( double f[nvel], double vx, double vy, double vz, double rho, double acc_x, double acc_y, double acc_z, 
				double tau, double op_col[nvel], LATTICE lattice )
{
    double f_eq[nvel];
    
    double S[nvel]{};

    double one_over_tau = 1.0 / tau;

    double vx_acc = vx + acc_x * 0.5;
    double vy_acc = vy + acc_y * 0.5;
    double vz_acc = vz + acc_z * 0.5;
    
    double Fx = acc_x * rho;
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;  

	if ( acc_x != 0 || acc_y != 0 || acc_z != 0 )
	{
		source ( Fx, Fy, Fz, vx_acc, vy_acc, vz_acc, rho, tau, S, lattice );  
	}
	
    dist_eq ( f_eq, vx_acc, vy_acc, vz_acc, rho, lattice );

    op_col[0] =  0.0;

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_minus = ( f[i] - f[i+1] ) * 0.5;

        double f_eq_minus = ( f_eq[i] - f_eq[i+1] ) * 0.5;
        
        double S_minus = ( S[i] - S[i+1] ) * 0.5;

        op_col[i] = ( f_eq_minus - f_minus ) * one_over_tau + S_minus;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_minus =  ( f[i] - f[i-1] ) * 0.5;

        double f_eq_minus = ( f_eq[i] - f_eq[i-1] ) * 0.5;
        
        double S_minus = ( S[i] - S[i-1] ) * 0.5;

        op_col[i] = ( f_eq_minus - f_minus ) * one_over_tau + S_minus;
    }

}

//====================================================================================================================//




//=============================== Calcula o operador BGK =============================================================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//====================================================================================================================//

#pragma acc routine seq
void op_bgk ( double f[nvel], double vx, double vy, double vz, double rho, double acc_x, double acc_y, double acc_z, 
				double tau, double op_col[nvel], LATTICE lattice )
{
    double f_eq[nvel];
    
    double S[nvel]{};

    double one_over_tau = 1.0 / tau;

    double vx_acc = vx + acc_x * 0.5;
    double vy_acc = vy + acc_y * 0.5;
    double vz_acc = vz + acc_z * 0.5;
    
    double Fx = acc_x * rho;
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;  

	if ( acc_x != 0 || acc_y != 0 || acc_z != 0 )
	{
		source ( Fx, Fy, Fz, vx_acc, vy_acc, vz_acc, rho, tau, S, lattice );  
	}

    dist_eq ( f_eq, vx_acc, vy_acc, vz_acc, rho, lattice );
	
	for ( int i = 0; i < nvel; i++ )	op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau + S[i];
}

//====================================================================================================================//



//=============================== Matrizes de transformacao do modelo MRT ============================================//
//
//  As matrizes sao constantes: ficam em escopo de arquivo e sao enviadas uma unica vez ao acelerador
//  por 'acc declare copyin'. Dentro de uma 'acc routine seq' o gcc rejeita dados estaticos declarados
//  na propria funcao ('requires a declare directive for use in a routine function'); alem disso, a
//  versao local reinicializava 722 doubles (D3Q19) a cada chamada da colisao na CPU.
//
//  Ordem dos momentos (D3Q19): rho, e, eps, jx, qx, jy, qy, jz, qz, 3pxx, 3pi_xx, p_ww, pi_ww,
//                              p_xy, p_yz, p_xz, mx, my, mz
//  Fonte: d'Humieres-Lallemand-Luo (2002), reordenada para a numeracao de def_lattice_d3q19.
//
//====================================================================================================================//

static const double MRT_D3Q19_M[19][19] = { 
{ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. },
{ -30., -11., -11., -11., -11., -11., -11., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8. },
{ 12., -4., -4., -4., -4., -4., -4., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. },
{ 0., 1., -1., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0. },
{ 0., -4., 4., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., 1., -1., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., -1., 1., -1., 1. },
{ 0., 0., 0., -4., 4., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., -1., 1., -1., 1. },
{ 0., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., 1., -1., -1., 1., -1., 1., 1., -1. },
{ 0., 0., 0., 0., 0., -4., 4., 0., 0., 0., 0., 1., -1., -1., 1., -1., 1., 1., -1. },
{ 0., 2., 2., -1., -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2. },
{ 0., -4., -4., 2., 2., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2. },
{ 0., 0., 0., 1., 1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., -2., -2., 2., 2., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1. },
{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., 1., -1., 1., -1., -1., 1., -1., 1., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., -1., 0., 0., 0., 0., -1., 1., -1., 1. },
{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 1., -1., -1., 1. }
};

static const double MRT_D3Q19_M_INV[19][19] = {
{ 1./19.,-5./399., 1./21., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 1./10.,-1./10., 0., 0., 0., 0., 1./18.,-1./18., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63.,-1./10., 1./10., 0., 0., 0., 0., 1./18.,-1./18., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0., 1./10.,-1./10., 0., 0.,-1./36., 1./36., 1./12.,-1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0.,-1./10., 1./10., 0., 0.,-1./36., 1./36., 1./12.,-1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0., 0., 0., 1./10.,-1./10.,-1./36., 1./36.,-1./12., 1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0., 0., 0.,-1./10., 1./10.,-1./36., 1./36.,-1./12., 1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40., 1./10., 1./40., 0., 0., 1./36., 1./72., 1./12., 1./24., 1./4., 0., 0., 1./8.,-1./8., 0. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40.,-1./10.,-1./40., 0., 0., 1./36., 1./72., 1./12., 1./24., 1./4., 0., 0.,-1./8., 1./8., 0. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40.,-1./10.,-1./40., 0., 0., 1./36., 1./72., 1./12., 1./24.,-1./4., 0., 0., 1./8., 1./8., 0. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40., 1./10., 1./40., 0., 0., 1./36., 1./72., 1./12., 1./24.,-1./4., 0., 0.,-1./8.,-1./8., 0. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40., 0., 0., 1./10., 1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0., 1./4.,-1./8., 0., 1./8. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40., 0., 0.,-1./10.,-1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0., 1./4., 1./8., 0.,-1./8. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40., 0., 0.,-1./10.,-1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0.,-1./4.,-1./8., 0.,-1./8. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40., 0., 0., 1./10., 1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0.,-1./4., 1./8., 0., 1./8. },
{ 1./19., 4./1197., 1./252., 0., 0.,-1./10.,-1./40.,-1./10.,-1./40.,-1./18.,-1./36., 0., 0., 0., 1./4., 0., 0.,-1./8., 1./8. },
{ 1./19., 4./1197., 1./252., 0., 0., 1./10., 1./40., 1./10., 1./40.,-1./18.,-1./36., 0., 0., 0., 1./4., 0., 0., 1./8.,-1./8. },
{ 1./19., 4./1197., 1./252., 0., 0.,-1./10.,-1./40., 1./10., 1./40.,-1./18.,-1./36., 0., 0., 0.,-1./4., 0., 0.,-1./8.,-1./8. },
{ 1./19., 4./1197., 1./252., 0., 0., 1./10., 1./40.,-1./10.,-1./40.,-1./18.,-1./36., 0., 0., 0.,-1./4., 0., 0., 1./8., 1./8. }
};

static const double MRT_D2Q9_M[9][9] = {
{ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },	// m0 = rho
{ -4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0 },	// m1 = e
{  4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0 },	// m2 = epsilon
{  0.0,  1.0, -1.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0 },	// m3 = jx
{  0.0, -2.0,  2.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0 },	// m4 = qx
{  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0 },	// m5 = jy
{  0.0,  0.0,  0.0, -2.0,  2.0,  1.0, -1.0, -1.0,  1.0 },	// m6 = qy
{  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0 },	// m7 = pxx
{  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0 }	// m8 = pxy
};

static const double MRT_D2Q9_M_INV[9][9] = {

{  1./9.,   -1./9.,   1./9.,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 },				// f0
{  1./9.,  -1./36., -1.0/18.0,  1./6.,  -1./6.,   0.0,   0.0,   1./4.,   0.0 },			// f1  
{  1./9.,  -1./36., -1.0/18.0, -1./6.,   1./6.,   0.0,   0.0,   1./4.,   0.0 },			// f2  
{  1./9.,  -1./36., -1.0/18.0,  0.0,   0.0,   1./6.,  -1./6.,  -1./4.,   0.0 },			// f3  
{  1./9.,  -1./36., -1.0/18.0,  0.0,   0.0,  -1./6.,   1./6.,  -1./4.,   0.0 },			// f4  
{  1./9.,   1./18.,   1./36.,   1./6.,   1./12.,   1./6.,   1./12.,  0.0,   1./4. },	// f5  
{  1./9.,   1./18.,   1./36.,  -1./6.,  -1./12.,  -1./6.,  -1./12.,  0.0,   1./4. },	// f6  
{  1./9.,   1./18.,   1./36.,   1./6.,   1./12.,  -1./6.,  -1./12.,  0.0,  -1./4. },	// f7  
{  1./9.,   1./18.,   1./36.,  -1./6.,  -1./12.,   1./6.,   1./12.,  0.0,  -1./4. }		// f8  
};

#pragma acc declare copyin( MRT_D3Q19_M, MRT_D3Q19_M_INV, MRT_D2Q9_M, MRT_D2Q9_M_INV )

//====================================================================================================================//


//=============================== Collision step using the MRT model =================================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_MRT ( double* f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{
	if ( dim == 3 )
	{
		const double s_nu = 1. / parameters.tau;
		
		const double s_q = 8. * ( 2. - s_nu ) / (8. - s_nu );
		
		const double s_e = 1.19;
		const double s_eps = 1.4; 	
		const double s_pi = 1.4;
		const double s_m = 1.98;
		
		double omega[19] = { 0., s_e, s_eps, 0., s_q, 0., s_q, 0., s_q, s_nu, s_pi, s_nu, s_pi, s_nu, s_nu, s_nu, 
									s_m, s_m, s_m };

		//static double omega[19] = { 0., 1.19, 1.4, 0., 1.2, 0., 1.2, 0., 1.2, s_nu, 1.4, s_nu, s_nu, s_nu, s_nu, s_nu, s_nu, 
		//						1.98, 1.98 };
		
		// Matriz de transformação G_M (moments  ←→ populations) para MRT-LBM D3Q19
		// Fonte: construção algébrica de d’Humières–Lallemand–Luo (2002)	
		// ordem dos momentos: ρ, e, ε, jx, qx, jy, qy, jz, qz, 3pxx, 3pi_ww, p_ww, pi_ww, p_xy, p_yz, p_xz, mx, my, mz
		
		if ( nvel == 19 )
		{

			// Inversa da matriz de transformação G_M  (MRT – D3Q19)

		 
			//------------ Distribuições de Equilíbrio -----------------------------------------------------------------------//
			
			double f_eq[nvel];
			
			double vx, vy, vz, rho;

			calcula ( f, vx, vy, vz, rho, lattice );

			dist_eq ( f_eq, vx, vy, vz, rho, lattice );
			
			//---------------- Termo de força (Guo) --------------------------------------------------------------------------//
			
			double S[nvel]{}; 
			
			if ( acc_x != 0 || acc_y != 0 || acc_z != 0 )
			{
					double one_ov_cs_sqd = lattice.one_over_c_s2;
					
					double v[3] = { vx, vy, vz };
					
					double F[3] = { acc_x * rho, acc_y * rho, acc_z * rho };
								
					double vF = dot_product ( F, v );
					
					double c_i[dim];

					for ( int i = 0; i < nvel; i++ )
					{
						c_i[0] = lattice.c_i[ i * dim + 0 ];
						c_i[1] = lattice.c_i[ i * dim + 1 ];
						c_i[2] = lattice.c_i[ i * dim + 2 ];           
						  
						double cF = dot_product ( c_i, F );
						
						double cv = dot_product ( c_i, v );

						S[i] = lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
					}
			}

			//---------------- Colisão no espaço dos momentos --------------------------------------------------------//

			double moment[nvel];
			double mom_eq[nvel];
			double F_mom[nvel];

			for ( int lin = 0; lin < nvel; lin++ )
			{
				moment[lin] = 0.0;
				mom_eq[lin] = 0.0;
				F_mom[lin] = 0.0;
				
				for ( int i = 0; i < nvel; i++ )
				{
					moment[lin] = moment[lin] + ( MRT_D3Q19_M[lin][i] ) * f[i];			
					mom_eq[lin] = mom_eq[lin] + ( MRT_D3Q19_M[lin][i] ) * f_eq[i];
					F_mom[lin]  = F_mom[lin]  + ( MRT_D3Q19_M[lin][i] ) * S[i];			
				}
				
				moment[lin] = moment[lin] - omega[lin] * ( moment[lin] - mom_eq[lin] ) + ( 1. - 0.5 * omega[lin] ) * F_mom[lin];
			}
			
			//---------------- Retorna do espaço dos momentos ----------------------------------------------------------------//
			
			for ( int i = 0; i < nvel; i++ )
			{
				
				f[i] = 0.;
				
				for ( int col = 0; col < nvel; col++ )
				{
					f[i] = f[i] + MRT_D3Q19_M_INV[i][col] * moment[col];
				}
			}
		}
	}
	else
	{
		// taxas de relaxação 
		const double s0 = 0.0;					// rho (conservado)
		const double s1 = 1.6;					// e
		const double s2 = 1.8;					// epsilon
		const double s3 = 0.0;					// jx (conservado)
		const double s4 = 1.9;					// qx
		const double s5 = 0.0;					// jy (conservado)
		const double s6 = 1.9;					// qy
		const double s7 = 1. / parameters.tau;	// pxx (shear)
		const double s8 = 1. / parameters.tau;	// pxy (shear)
		 
		const double Omega[9][9] = {
			{ s0, 0., 0., 0., 0., 0., 0., 0., 0. },	
			{ 0., s1, 0., 0., 0., 0., 0., 0., 0. },	
			{ 0., 0., s2, 0., 0., 0., 0., 0., 0. },	
			{ 0., 0., 0., s3, 0., 0., 0., 0., 0. },	
			{ 0., 0., 0., 0., s4, 0., 0., 0., 0. },	
			{ 0., 0., 0., 0., 0., s5, 0., 0., 0. },	
			{ 0., 0., 0., 0., 0., 0., s6, 0., 0. },	
			{ 0., 0., 0., 0., 0., 0., 0., s7, 0. },	
			{ 0., 0., 0., 0., 0., 0., 0., 0., s8 }	
		};
	
		// Matriz de transformação 

		// Inversa da matriz de transformação 
			
		//------------ Distribuições de Equilíbrio -------------------------------------------------------------------//
		
		double f_eq[nvel];
		
		double vx, vy, vz, rho;

		calcula ( f, vx, vy, vz, rho, lattice );

		dist_eq ( f_eq, vx, vy, vz, rho, lattice );
		
		//------------ Termo de força --------------------------------------------------------------------------------//

		double S_force[nvel] = {0.0};
		
		if (acc_x != 0.0 || acc_y != 0.0 || acc_z != 0.0)
		{
			const double one_ov_cs_sqd = lattice.one_over_c_s2;

			double v[2] = { vx, vy };
			double F[2] = { acc_x * rho, acc_y * rho };

			const double vF = dot_product(F, v);

			double c_i[2];
			
			for (int i = 0; i < nvel; ++i)
			{
				c_i[0] = lattice.c_i[i*dim + 0];
				c_i[1] = lattice.c_i[i*dim + 1];

				const double cF = dot_product(c_i, F);
				const double cv = dot_product(c_i, v);

				S_force[i] = lattice.w[i] * one_ov_cs_sqd *
							 ( cF + one_ov_cs_sqd * cF * cv - vF );
			}
		}
		
		//------------- Transforma para o espaço dos momentos --------------------------------------------------------//

		double moment[nvel];		
		double mom_eq[nvel];		
		double F_mom[nvel];

		for (int lin = 0; lin < nvel; ++lin)
		{
			moment[lin] = 0.0;			
			mom_eq[lin] = 0.0;			
			F_mom[lin]  = 0.0;
			
			for (int col = 0; col < nvel; ++col)
			{
				moment[lin] = moment[lin] + MRT_D2Q9_M[lin][col] * f[col];
				mom_eq[lin] = mom_eq[lin] + MRT_D2Q9_M[lin][col] * f_eq[col];
				F_mom[lin]  = F_mom[lin]  + MRT_D2Q9_M[lin][col] * S_force[col];
			}
		}

		//---------------- Colisão em espaço de momentos -------------------------------------------------------------//

		for (int lin = 0; lin < nvel; ++lin)
		{
			const double s = Omega[lin][lin];   // matriz diagonal → pega só o elemento da diagonal

			moment[lin] = moment[lin] - s * (moment[lin] - mom_eq[lin])  + (1.0 - 0.5 * s) * F_mom[lin];
		}

		//---------------- Volta para o espaço das distribuições -----------------------------------------------------//

		for (int lin = 0; lin < nvel; ++lin)
		{
			f[lin] = 0.0;
			
			for (int col = 0; col < nvel; ++col)
			{
				f[lin] = f[lin] + MRT_D2Q9_M_INV[lin][col] * moment[col];
			}
		}
	}

}

//====================================================================================================================//




//=============================== Collision step using the MRT model for Stokes equation =============================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_MRT_stokes ( double* f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{
	
	const double s_nu = 1. / parameters.tau;
	
	const double s_q = 8. * ( 2. - s_nu ) / (8. - s_nu );
	
	const double s_e = 1.19;
	const double s_eps = 1.4; 	
	const double s_pi = 1.4;
	const double s_m = 1.98;
	
	double omega[19] = { 0., s_e, s_eps, 0., s_q, 0., s_q, 0., s_q, s_nu, s_pi, s_nu, s_pi, s_nu, s_nu, s_nu, 
								s_m, s_m, s_m };

	//static double omega[19] = { 0., 1.19, 1.4, 0., 1.2, 0., 1.2, 0., 1.2, s_nu, 1.4, s_nu, s_nu, s_nu, s_nu, s_nu, s_nu, 
	//						1.98, 1.98 };
	
	// Matriz de transformação G_M (moments  ←→ populations) para MRT-LBM D3Q19
	// Fonte: construção algébrica de d’Humières–Lallemand–Luo (2002)	
	// ordem dos momentos: ρ, e, ε, jx, qx, jy, qy, jz, qz, 3pxx, 3pi_ww, p_ww, pi_ww, p_xy, p_yz, p_xz, mx, my, mz
	
	if ( nvel == 19 )
	{

		// Inversa da matriz de transformação G_M  (MRT – D3Q19)

	 
		//------------ momentos de Equilíbrio -----------------------------------------------------------------------//
			
		double vx, vy, vz, rho;

		calcula ( f, vx, vy, vz, rho, lattice );
		
		double j_x = rho * vx;
		double j_y = rho * vy;
		double j_z = rho * vz;
		
		double mom_eq[nvel];
		
		mom_eq[0] = rho;						// 	ρ
		mom_eq[1] = -11. * rho;					//	e
		mom_eq[2] = 3. * rho;					//	ε
		mom_eq[3] =	j_x;						//	jx 
		mom_eq[4] =	-(2./3.) * j_x;				//	qx
		mom_eq[5] = j_y;						//	jy
		mom_eq[6] = -(2./3.) * j_y;				//	qy
		mom_eq[7] = j_z;						// 	jz
		mom_eq[8] = -(2./3.) * j_z;				// 	qz
		mom_eq[9] = 0.;							// 	3pxx
		mom_eq[10] = 0.;						// 	3pi_ww
		mom_eq[11] = 0.;						// 	p_ww
		mom_eq[12] = 0.;						// 	pi_ww
		mom_eq[13] = 0.;						// 	p_xy
		mom_eq[14] = 0.;						// 	p_yz
		mom_eq[15] = 0.;						// 	p_xz
		mom_eq[16] = 0.;						// 	mx
		mom_eq[17] = 0.;						// 	my
		mom_eq[18] = 0.;						// 	mz
			
		//---------------- Termo de força (Guo) --------------------------------------------------------------------------//
		
		double S[nvel]{}; 
		
		if ( acc_x != 0 || acc_y != 0 || acc_z != 0 )
		{
			double one_ov_cs_sqd = lattice.one_over_c_s2;
			
			double v[3] = { vx, vy, vz };
			
			double F[3] = { acc_x * rho, acc_y * rho, acc_z * rho };
						
			double vF = dot_product ( F, v );
			
			double c_i[dim];

			for ( int i = 0; i < nvel; i++ )
			{
				c_i[0] = lattice.c_i[ i * dim + 0 ];
				c_i[1] = lattice.c_i[ i * dim + 1 ];
				c_i[2] = lattice.c_i[ i * dim + 2 ];           
				  
				double cF = dot_product ( c_i, F );
				
				double cv = dot_product ( c_i, v );

				S[i] = lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
			}
		}

		//---------------- Colisão no espaço dos momentos ----------------------------------------------------------------//

		double moment[nvel];
		
		double F_mom[nvel];

		// Atencao: mom_eq[] ja contem os momentos de equilibrio linearizados (Stokes) calculados acima.
		// Zera-lo dentro deste laco relaxaria e, epsilon e q_alpha para zero em vez de -11 rho, 3 rho
		// e -(2/3) j_alpha, destruindo a parte de pressao/energia do modelo.

		for ( int lin = 0; lin < nvel; lin++ )
		{
			moment[lin] = 0.0;
			F_mom[lin] = 0.0;

			for ( int i = 0; i < nvel; i++ )
			{
				moment[lin] = moment[lin] + ( MRT_D3Q19_M[lin][i] ) * f[i];
				F_mom[lin]  = F_mom[lin]  + ( MRT_D3Q19_M[lin][i] ) * S[i];
			}
			
			moment[lin] = moment[lin] - omega[lin] * ( moment[lin] - mom_eq[lin] ) + ( 1. - 0.5 * omega[lin] ) * F_mom[lin];
		}
		
		//---------------- Retorna do espaço dos momentos ----------------------------------------------------------------//
		
		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = 0.;
			
			for ( int col = 0; col < nvel; col++ )
			{
				f[i] = f[i] + MRT_D3Q19_M_INV[i][col] * moment[col];
			}
		}
	}
	
}

//====================================================================================================================//




//========= Regularized BGK collision - J. Latt, B. Chopard / Mathematics and Computers in Simulation 72 (2006) ======//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//


#pragma acc routine seq
void coll_Reg ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   	
	double tau = parameters.tau;
	 
    double S[nvel];

    double vx, vy, vz, rho;

    double one_over_tau = 1.0 / tau;

    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    if ( dim == 2 )
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    

		double Pi_xy = 0.0;    

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 3 ] * Pi_yy 
					 + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy );         
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - one_over_tau ) * f_1[i] + S[i];
	}
	
	else
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    
		double Pi_zz = 0.0;

		double Pi_xy = 0.0;    
		double Pi_xz = 0.0;    
		double Pi_yz = 0.0;

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			Pi_zz = Pi_zz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 2 ] * lattice.c_i[ i*dim + 2 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
			Pi_xz = Pi_xz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 2 ];        
			Pi_yz = Pi_yz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 2 ];
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 4 ] * Pi_yy 
					+ lattice.Q_i[ ind + 8 ] * Pi_zz + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy 
					+ 2 * lattice.Q_i[ ind + 2 ] * Pi_xz + 2 * lattice.Q_i[ ind + 5 ] * Pi_yz );         
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - one_over_tau ) * f_1[i] + S[i];
	}
	
}

//====================================================================================================================//




//=========================== Regularized TRT collision ==============================================================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//


#pragma acc routine seq
void coll_Reg_TRT ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   	
	double tau = parameters.tau;
	
	const double inv_tau = 1.0 / tau;

    const double tau_anti = ( ( 8.0 - inv_tau ) / ( 8.0 * (  2.0 - inv_tau ) ) ); 
    
    const double inv_tau_anti = 1.0 / tau_anti;
	 
    double S[nvel];
    
    double op_col_simm[nvel];
    double op_col_anti[nvel];    

    double vx, vy, vz, rho;

    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    if ( dim == 2 )
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    

		double Pi_xy = 0.0;    

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 3 ] * Pi_yy 
					 + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy );         
		}
		
		//---------------- Computes the simmetrical and anti-simmetrical terms --------------------------//
		
		op_col_simm[0] = f_eq[0] + ( 1.0 - inv_tau ) * f_1[0] + S[0];
		op_col_anti[0] = 0;

		for ( int i = 1; i < nvel; i = i + 2 )
		{
			double f_eq_simm = ( f_eq[i] + f_eq[i+1] ) * 0.5;
			double f_eq_anti = ( f_eq[i] - f_eq[i+1] ) * 0.5;
			
			double f_1_simm = ( f_1[i] + f_1[i+1] ) * 0.5;
			double f_1_anti = ( f_1[i] - f_1[i+1] ) * 0.5;
			
			double S_simm = ( S[i] + S[i+1] ) * 0.5;
			double S_anti = ( S[i] - S[i+1] ) * 0.5;

			op_col_simm[i] = f_eq_simm + ( 1.0 - inv_tau ) * f_1_simm + S_simm;
			op_col_anti[i] = f_eq_anti + ( 1.0 - inv_tau_anti ) * f_1_anti + S_anti;
		}

		for ( int i = 2; i < nvel; i = i + 2 )
		{
			double f_eq_simm = ( f_eq[i] + f_eq[i-1] ) * 0.5;
			double f_eq_anti = ( f_eq[i] - f_eq[i-1] ) * 0.5;
			
			double f_1_simm = ( f_1[i] + f_1[i-1] ) * 0.5;
			double f_1_anti = ( f_1[i] - f_1[i-1] ) * 0.5;
			
			double S_simm = ( S[i] + S[i-1] ) * 0.5;
			double S_anti = ( S[i] - S[i-1] ) * 0.5;
			
			op_col_simm[i] = f_eq_simm + ( 1.0 - inv_tau ) * f_1_simm + S_simm;
			op_col_anti[i] = f_eq_anti + ( 1.0 - inv_tau_anti ) * f_1_anti + S_anti;
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )	f[i] = op_col_simm[i] + op_col_anti[i];
	}
	
	else
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    
		double Pi_zz = 0.0;

		double Pi_xy = 0.0;    
		double Pi_xz = 0.0;    
		double Pi_yz = 0.0;

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			Pi_zz = Pi_zz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 2 ] * lattice.c_i[ i*dim + 2 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
			Pi_xz = Pi_xz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 2 ];        
			Pi_yz = Pi_yz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 2 ];
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 4 ] * Pi_yy 
					+ lattice.Q_i[ ind + 8 ] * Pi_zz + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy 
					+ 2 * lattice.Q_i[ ind + 2 ] * Pi_xz + 2 * lattice.Q_i[ ind + 5 ] * Pi_yz );         
		}
		
		//---------------- Computes the simmetrical and anti-simmetrical terms --------------------------//
		
		op_col_simm[0] = f_eq[0] + ( 1.0 - inv_tau ) * f_1[0] + S[0];
		op_col_anti[0] = 0;

		for ( int i = 1; i < nvel; i = i + 2 )
		{
			double f_eq_simm = ( f_eq[i] + f_eq[i+1] ) * 0.5;
			double f_eq_anti = ( f_eq[i] - f_eq[i+1] ) * 0.5;
			
			double f_1_simm = ( f_1[i] + f_1[i+1] ) * 0.5;
			double f_1_anti = ( f_1[i] - f_1[i+1] ) * 0.5;
			
			double S_simm = ( S[i] + S[i+1] ) * 0.5;
			double S_anti = ( S[i] - S[i+1] ) * 0.5;

			op_col_simm[i] = f_eq_simm + ( 1.0 - inv_tau ) * f_1_simm + S_simm;
			op_col_anti[i] = f_eq_anti + ( 1.0 - inv_tau_anti ) * f_1_anti + S_anti;
		}

		for ( int i = 2; i < nvel; i = i + 2 )
		{
			double f_eq_simm = ( f_eq[i] + f_eq[i-1] ) * 0.5;
			double f_eq_anti = ( f_eq[i] - f_eq[i-1] ) * 0.5;
			
			double f_1_simm = ( f_1[i] + f_1[i-1] ) * 0.5;
			double f_1_anti = ( f_1[i] - f_1[i-1] ) * 0.5;
			
			double S_simm = ( S[i] + S[i-1] ) * 0.5;
			double S_anti = ( S[i] - S[i-1] ) * 0.5;
			
			op_col_simm[i] = f_eq_simm + ( 1.0 - inv_tau ) * f_1_simm + S_simm;
			op_col_anti[i] = f_eq_anti + ( 1.0 - inv_tau_anti ) * f_1_anti + S_anti;
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )	f[i] = op_col_simm[i] + op_col_anti[i];
	}
	
}

//====================================================================================================================//



//============================ Over Regularized BGK collision =======================================================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//


#pragma acc routine seq
void coll_Ext_Reg ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   
	double tau = parameters.tau;
	
	double tau_gt1 = parameters.tau_nh;
	 
    double S[nvel];

    double vx, vy, vz, rho;
    
    double omega_1 = 1.0 - 1.0 / tau;
        
    double omega_gt1 = 1.0 - 1.0 / tau_gt1;
    
    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    if ( dim == 2 )
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    

		double Pi_xy = 0.0;    

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 3 ] * Pi_yy 
					 + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy ); 
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Collision ---------------------------------------------------------------------//
	
		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = f_eq[i] + omega_1 * f_1[i] + omega_gt1 * f_gt1[i] + S[i];
		}
	}
	else
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    
		double Pi_zz = 0.0;

		double Pi_xy = 0.0;    
		double Pi_xz = 0.0;    
		double Pi_yz = 0.0;

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			Pi_zz = Pi_zz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 2 ] * lattice.c_i[ i*dim + 2 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
			Pi_xz = Pi_xz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 2 ];        
			Pi_yz = Pi_yz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 2 ];
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 4 ] * Pi_yy 
					+ lattice.Q_i[ ind + 8 ] * Pi_zz + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy 
					+ 2 * lattice.Q_i[ ind + 2 ] * Pi_xz + 2 * lattice.Q_i[ ind + 5 ] * Pi_yz );  
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Collision ---------------------------------------------------------------------//
	
		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = f_eq[i] + omega_1 * f_1[i] + omega_gt1 * f_gt1[i] + S[i];
		}
	}	
}

//====================================================================================================================//




//============================ Over Regularized BGK collision =======================================================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//


#pragma acc routine seq
void op_Ext_Reg ( double f[nvel], double vx, double vy, double vz, double rho, double acc_x, double acc_y, double acc_z, 
					double tau, double op_col[nvel], LATTICE lattice, PARAMETERS parameters )
{   	
	double tau_gt1 = parameters.tau_nh;
	 
    double S[nvel];
    
    double one_over_tau = 1.0 / tau;
        
    double one_over_tau_nh = 1.0 / tau_gt1;
    
    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    if ( dim == 2 )
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    

		double Pi_xy = 0.0;    

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 3 ] * Pi_yy 
					 + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy ); 
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Compute the collision operator  -------------------------------------------------------------//
	
		for ( int i = 0; i < nvel; i++ ) op_col[i] = S[i] - ( one_over_tau * f_1[i] + one_over_tau_nh * f_gt1[i] );
		
	}
	else
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    
		double Pi_zz = 0.0;

		double Pi_xy = 0.0;    
		double Pi_xz = 0.0;    
		double Pi_yz = 0.0;

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			Pi_zz = Pi_zz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 2 ] * lattice.c_i[ i*dim + 2 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
			Pi_xz = Pi_xz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 2 ];        
			Pi_yz = Pi_yz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 2 ];
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 4 ] * Pi_yy 
					+ lattice.Q_i[ ind + 8 ] * Pi_zz + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy 
					+ 2 * lattice.Q_i[ ind + 2 ] * Pi_xz + 2 * lattice.Q_i[ ind + 5 ] * Pi_yz );  
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Compute the collision operator  -------------------------------------------------------------//
	
		for ( int i = 0; i < nvel; i++ ) op_col[i] = S[i] - ( one_over_tau * f_1[i] + one_over_tau_nh * f_gt1[i] );
		
	}	
}

//====================================================================================================================//




//============================ Over Regularized BGK collision =======================================================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//


#pragma acc routine seq
void coll_Ext_Reg_fourth ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   
	double tau = parameters.tau;
	
	double tau_gt1 = parameters.tau_nh;
	
    double S[nvel];

    double vx, vy, vz, rho;
    
    double omega_1 = 1.0 - 1.0 / tau;
        
    double omega_gt1 = 1.0 - 1.0 / tau_gt1;
    
    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq_fourth ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    if ( dim == 2 )
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    

		double Pi_xy = 0.0;    

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 3 ] * Pi_yy 
					 + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy ); 
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = f_eq[i] + omega_1 * f_1[i] + omega_gt1 * f_gt1[i] + S[i];
		}
	}
	else
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    
		double Pi_zz = 0.0;

		double Pi_xy = 0.0;    
		double Pi_xz = 0.0;    
		double Pi_yz = 0.0;

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			Pi_zz = Pi_zz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 2 ] * lattice.c_i[ i*dim + 2 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
			Pi_xz = Pi_xz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 2 ];        
			Pi_yz = Pi_yz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 2 ];
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 4 ] * Pi_yy 
					+ lattice.Q_i[ ind + 8 ] * Pi_zz + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy 
					+ 2 * lattice.Q_i[ ind + 2 ] * Pi_xz + 2 * lattice.Q_i[ ind + 5 ] * Pi_yz );  
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = f_eq[i] + omega_1 * f_1[i] + omega_gt1 * f_gt1[i] + S[i];
		}
	
	}
	
}

//====================================================================================================================//



//============================ Extended Regularized collision ========================================================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//


#pragma acc routine seq
void coll_Ext_Reg_sixth ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   
	double tau = parameters.tau;
	
	double tau_gt1 = parameters.tau_nh;
	
    double S[nvel];

    double vx, vy, vz, rho;
    
    double omega_1 = 1.0 - 1.0 / tau;
        
    double omega_gt1 = 1.0 - 1.0 / tau_gt1;
    
    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq_sixth ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    if ( dim == 2 )
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    

		double Pi_xy = 0.0;    

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 3 ] * Pi_yy 
					 + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy ); 
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = f_eq[i] + omega_1 * f_1[i] + omega_gt1 * f_gt1[i] + S[i];
		}
	}
	else
	{
		//----------- Calculates the non-equilibrium flow tensor ------------------------------------//

		double Pi_xx = 0.0;
		double Pi_yy = 0.0;    
		double Pi_zz = 0.0;

		double Pi_xy = 0.0;    
		double Pi_xz = 0.0;    
		double Pi_yz = 0.0;

		for ( int i = 0; i < nvel; i++ ) 
		{
			Pi_xx = Pi_xx + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim     ];
			Pi_yy = Pi_yy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 1 ];
			Pi_zz = Pi_zz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 2 ] * lattice.c_i[ i*dim + 2 ];
			
			Pi_xy = Pi_xy + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 1 ];        
			Pi_xz = Pi_xz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim     ] * lattice.c_i[ i*dim + 2 ];        
			Pi_yz = Pi_yz + ( f[i] - f_eq[i] ) * lattice.c_i[ i*dim + 1 ] * lattice.c_i[ i*dim + 2 ];
		}

		//----------- Calculates f(1) -------------------------------------------------------------------//

		const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

		double f_1[nvel];

		double f_gt1[nvel];	

		for ( int i = 0; i < nvel; i ++ )
		{
			int ind = i * dim * dim;
										
			f_1[i] = one_over_2cs4 * lattice.w[i] * ( lattice.Q_i[ ind ] * Pi_xx + lattice.Q_i[ ind + 4 ] * Pi_yy 
					+ lattice.Q_i[ ind + 8 ] * Pi_zz + 2 * lattice.Q_i[ ind + 1 ] * Pi_xy 
					+ 2 * lattice.Q_i[ ind + 2 ] * Pi_xz + 2 * lattice.Q_i[ ind + 5 ] * Pi_yz );  
			
			f_gt1[i] = f[i] - f_eq[i] - f_1[i];
		}
		
		//-------------- Collision ---------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			
			f[i] = f_eq[i] + omega_1 * f_1[i] + omega_gt1 * f_gt1[i] + S[i];
		}
	
	}
	
}

//====================================================================================================================//



//========= Onsager Regularized BGK collision - Jonnalagadda, Sharma & Agrawa (PHYS. REV. E 104, 015313 (2021) =======//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_Onsager_Reg ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   
	double tau = parameters.tau;
		 
    double S[nvel];

    double vx, vy, vz, rho;

    double one_over_tau = 1.0 / tau;

    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    //----------- Calculates the non-equilibrium distribution and the peculiar velocity ----------//

    double f_neq[nvel];

	double C_i[ dim * nvel ]; 
	
	if ( dim == 3 )
	{
		for ( int i = 0; i < nvel; i++ )
		{
			f_neq[i] = f[i] - f_eq[i];
			
			C_i[ i * dim + 0 ] = lattice.c_i[ i * dim + 0 ] - vx_alt;
			C_i[ i * dim + 1 ] = lattice.c_i[ i * dim + 1 ] - vy_alt;
			C_i[ i * dim + 2 ] = lattice.c_i[ i * dim + 2 ] - vz_alt;
		}
	}
	
	if ( dim == 2 )
	{
		for ( int i = 0; i < nvel; i++ )
		{
			f_neq[i] = f[i] - f_eq[i];
			
			C_i[ i * dim + 0 ] = lattice.c_i[ i * dim + 0 ] - vx_alt;
			C_i[ i * dim + 1 ] = lattice.c_i[ i * dim + 1 ] - vy_alt;
		}
	}
	
    //----------- Calculates f(1) -------------------------------------------------------------------//

	const double gamma = 1.0; // adiabatic index
	
    double factor = 0.5 * gamma * lattice.one_over_c_s2 * lattice.one_over_c_s2 / rho ;    
            
    double ci_ci[dim * dim ];
    
    //  CORRECAO: e um acumulador (sum_ci_ci[j] = sum_ci_ci[j] + ...) e precisa comecar em zero.
    //  Sem a inicializacao o operador partia do lixo da pilha: no equilibrio dava f != f_eq e
    //  a colisao deixava de conservar massa e quantidade de movimento.

    double sum_ci_ci[dim * dim ] = {};
    
    double GAMMA_i[dim * dim * nvel ];
	
	double Ci_Ci[dim * dim ];
         
	for ( int i = 0; i < nvel; i ++ )
	{
		outer_product( &C_i[ i * dim ], &C_i[ i * dim ], Ci_Ci );
		
		double Ci_2 = dot_product( &C_i[ i * dim ], &C_i[ i * dim ] );
		
		outer_product( &lattice.c_i[i * dim], &lattice.c_i[i * dim], ci_ci );
		
		double ci_dot_ci = dot_product( &lattice.c_i[i * dim],  &lattice.c_i[i * dim] );
		
		for ( int j = 0; j < dim; j++ )
		{
			for ( int k = 0; k < dim; k++ )
			{
				sum_ci_ci[j + k * dim] = sum_ci_ci[j + k * dim] 
										+  ( ci_ci[j + k * dim] - ci_dot_ci * ( j == k ) / dim ) * f_neq[i];
										
				GAMMA_i[ j + k * dim + i * dim * dim ] = Ci_Ci[ j + k * dim ] - 0.5 * ( gamma - 1 ) * Ci_2 * ( j == k );
			}
		}		
	}
	
    double f_1[nvel];

    for ( int i = 0; i < nvel; i ++ )
    {		
		f_1[i] =  factor * f_eq[i] * dot_dot_product( &GAMMA_i[i * dim * dim], sum_ci_ci );         
    }

    //-------------- Collision ---------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - one_over_tau ) * f_1[i] + S[i];
    
    //----------------------------------------------------------------------------------------------//
}
//====================================================================================================================//




//========= Onsager Regularized BGK collision - Jonnalagadda, Sharma & Agrawa (PHYS. REV. E 104, 015313 (2021) =======//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_Onsager_ExtReg ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{   
	double tau = parameters.tau;
	
	double tau_gt1 = parameters.tau_nh;
	 	 
    double S[nvel];

    double vx, vy, vz, rho;

    double one_over_tau = 1.0 / tau;
    
    double omega_gt1 = 1.0 - 1.0 / tau_gt1;

    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    //----------- Calculates the non-equilibrium distribution and the peculiar velocity ----------//

    double f_neq[nvel];

	double C_i[ dim * nvel ]; 

	if ( dim == 3 )
	{
		for ( int i = 0; i < nvel; i++ )
		{
			f_neq[i] = f[i] - f_eq[i];
			
			C_i[ i * dim + 0 ] = lattice.c_i[ i * dim + 0 ] - vx_alt;
			C_i[ i * dim + 1 ] = lattice.c_i[ i * dim + 1 ] - vy_alt;
			C_i[ i * dim + 2 ] = lattice.c_i[ i * dim + 2 ] - vz_alt;
		}
	}
	
	if ( dim == 2 )
	{
		for ( int i = 0; i < nvel; i++ )
		{
			f_neq[i] = f[i] - f_eq[i];
			
			C_i[ i * dim + 0 ] = lattice.c_i[ i * dim + 0 ] - vx_alt;
			C_i[ i * dim + 1 ] = lattice.c_i[ i * dim + 1 ] - vy_alt;
		}
	}
		
    //----------- Calculates f(1) -------------------------------------------------------------------//

	const double gamma = 1.0; // adiabatic index
	
    double factor = 0.5 * gamma * lattice.one_over_c_s2 * lattice.one_over_c_s2 / rho ;    
            
    double ci_ci[dim * dim ];
    
    //  CORRECAO: e um acumulador (sum_ci_ci[j] = sum_ci_ci[j] + ...) e precisa comecar em zero.
    //  Sem a inicializacao o operador partia do lixo da pilha: no equilibrio dava f != f_eq e
    //  a colisao deixava de conservar massa e quantidade de movimento.

    double sum_ci_ci[dim * dim ] = {};
    
    double GAMMA_i[dim * dim * nvel ];
	
	double Ci_Ci[dim * dim ];
         
	for ( int i = 0; i < nvel; i ++ )
	{
		outer_product( &C_i[ i * dim ], &C_i[ i * dim ], Ci_Ci );
		
		double Ci_2 = dot_product( &C_i[ i * dim ], &C_i[ i * dim ] );
		
		outer_product( &lattice.c_i[i * dim], &lattice.c_i[i * dim], ci_ci );
		
		double ci_dot_ci = dot_product( &lattice.c_i[i * dim],  &lattice.c_i[i * dim] );
		
		for ( int j = 0; j < dim; j++ )
		{
			for ( int k = 0; k < dim; k++ )
			{
				sum_ci_ci[j + k * dim] = sum_ci_ci[j + k * dim] 
										+  ( ci_ci[j + k * dim] - ci_dot_ci * ( j == k ) / dim ) * f_neq[i];
										
				GAMMA_i[ j + k * dim + i * dim * dim ] = Ci_Ci[ j + k * dim ] - 0.5 * ( gamma - 1 ) * Ci_2 * ( j == k );
			}
		}		
	}
	
    double f_1[nvel];
    
    double f_gt1[nvel];	

    for ( int i = 0; i < nvel; i ++ )
    {		
		f_1[i] =  factor * f_eq[i] * dot_dot_product( &GAMMA_i[i * dim * dim], sum_ci_ci );     
		
		f_gt1[i] = f[i] - f_eq[i] - f_1[i];
    }

    //-------------- Collision ---------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - one_over_tau ) * f_1[i] + omega_gt1 * f_gt1[i]  + S[i];
    
    //----------------------------------------------------------------------------------------------//
}

//====================================================================================================================//





//========= Onsager Regularized BGK collision - Jonnalagadda, Sharma & Agrawa (PHYS. REV. E 104, 015313 (2021) =======//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_Onsager_ExtReg_fourth ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, 
									PARAMETERS parameters )
{   
	double tau = parameters.tau;
	
	double tau_gt1 = parameters.tau_nh;
	
    double S[nvel];

    double vx, vy, vz, rho;

    double one_over_tau = 1.0 / tau;
    
    double omega_gt1 = 1.0 - 1.0 / tau_gt1;

    calcula ( f, vx, vy, vz, rho, lattice );

    double vx_alt = vx + 0.5 * acc_x;
    double vy_alt = vy + 0.5 * acc_y;
    double vz_alt = vz + 0.5 * acc_z;
    
    double Fx = acc_x * rho;    
    double Fy = acc_y * rho;
    double Fz = acc_z * rho;
    
    source ( Fx, Fy,Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

    //----------- Calculates the equilibrium distribution ----------------------------------------//

    double f_eq[nvel];

    dist_eq_fourth( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
    
    //----------- Calculates the non-equilibrium distribution and the peculiar velocity ----------//

    double f_neq[nvel];

	double C_i[ dim * nvel ]; 

	if ( dim == 3 )
	{
		for ( int i = 0; i < nvel; i++ )
		{
			f_neq[i] = f[i] - f_eq[i];
			
			C_i[ i * dim + 0 ] = lattice.c_i[ i * dim + 0 ] - vx_alt;
			C_i[ i * dim + 1 ] = lattice.c_i[ i * dim + 1 ] - vy_alt;
			C_i[ i * dim + 2 ] = lattice.c_i[ i * dim + 2 ] - vz_alt;
		}
	}
	
	if ( dim == 2 )
	{
		for ( int i = 0; i < nvel; i++ )
		{
			f_neq[i] = f[i] - f_eq[i];
			
			C_i[ i * dim + 0 ] = lattice.c_i[ i * dim + 0 ] - vx_alt;
			C_i[ i * dim + 1 ] = lattice.c_i[ i * dim + 1 ] - vy_alt;
		}
	}
	
    //----------- Calculates f(1) -------------------------------------------------------------------//

	const double gamma = 1.0; // adiabatic index
	
    double factor = 0.5 * gamma * lattice.one_over_c_s2 * lattice.one_over_c_s2 / rho ;    
            
    double ci_ci[dim * dim ];
    
    //  CORRECAO: e um acumulador (sum_ci_ci[j] = sum_ci_ci[j] + ...) e precisa comecar em zero.
    //  Sem a inicializacao o operador partia do lixo da pilha: no equilibrio dava f != f_eq e
    //  a colisao deixava de conservar massa e quantidade de movimento.

    double sum_ci_ci[dim * dim ] = {};
    
    double GAMMA_i[dim * dim * nvel ];
	
	double Ci_Ci[dim * dim ];
         
	for ( int i = 0; i < nvel; i ++ )
	{
		outer_product( &C_i[ i * dim ], &C_i[ i * dim ], Ci_Ci );
		
		double Ci_2 = dot_product( &C_i[ i * dim ], &C_i[ i * dim ] );
		
		outer_product( &lattice.c_i[i * dim], &lattice.c_i[i * dim], ci_ci );
		
		double ci_dot_ci = dot_product( &lattice.c_i[i * dim],  &lattice.c_i[i * dim] );
		
		for ( int j = 0; j < dim; j++ )
		{
			for ( int k = 0; k < dim; k++ )
			{
				sum_ci_ci[j + k * dim] = sum_ci_ci[j + k * dim] 
										+  ( ci_ci[j + k * dim] - ci_dot_ci * ( j == k ) / dim ) * f_neq[i];
										
				GAMMA_i[ j + k * dim + i * dim * dim ] = Ci_Ci[ j + k * dim ] - 0.5 * ( gamma - 1 ) * Ci_2 * ( j == k );
			}
		}		
	}
	
    double f_1[nvel];
    
    double f_gt1[nvel];	

    for ( int i = 0; i < nvel; i ++ )
    {		
		f_1[i] =  factor * f_eq[i] * dot_dot_product( &GAMMA_i[i * dim * dim], sum_ci_ci );     
		
		f_gt1[i] = f[i] - f_eq[i] - f_1[i];
    }

    //-------------- Collision ---------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - one_over_tau ) * f_1[i] + omega_gt1 * f_gt1[i]  + S[i];
    
    //----------------------------------------------------------------------------------------------//
}

//====================================================================================================================//



//=============================== Collision step using the Santos, Facin & Philippi model (immicisble) ===============//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_SFP ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, double gx_B, double gy_B,
                double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int = parameters.A_fact;
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;

    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );
   
    double vx_alt_R = vx_R + fat_int * vx_m_uni;
    double vy_alt_R = vy_R + fat_int * vy_m_uni;
    double vz_alt_R = vz_R + fat_int * vz_m_uni;

    double vx_alt_B = vx_B - fat_int * vx_m_uni;
    double vy_alt_B = vy_B - fat_int * vy_m_uni;
    double vz_alt_B = vz_B - fat_int * vz_m_uni;
    
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_bgk ( f_R, vx_alt_B, vy_alt_B, vz_alt_B, rho_R, gx_R, gy_R, gz_R, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_alt_R, vy_alt_R, vz_alt_R, rho_B, gx_B, gy_B, gz_B, tau_m, op_col_BR, lattice );

    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, gx_R, gy_R, gz_R, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, gx_B, gy_B, gz_B, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
}

//====================================================================================================================//



//=============================== Collision step using the Santos, Facin & Philippi model (immicisble) ===============//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_SFP_recoll ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, double gx_B, 
						double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int = parameters.A_fact;
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];
    
    double f[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;

    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );
   
    double vx_alt_R = vx_R + fat_int * vx_m_uni;
    double vy_alt_R = vy_R + fat_int * vy_m_uni;
    double vz_alt_R = vz_R + fat_int * vz_m_uni;

    double vx_alt_B = vx_B - fat_int * vx_m_uni;
    double vy_alt_B = vy_B - fat_int * vy_m_uni;
    double vz_alt_B = vz_B - fat_int * vz_m_uni;
    
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_bgk ( f_R, vx_alt_B, vy_alt_B, vz_alt_B, rho_R, gx_R, gy_R, gz_R, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_alt_R, vy_alt_R, vz_alt_R, rho_B, gx_B, gy_B, gz_B, tau_m, op_col_BR, lattice );

    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, gx_R, gy_R, gz_R, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, gx_B, gy_B, gz_B, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
        
        f[i] = f_R[i] + f_B[i];
    }
    
     //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {vx_m, vy_m, vz_m};   
    
    double mod_grad = sqrt( vx_m * vx_m + vy_m * vy_m + vz_m * vz_m );  

	double beta = parameters.recoll;
	
    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
    
    //----------------------------------------------------------------------------------------------//
}

//====================================================================================================================//




//=============================== Collision step using the Santos model (immicisble) =================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, double gx_B, 
						double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;
    
    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );   

	double acc_B_Rx = gx_B + fat_int_RB * vx_m_uni + fat_int_BB * vx_m_uni;
    double acc_B_Ry = gy_B + fat_int_RB * vy_m_uni + fat_int_BB * vy_m_uni;
    double acc_B_Rz = gz_B + fat_int_RB * vz_m_uni + fat_int_BB * vz_m_uni;
    
    double acc_R_Bx = gx_R - fat_int_RB * vx_m_uni - fat_int_RR * vx_m_uni;
    double acc_R_By = gy_R - fat_int_RB * vy_m_uni - fat_int_RR * vy_m_uni;
    double acc_R_Bz = gz_R - fat_int_RB * vz_m_uni - fat_int_RR * vz_m_uni;    
    
    //------------ Colisão monofásica e bifásica --------------------------------------------------//
    
    op_bgk ( f_R, vx_B, vy_B, vz_B, rho_R, acc_R_Bx, acc_R_By, acc_R_Bz, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_R, vy_R, vz_R, rho_B, acc_B_Rx, acc_B_Ry, acc_B_Rz, tau_m, op_col_BR, lattice );
    
    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, gx_R, gy_R, gz_R, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, gx_B, gy_B, gz_B, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
}

//====================================================================================================================//




//=============================== Colisao do modelo de Santos imiscivel - forma do DSFD2026 ==========================//
//
//  Mesmo modelo de coll_immsantos_BGK, escrito exatamente como nos slides DSFD2026:
//
//      Omega_i^{rb} = (rho^b/rho) [ f^eq_i( rho^r , u^b ) - f^r_i ] / tau^m  +  (rho^b/rho) S^r_i
//      Omega_i^{br} = (rho^r/rho) [ f^eq_i( rho^b , u^r ) - f^b_i ] / tau^m  +  (rho^r/rho) S^b_i
//
//      S^r_i = G_i(  lambda rho^r n , u^b )        S^b_i = G_i( -lambda rho^b n , u^r )
//
//      G_i(F,v) = w_i ( 1 - 1/(2 tau^m) ) [ (c_i - v)/cs^2 + (c_i.v) c_i / cs^4 ] . F
//
//  ESTA FUNCAO E APENAS A FORMA DOS SLIDES.  NAO E A RECOMENDADA -- ver abaixo.
//
//  DIFERENCA em relacao a coll_immsantos_BGK: aquela passa a aceleracao por op_bgk(), que aplica o
//  deslocamento de meio passo de Guo  v -> v + a/2  tanto ao equilibrio quanto a fonte.  Como no
//  termo cruzado a velocidade e da OUTRA especie e a forca e desta, o deslocamento faz duas coisas:
//
//    1) injeta  rho^r a / (2 tau^m)  de momento a mais, o que equivale a  lambda / (1 - 1/(2 tau^m))
//       no balanco de anti-difusao -- ou seja, afina a interface;
//
//    2) o equilibrio deslocado e quadratico na velocidade, entao carrega tambem um termo  a a / 4 ,
//       isto e uma tensao  ( lambda^2 / 4 tau^m ) rho phi (1-phi) n n  localizada na interface.
//       Uma tensao com essa forma E tensao interfacial ( tipo Korteweg ), e ela chega pelo
//       equilibrio, que a rede trata com isotropia de segunda ordem -- e nao pela fonte discreta,
//       que e a origem classica das correntes espurias.
//
//  Medido ( gota estatica 2D, tau_R = tau_B = 1, tau_m = 1.5, lei de Laplace verificada ):
//
//      sigma = 0.97 lambda   com coll_immsantos_BGK        sigma = 0.34 lambda   com esta funcao
//
//  e, comparando A MESMA tensao interfacial, coll_immsantos_BGK da correntes espurias de 1.5 a 3.6
//  vezes MENORES, com a vantagem crescendo com sigma.  Ou seja: o deslocamento nao e um defeito, e
//  o que da ao modelo tensao interfacial barata.  Use coll_immsantos_BGK.
//
//  Esta funcao serve para: reproduzir literalmente o que esta nos slides, e isolar a contribuicao
//  do termo G_i sozinho.  Ela reproduz  |grad phi| = (lambda/cs^2) phi (1-phi)  com erro de 0.3% a
//  0.8% para qualquer tau^m, enquanto coll_immsantos_BGK afina a interface pelo fator do item (1).
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK_dsfd ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R,
						double gx_B, double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;

	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;

	//--------------- Calcula vel. e densidade das particulas ------------------------------------//

	double vx_R, vy_R, vz_R, rho_R;
	double vx_B, vy_B, vz_B, rho_B;

	calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
	calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

	double rho = rho_R + rho_B;
	double conc_R = rho_R / rho;
	double conc_B = 1.0 - conc_R;

	//------------- Direcao do gradiente, pelos mediadores ---------------------------------------//
	//
	//  soma_i M_i c_i ~ - delta cs^2 grad(phi)  =>  n = - unit( soma_i M_i c_i )

	double vx_m, vy_m, vz_m, rho_m;

	calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );

	double ex, ey, ez;

	unit_vector ( vx_m, vy_m, vz_m, ex, ey, ez );

	double n_x = - ex;
	double n_y = - ey;
	double n_z = - ez;

	//------------- Forcas de segregacao ---------------------------------------------------------//

	double lam_R = fat_int_RB + fat_int_RR;
	double lam_B = fat_int_RB + fat_int_BB;

	double Fx_R =   lam_R * rho_R * n_x;
	double Fy_R =   lam_R * rho_R * n_y;
	double Fz_R =   lam_R * rho_R * n_z;

	double Fx_B = - lam_B * rho_B * n_x;
	double Fy_B = - lam_B * rho_B * n_y;
	double Fz_B = - lam_B * rho_B * n_z;

	//------------- Equilibrios ------------------------------------------------------------------//

	double f_eq_RR[nvel], f_eq_BB[nvel], f_eq_RB[nvel], f_eq_BR[nvel];

	//  Auto-colisao: a aceleracao de corpo e desta especie, entao o deslocamento de Guo se aplica.

	double vxR_g = vx_R + 0.5 * gx_R;
	double vyR_g = vy_R + 0.5 * gy_R;
	double vzR_g = vz_R + 0.5 * gz_R;

	double vxB_g = vx_B + 0.5 * gx_B;
	double vyB_g = vy_B + 0.5 * gy_B;
	double vzB_g = vz_B + 0.5 * gz_B;

	dist_eq ( f_eq_RR, vxR_g, vyR_g, vzR_g, rho_R, lattice );
	dist_eq ( f_eq_BB, vxB_g, vyB_g, vzB_g, rho_B, lattice );

	//  Colisao cruzada: equilibrio construido sobre a velocidade da OUTRA especie, sem deslocamento.

	dist_eq ( f_eq_RB, vx_B, vy_B, vz_B, rho_R, lattice );
	dist_eq ( f_eq_BR, vx_R, vy_R, vz_R, rho_B, lattice );

	//------------- Termos de fonte --------------------------------------------------------------//

	double S_R[nvel]{};
	double S_B[nvel]{};

	double G_R[nvel]{};
	double G_B[nvel]{};

	source ( Fx_R, Fy_R, Fz_R, vx_B, vy_B, vz_B, rho_R, tau_m, S_R, lattice );
	source ( Fx_B, Fy_B, Fz_B, vx_R, vy_R, vz_R, rho_B, tau_m, S_B, lattice );

	if ( gx_R != 0 || gy_R != 0 || gz_R != 0 )
	{
		source ( gx_R * rho_R, gy_R * rho_R, gz_R * rho_R, vxR_g, vyR_g, vzR_g, rho_R, tau_R, G_R, lattice );
	}

	if ( gx_B != 0 || gy_B != 0 || gz_B != 0 )
	{
		source ( gx_B * rho_B, gy_B * rho_B, gz_B * rho_B, vxB_g, vyB_g, vzB_g, rho_B, tau_B, G_B, lattice );
	}

	//------------- Colisao ----------------------------------------------------------------------//

	double one_over_tau_R = 1.0 / tau_R;
	double one_over_tau_B = 1.0 / tau_B;
	double one_over_tau_m = 1.0 / tau_m;

	for ( int i = 0; i < nvel; i++ )
	{
		double op_RR = ( f_eq_RR[i] - f_R[i] ) * one_over_tau_R + G_R[i];
		double op_BB = ( f_eq_BB[i] - f_B[i] ) * one_over_tau_B + G_B[i];

		double op_RB = ( f_eq_RB[i] - f_R[i] ) * one_over_tau_m + S_R[i];
		double op_BR = ( f_eq_BR[i] - f_B[i] ) * one_over_tau_m + S_B[i];

		f_R[i] = f_R[i] + conc_R * op_RR + conc_B * op_RB;
		f_B[i] = f_B[i] + conc_B * op_BB + conc_R * op_BR;
	}
}

//====================================================================================================================//




//=============================== Collision step using the Santos model (immicisble) =================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK_Guo ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, double gx_B, 
						double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	//double fat_int_RR = parameters.A_fact_R;
	//double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	    
    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;
    
    //------------- Calcula a direção do gradiente ------------------------------------------------//

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );   
    
    //------------- Calcula o termo de força do Guo -----------------------------------------------//
    
    double G[nvel];
    
    double Fx =  fat_int_RB * vx_m_uni * rho_R * rho_B / rho;
    double Fy =  fat_int_RB * vy_m_uni * rho_R * rho_B / rho;
    double Fz =  fat_int_RB * vz_m_uni * rho_R * rho_B / rho;
    
    double ux = conc_R * vx_R + conc_B * vx_B;
    double uy = conc_R * vy_R + conc_B * vy_B;
    double uz = conc_R * vz_R + conc_B * vz_B;
    
    force_Guo ( Fx, Fy, Fz, ux, uy, uz, tau_m, G, lattice );
    
    //------------ Colisão monofásica e bifásica --------------------------------------------------//
    
    double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];
    
    op_bgk ( f_R, vx_B, vy_B, vz_B, rho_R, gx_R, gy_R, gz_R, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_R, vy_R, vz_R, rho_B, gx_B, gy_B, gz_B, tau_m, op_col_BR, lattice );
    
    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, gx_R, gy_R, gz_R, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, gx_B, gy_B, gz_B, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {        
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i] - G[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i] + G[i];
    }
}

//====================================================================================================================//




//=============================== Collision step using the Santos model (immicisble) =================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK_recoll ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, 
						double gx_B, double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];
    
    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;
    
	//------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );   

	double acc_B_Rx = gx_B + fat_int_RB * vx_m_uni;
    double acc_B_Ry = gy_B + fat_int_RB * vy_m_uni;
    double acc_B_Rz = gz_B + fat_int_RB * vz_m_uni;
    
    double acc_R_Bx = gx_R - fat_int_RB * vx_m_uni;
    double acc_R_By = gy_R - fat_int_RB * vy_m_uni;
    double acc_R_Bz = gz_R - fat_int_RB * vz_m_uni;    
    
    double acc_B_Bx = gx_B + fat_int_BB * vx_m_uni;
    double acc_B_By = gy_B + fat_int_BB * vy_m_uni;
    double acc_B_Bz = gz_B + fat_int_BB * vz_m_uni;
    
    double acc_R_Rx = gx_R - fat_int_RR * vx_m_uni;
    double acc_R_Ry = gy_R - fat_int_RR * vy_m_uni;
    double acc_R_Rz = gz_R - fat_int_RR * vz_m_uni;    
    
    //------------ Colisão monofásica e bifásica --------------------------------------------------//
    
    op_bgk ( f_R, vx_B, vy_B, vz_B, rho_R, acc_R_Bx, acc_R_By, acc_R_Bz, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_R, vy_R, vz_R, rho_B, acc_B_Rx, acc_B_Ry, acc_B_Rz, tau_m, op_col_BR, lattice );
    
    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, acc_R_Rx, acc_R_Ry, acc_R_Rz, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, acc_B_Bx, acc_B_By, acc_B_Bz, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
    
    //--------------- Etapa de recoloração ----------------------------------------------//
    
    double f[nvel];
    
    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f_R[i] + f_B[i];
	}
	
	double mod_grad = sqrt( vx_m * vx_m + vy_m * vy_m + vz_m * vz_m);
    
    double grad[3] = {vx_m, vy_m, vz_m};    

	double beta = parameters.recoll;
	
    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//



//=============================== Collision step using the Santos model (immicisble with gradients ) =================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK ( double* f_R, double* f_B, int x, int y, int z, double gx_R, double gy_R, double gz_R, 
							double gx_B, double gy_B, double gz_B, GEOMETRY geometry, LATTICE lattice, 
							PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;
    
    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
	
    double acc_B_Rx = gx_B - fat_int_RB * grad_uni_x;
    double acc_B_Ry = gy_B - fat_int_RB * grad_uni_y;
    double acc_B_Rz = gz_B - fat_int_RB * grad_uni_z;
    
    double acc_R_Bx = gx_R + fat_int_RB * grad_uni_x;
    double acc_R_By = gy_R + fat_int_RB * grad_uni_y;
    double acc_R_Bz = gz_R + fat_int_RB * grad_uni_z;    
    
    double acc_B_Bx = gx_B - fat_int_BB * grad_uni_x;
    double acc_B_By = gy_B - fat_int_BB * grad_uni_y;
    double acc_B_Bz = gz_B - fat_int_BB * grad_uni_z;
    
    double acc_R_Rx = gx_R + fat_int_RR * grad_uni_x;
    double acc_R_Ry = gy_R + fat_int_RR * grad_uni_y;
    double acc_R_Rz = gz_R + fat_int_RR * grad_uni_z;  
       
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_bgk ( f_R, vx_B, vy_B, vz_B, rho_R, acc_R_Bx, acc_R_By, acc_R_Bz, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_R, vy_R, vz_R, rho_B, acc_B_Rx, acc_B_Ry, acc_B_Rz, tau_m, op_col_BR, lattice );
    
    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, acc_R_Rx, acc_R_Ry, acc_R_Rz, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, acc_B_Bx, acc_B_By, acc_B_Bz, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
}

//====================================================================================================================//







//=============================== Collision step using the Santos model (immicisble with gradients ) =================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK_LG ( double* f_R, double* f_B, int x, int y, int z, double gx_R, double gy_R, double gz_R, 
							double gx_B, double gy_B, double gz_B, GEOMETRY geometry, LATTICE lattice, 
							PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;
    
    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
	
    double acc_Bx = gx_B - fat_int_RB * grad_uni_x;
    double acc_By = gy_B - fat_int_RB * grad_uni_y;
    double acc_Bz = gz_B - fat_int_RB * grad_uni_z;
    
    double acc_Rx = gx_R + fat_int_RB * grad_uni_x;
    double acc_Ry = gy_R + fat_int_RB * grad_uni_y;
    double acc_Rz = gz_R + fat_int_RB * grad_uni_z;   
    
    //------------ Interação interna do fluido Red-Red --------------------------------------------//
    
    double rho_std = 1.;
    
    double alpha_2 = ( ( rho - rho_std ) / rho_std ) * ( ( rho - rho_std ) / rho_std );
    
    double fact_dens = exp( - 10 * alpha_2 );
    
    double acc_B_Bx = gx_B - fat_int_BB * grad_uni_x;
    double acc_B_By = gy_B - fat_int_BB * grad_uni_y;
    double acc_B_Bz = gz_B - fat_int_BB * grad_uni_z;
    
    double acc_R_Rx = gx_R + fat_int_RR * grad_uni_x + fat_int_RR * fact_dens * grad_uni_x; 
    double acc_R_Ry = gy_R + fat_int_RR * grad_uni_y + fat_int_RR * fact_dens * grad_uni_y;
    double acc_R_Rz = gz_R + fat_int_RR * grad_uni_z + fat_int_RR * fact_dens * grad_uni_z;  
    
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_bgk ( f_R, vx_B, vy_B, vz_B, rho_R, acc_Rx, acc_Ry, acc_Rz, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_R, vy_R, vz_R, rho_B, acc_Bx, acc_By, acc_Bz, tau_m, op_col_BR, lattice );
    
    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, acc_R_Rx, acc_R_Ry, acc_R_Rz, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, acc_B_Bx, acc_B_By, acc_B_Bz, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
}

//====================================================================================================================//




//=============================== Collision step using the Santos model (immicisble with gradients ) =================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_BGK_LG ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, double gx_B, 
						double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;
    
    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );   

	double acc_B_Rx = gx_B + fat_int_RB * vx_m_uni + fat_int_BB * vx_m_uni;
    double acc_B_Ry = gy_B + fat_int_RB * vy_m_uni + fat_int_BB * vy_m_uni;
    double acc_B_Rz = gz_B + fat_int_RB * vz_m_uni + fat_int_BB * vz_m_uni;
    
    double acc_R_Bx = gx_R - fat_int_RB * vx_m_uni - fat_int_RR * vx_m_uni;
    double acc_R_By = gy_R - fat_int_RB * vy_m_uni - fat_int_RR * vy_m_uni;
    double acc_R_Bz = gz_R - fat_int_RB * vz_m_uni - fat_int_RR * vz_m_uni;  
    
    //------------ Interação interna do fluido Red-Red --------------------------------------------//
    
    //double rho_std = 100.;
    
    //double alpha_2 = ( ( rho - rho_std ) / rho_std ) * ( ( rho - rho_std ) / rho_std );
    
    //double fact_dens = exp( - 1. * alpha_2 );
    double fact_dens = 1.;
    
    if ( rho > 50 ) fact_dens = 0.;
    
    double acc_R_Rx = gx_R - fat_int_RR * fact_dens * vx_m_uni;
    double acc_R_Ry = gy_R - fat_int_RR * fact_dens * vy_m_uni;
    double acc_R_Rz = gz_R - fat_int_RR * fact_dens * vz_m_uni;  
    
    //------------ Colisão monofásica e bifásica --------------------------------------------------//
    
    op_bgk ( f_R, vx_B, vy_B, vz_B, rho_R, acc_R_Bx, acc_R_By, acc_R_Bz, tau_m, op_col_RB, lattice );
    op_bgk ( f_B, vx_R, vy_R, vz_R, rho_B, acc_B_Rx, acc_B_Ry, acc_B_Rz, tau_m, op_col_BR, lattice );
    
    op_bgk ( f_R, vx_R, vy_R, vz_R, rho_R, acc_R_Rx, acc_R_Ry, acc_R_Rz, tau_R, op_col_R, lattice );
    op_bgk ( f_B, vx_B, vy_B, vz_B, rho_B, gx_B, gy_B, gz_B, tau_B, op_col_B, lattice );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
}

//====================================================================================================================//




//=============================== Collision step using the Santos model (immicisble) =================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_phasetrans_BGK ( double* f, double* f_m, double gx, double gy, double gz, LATTICE lattice, 
							PARAMETERS parameters )
{
	double fat_int = parameters.A_fact_R;
	
	double tau = parameters.tau;
	
	double op_col[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx, vy, vz, rho;
 
    calcula ( f, vx, vy, vz, rho, lattice );

    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );  
    
    double rho_std = 1.;
    
    double alpha_2 = ( ( rho - rho_std ) / rho_std ) * ( ( rho - rho_std ) / rho_std );
    
    double fact_dens = exp( - 10 * alpha_2 );
   
    double acc_x = gx - fat_int * fact_dens * vx_m_uni;
    double acc_y = gy - fat_int * fact_dens * vy_m_uni;
    double acc_z = gz - fat_int * fact_dens * vz_m_uni;
    
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_bgk ( f, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau, op_col, lattice );   

    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f[i] + op_col[i];
    }
}

//====================================================================================================================//




//=============================== Collision step using the Santos model (phase transitions) ==========================//
//
//      Input: position, acelerations, geometry, lattice, parameters
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_phasetrans_BGK ( double *f, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry, 
							LATTICE lattice, PARAMETERS parameters )
{
	
	double fat_int = parameters.A_fact_R;
	
	double tau = parameters.tau;
	
	double op_col[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx, vy, vz, rho;
 
    calcula ( f, vx, vy, vz, rho, lattice );

    //------------- Calcula a direção do gradiente ------------------------------------------------//    
    
    double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );  
    
    double rho_std = 1.;
    
    double alpha_2 = ( ( rho - rho_std ) / rho_std ) * ( ( rho - rho_std ) / rho_std );
    
    double fact_dens = exp( - 10 * alpha_2 );
   
    double acc_x = gx + fat_int * fact_dens * grad_uni_x;
    double acc_y = gy + fat_int * fact_dens * grad_uni_y;
    double acc_z = gz + fat_int * fact_dens * grad_uni_z;
    
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_bgk ( f, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau, op_col, lattice );   

    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f[i] + op_col[i];
    }
}

//====================================================================================================================//



//=============================== Collision step using the Shan Chen model (phase transitions) =======================//
//
//      Input: position, acelerations, geometry, lattice, parameters
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_phasetrans_BGK_SC ( double *f, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry, 
							LATTICE lattice, PARAMETERS parameters )
{
	double G = parameters.A_fact_R;
	
	double tau = parameters.tau;
	
	double op_col[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx, vy, vz, rho;
 
    calcula ( f, vx, vy, vz, rho, lattice );

    //------------- Calcula a direção do gradiente ------------------------------------------------//    
    
    double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    double *psi = ini_psi + x + y * geometry.nx + z * geometry.nx * geometry.ny;
       
    double acc_x = gx - G * psi[0] * grad_x;
    double acc_y = gy - G * psi[0] * grad_y;
    double acc_z = gz - G * psi[0] * grad_z;
    
    //------------ Colisão monofásica -------------------------------------------------------------//
    
    op_bgk ( f, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau, op_col, lattice );   

    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f[i] + op_col[i];
    }
}

//====================================================================================================================//



//=============================== Collision step using the Santos model (immicisble) =================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_immsantos_ExtReg ( double* f_R, double* f_B, double* f_m, double gx_R, double gy_R, double gz_R, double gx_B, 
						double gy_B, double gz_B, LATTICE lattice, PARAMETERS parameters )
{
	double fat_int_RB = parameters.A_fact;
	double fat_int_RR = parameters.A_fact_R;
	double fat_int_BB = parameters.A_fact_B;	
	
	double tau_m = parameters.tau_m;
	double tau_R = parameters.tau_R;
	double tau_B = parameters.tau_B;
	
	double op_col_R[nvel];
    double op_col_B[nvel];

    double op_col_RB[nvel];
    double op_col_BR[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;
    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    double rho = rho_R + rho_B;
    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;

    //------------- Calcula a direção do gradiente ------------------------------------------------/

    double vx_m, vy_m, vz_m, rho_m;
    
    calcula ( f_m, vx_m, vy_m, vz_m, rho_m, lattice );
    
    double vx_m_uni, vy_m_uni, vz_m_uni;

    unit_vector ( vx_m, vy_m, vz_m, vx_m_uni, vy_m_uni, vz_m_uni );   

	double acc_B_Rx = gx_B + fat_int_RB * vx_m_uni;
    double acc_B_Ry = gy_B + fat_int_RB * vy_m_uni;
    double acc_B_Rz = gz_B + fat_int_RB * vz_m_uni;
    
    double acc_R_Bx = gx_R - fat_int_RB * vx_m_uni;
    double acc_R_By = gy_R - fat_int_RB * vy_m_uni;
    double acc_R_Bz = gz_R - fat_int_RB * vz_m_uni;
    
    double acc_B_Bx = gx_B + fat_int_BB * vx_m_uni;
    double acc_B_By = gy_B + fat_int_BB * vy_m_uni;
    double acc_B_Bz = gz_B + fat_int_BB * vz_m_uni;
    
    double acc_R_Rx = gx_R - fat_int_RR * vx_m_uni;
    double acc_R_Ry = gy_R - fat_int_RR * vy_m_uni;
    double acc_R_Rz = gz_R - fat_int_RR * vz_m_uni;      
 
    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    op_Ext_Reg ( f_R, vx_B, vy_B, vz_B, rho_R, acc_R_Bx, acc_R_By, acc_R_Bz, tau_m, op_col_RB, lattice, parameters );
    op_Ext_Reg ( f_B, vx_R, vy_R, vz_R, rho_B, acc_B_Rx, acc_B_Ry, acc_B_Rz, tau_m, op_col_BR, lattice, parameters );

    op_Ext_Reg ( f_R, vx_R, vy_R, vz_R, rho_R, acc_R_Rx, acc_R_Ry, acc_R_Rz, tau_R, op_col_R, lattice, parameters );
    op_Ext_Reg ( f_B, vx_B, vy_B, vz_B, rho_B, acc_B_Bx, acc_B_By, acc_B_Bz, tau_B, op_col_B, lattice, parameters );

    for ( int i = 0; i < nvel; i++ )
    {
        f_R[i] = f_R[i] + conc_R * op_col_R[i] + conc_B * op_col_RB[i];
        f_B[i] = f_B[i] + conc_B * op_col_B[i] + conc_R * op_col_BR[i];
    }
}

//====================================================================================================================//



//=============================== Modelo color-gradient 3D melhorado - Wen, Li, Yu & Luo (2019) ======================//
//
//  Z. X. Wen, Q. Li, Y. Yu e Kai H. Luo, "Improved three-dimensional color-gradient lattice Boltzmann
//  model for immiscible two-phase flows", Phys. Rev. E 100, 023301 (2019).
//
//  Modelo de gradiente de cor com tres operadores, colisao MRT:
//
//      Omega_i = ( Omega_i )^(3) [ ( Omega_i )^(1) + ( Omega_i )^(2) ]
//
//      (1) colisao monofasica MRT com termo de correcao   -- eq. (44)
//      (2) operador de perturbacao ( tensao interfacial ) -- eqs. (14) e (60)
//      (3) recoloracao de Latva-Kokko                     -- eqs. (9) e (10)
//
//  O QUE E "MELHORADO":  nos modelos anteriores a pressao entra no equilibrio como p_k = rho_k cs_k^2
//  no SEGUNDO momento, mas o TERCEIRO momento continua carregando rho_k c^2/3.  Essa inconsistencia
//  gera termos de erro na equacao de momento -- eq. (17) do artigo -- que quebram a invariancia
//  galileana quando as densidades sao diferentes.  O modelo corrige isso de duas maneiras:
//
//    a) um termo de ordem alta no equilibrio, eq. (22), que poe p_k tambem nos elementos FORA da
//       diagonal do terceiro momento;
//
//    b) um termo de correcao G_i na colisao, eqs. (24) e (52), que trata os elementos DA diagonal,
//       que a simetria da D3Q19 nao permite corrigir pelo equilibrio.
//
//  Equilibrio ( eq. 22, com c = 1 ):
//
//      f^eq_i = rho_k [ phi^k_i + w_i ( 3 (e.u) + 9/2 (e.u)^2 - 3/2 |u|^2
//                                       + 3 (e.u) ( 3 cs_k^2 - 1 ) ( 3 |e_i|^2 - 5 ) ) ]
//
//      phi^k_0 = alpha_k ,  phi^k_{1..6} = (1-alpha_k)/12 ,  phi^k_{7..18} = (1-alpha_k)/24
//      cs_k^2  = ( 1 - alpha_k ) / 2       p_k = rho_k cs_k^2
//
//  NOTA sobre o coeficiente do termo de ordem alta:  o artigo o imprime como 3(e.u)/(2c^2)(...)(...),
//  mas esse valor nao satisfaz a eq. (23), que e o proprio objetivo do termo.  Impondo a eq. (23) --
//  terceiro momento fora da diagonal igual a p_k, e nao a rho_k c^2/3 -- o coeficiente sai unico e
//  vale 3, nao 3/2.  Com 3 as quinze condicoes das eqs. (20) e (23) sao satisfeitas exatamente
//  ( verificado em algebra simbolica ); com 3/2, seis delas falham.  Ver LI_YU_LUO.md.
//
//  Correcao G_i, em espaco de momentos ( eq. 52 ):  so os momentos 4, 5 e 6 sao nao nulos,
//
//      C_4 = Qx + Qy + Qz ,   C_5 = 2 Qx - Qy - Qz ,   C_6 = Qy - Qz
//
//      Q_alpha = d_alpha [ rho_k u_alpha ( c^2 - 3 cs_k^2 ) ]        eqs. (39)-(41)
//
//  Q_alpha e calculado fora desta funcao ( precisa dos vizinhos ) e chega pelos vetores Q_R e Q_B.
//  Repare que  c^2 - 3 cs_k^2 = 0  quando alpha_k = 1/3: com densidades iguais a correcao some, que
//  e exatamente o que o artigo diz na secao II C.
//
//  Perturbacao ( eq. 14 ), com n = grad(rho^N) / |grad(rho^N)| :
//
//      ( Omega^k_i )^(2) = ( A_k / 2 ) |grad rho^N| [ w_i (e_i.n)^2 - B_i ]
//      B_0 = -1/3 ,  B_{1..6} = 1/18 ,  B_{7..18} = 1/36
//
//  Em espaco de momentos esse vetor tem forma fechada -- nao e preciso multiplicar por M:
//
//      mP_4 = -4/9                mP_7 = 2 nx ny / 9        mP_16 = - nz^2 / 9
//      mP_5 = 2 ( 3 nx^2 - 1 )/9  mP_8 = 2 nx nz / 9        mP_17 = - ny^2 / 9
//      mP_6 = 2 ( ny^2 - nz^2 )/9 mP_9 = 2 ny nz / 9        mP_18 = - nx^2 / 9
//
//  e todos os outros sao zero -- em particular mP_0 = mP_1 = mP_2 = mP_3 = 0, ou seja a perturbacao
//  conserva massa e quantidade de movimento exatamente.  Seguindo a eq. (60), ela e multiplicada
//  pela matriz de relaxacao, o que torna a tensao interfacial independente de tau:
//
//      sigma = 2 ( A_R + A_B ) c^4 dt / 9
//
//  Recoloracao de Latva-Kokko ( eqs. 9 e 10 ):
//
//      f^{R,+}_i = (rho_R/rho) f*_i + beta (rho_R rho_B / rho^2) cos(phi_i) soma_k rho_k phi^k_i
//      f^{B,+}_i = (rho_B/rho) f*_i - beta (rho_R rho_B / rho^2) cos(phi_i) soma_k rho_k phi^k_i
//
//      cos(phi_i) = ( e_i . grad rho^N ) / ( |e_i| |grad rho^N| )
//
//  Parametros, lidos de PARAMETERS:
//
//      tau_R , tau_B    tau_v de cada fluido puro     ( nu_k = cs_k^2 ( tau_v^k - 1/2 ) )
//      tau_m            tau_e = tau_q = tau_pi        ( o artigo usa 1.0 )
//      A_fact_R/_B      A_R e A_B da tensao interfacial
//      recoll           beta da recoloracao
//      rho_ini_R/_B     rho_R^in e rho_B^in           ( definicao de rho^N, eq. 12 )
//      alfa_R           alpha_R;  alpha_B sai da condicao de equilibrio de pressao
//
//====================================================================================================================//

//  As matrizes MRT e a colisao de Li-Yu-Luo sao especificas do D3Q19.  Num programa D3Q27
//  ( Saito ) elas nao existem -- e indexar 19x19 com nvel = 27 seria erro de memoria.

#if nvel == 19

static const double MRT_LYL_M[19][19] =
{
	{        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1. },
	{        0.,        1.,       -1.,        0.,        0.,        0.,        0.,        1.,       -1.,        1.,       -1.,        1.,       -1.,        1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        1.,       -1.,        0.,        0.,        1.,       -1.,       -1.,        1.,        0.,        0.,        0.,        0.,       -1.,        1.,       -1.,        1. },
	{        0.,        0.,        0.,        0.,        0.,        1.,       -1.,        0.,        0.,        0.,        0.,        1.,       -1.,       -1.,        1.,       -1.,        1.,        1.,       -1. },
	{        0.,        1.,        1.,        1.,        1.,        1.,        1.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2. },
	{        0.,        2.,        2.,       -1.,       -1.,       -1.,       -1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,       -2.,       -2.,       -2.,       -2. },
	{        0.,        0.,        0.,        1.,        1.,       -1.,       -1.,        1.,        1.,        1.,        1.,       -1.,       -1.,       -1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,       -1.,       -1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,       -1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,       -1.,       -1. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,       -1.,        1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,        1.,       -1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,       -1.,        1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,        1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,       -1.,        1.,        1.,       -1. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,       -1.,        1.,       -1.,        1. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1.,        1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1.,        1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1.,        1. },
};

static const double MRT_LYL_M_INV[19][19] =
{
	{        1.,        0.,        0.,        0.,       -1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1. },
	{        0.,     1./2.,        0.,        0.,     1./6.,     1./6.,        0.,        0.,        0.,        0.,        0.,    -1./2.,        0.,    -1./2.,        0.,        0.,    -1./2.,    -1./2.,        0. },
	{        0.,    -1./2.,        0.,        0.,     1./6.,     1./6.,        0.,        0.,        0.,        0.,        0.,     1./2.,        0.,     1./2.,        0.,        0.,    -1./2.,    -1./2.,        0. },
	{        0.,        0.,     1./2.,        0.,     1./6.,   -1./12.,     1./4.,        0.,        0.,        0.,    -1./2.,        0.,        0.,        0.,        0.,    -1./2.,    -1./2.,        0.,    -1./2. },
	{        0.,        0.,    -1./2.,        0.,     1./6.,   -1./12.,     1./4.,        0.,        0.,        0.,     1./2.,        0.,        0.,        0.,        0.,     1./2.,    -1./2.,        0.,    -1./2. },
	{        0.,        0.,        0.,     1./2.,     1./6.,   -1./12.,    -1./4.,        0.,        0.,        0.,        0.,        0.,    -1./2.,        0.,    -1./2.,        0.,        0.,    -1./2.,    -1./2. },
	{        0.,        0.,        0.,    -1./2.,     1./6.,   -1./12.,    -1./4.,        0.,        0.,        0.,        0.,        0.,     1./2.,        0.,     1./2.,        0.,        0.,    -1./2.,    -1./2. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,     1./4.,     1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,    -1./4.,    -1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,    -1./4.,     1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,     1./4.,    -1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,     1./4.,     1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,    -1./4.,    -1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,    -1./4.,     1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,     1./4.,    -1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,        0.,    -1./4.,    -1./4.,        0.,        0.,     1./4. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,        0.,     1./4.,     1./4.,        0.,        0.,     1./4. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,        0.,     1./4.,    -1./4.,        0.,        0.,     1./4. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,        0.,    -1./4.,     1./4.,        0.,        0.,     1./4. },
};

#pragma acc declare copyin( MRT_LYL_M, MRT_LYL_M_INV )

#endif   // nvel == 19

//------------------ Equilibrio do modelo, eq. (22) -----------------------------------------------------------------//

#pragma acc routine seq
void dist_eq_LYL ( double *f_eq, double vx, double vy, double vz, double rho_k, double alfa_k,
                   double hi_order, LATTICE lattice )
{
	const double cs2_k = 0.5 * ( 1.0 - alfa_k );

	//  hi_order = 1 liga o termo de ordem alta da eq. (22);  0 recai no equilibrio da eq. (6),
	//  que e o "modelo original" com que o artigo compara.

	const double fator = hi_order * 3.0 * ( 3.0 * cs2_k - 1.0 );

	const double u2 = vx * vx + vy * vy + vz * vz;

	for ( int i = 0; i < nvel; i++ )
	{
		const double cx = lattice.c_i[ i * dim + 0 ];
		const double cy = lattice.c_i[ i * dim + 1 ];
		const double cz = ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0;

		const double eu = cx * vx + cy * vy + cz * vz;
		const double e2 = cx * cx + cy * cy + cz * cz;

		//  phi^k_i : alpha_k no repouso, (1-alpha_k)/12 nos eixos, (1-alpha_k)/24 nas diagonais

		double phi_i;

		if      ( i == 0 ) phi_i = alfa_k;
		else if ( i <  7 ) phi_i = ( 1.0 - alfa_k ) / 12.0;
		else               phi_i = ( 1.0 - alfa_k ) / 24.0;

		f_eq[i] = rho_k * ( phi_i + lattice.w[i] * ( 3.0 * eu + 4.5 * eu * eu - 1.5 * u2
		                                             + fator * eu * ( 3.0 * e2 - 5.0 ) ) );
	}
}

//------------------ Momentos de equilibrio, m^eq = M f^eq ----------------------------------------------------------//
//
//  Forma fechada, obtida simbolicamente de M f^eq e conferida contra ela em tempo de execucao por
//  check_LYL().  ATENCAO: m4 e m16..m18 NAO coincidem com a lista impressa no apendice do artigo --
//  m4 = 3 p_k + rho_k |u|^2 e imposto pela propria eq. (20) do artigo, e m16..m18 sao momentos de
//  quarta ordem que nao afetam a hidrodinamica.  Ver LI_YU_LUO.md.

#pragma acc routine seq
void mom_eq_LYL ( double *m_eq, double vx, double vy, double vz, double rho_k, double alfa_k )
{
	const double cs2_k = 0.5 * ( 1.0 - alfa_k );

	const double p_k = rho_k * cs2_k;

	const double ux2 = vx * vx, uy2 = vy * vy, uz2 = vz * vz;

	const double u2 = ux2 + uy2 + uz2;

	m_eq[ 0] = rho_k;
	m_eq[ 1] = rho_k * vx;
	m_eq[ 2] = rho_k * vy;
	m_eq[ 3] = rho_k * vz;

	m_eq[ 4] = 3.0 * p_k + rho_k * u2;
	m_eq[ 5] = rho_k * ( 2.0 * ux2 - uy2 - uz2 );
	m_eq[ 6] = rho_k * ( uy2 - uz2 );

	m_eq[ 7] = rho_k * vx * vy;
	m_eq[ 8] = rho_k * vx * vz;
	m_eq[ 9] = rho_k * vy * vz;

	m_eq[10] = p_k * vy;
	m_eq[11] = p_k * vx;
	m_eq[12] = p_k * vz;
	m_eq[13] = p_k * vx;
	m_eq[14] = p_k * vz;
	m_eq[15] = p_k * vy;

	m_eq[16] = rho_k * ( 1.0 - alfa_k + 2.0 * ux2 + 2.0 * uy2 - uz2 ) / 6.0;
	m_eq[17] = rho_k * ( 1.0 - alfa_k + 2.0 * ux2 - uy2 + 2.0 * uz2 ) / 6.0;
	m_eq[18] = rho_k * ( 1.0 - alfa_k - ux2 + 2.0 * uy2 + 2.0 * uz2 ) / 6.0;
}

//------------------ Tempo de relaxacao suavizado atraves da interface ( Grunau, eq. 59 ) ---------------------------//
//
//  Interpolacao parabolica de 1/tau entre os dois fluidos, com derivada nula nas emendas em
//  rho^N = +-delta, de modo que tau_v seja continuo e suave atraves da interface.

#pragma acc routine seq
double tau_v_LYL ( double rho_N, double tau_R, double tau_B, double delta )
{
	if ( rho_N >  delta ) return tau_R;
	if ( rho_N < -delta ) return tau_B;

	const double s_R = 1.0 / tau_R;
	const double s_B = 1.0 / tau_B;

	const double chi = 2.0 * s_R * s_B / ( s_R + s_B );

	double s;

	if ( rho_N > 0.0 )
	{
		const double eta = 2.0 * ( s_R - chi ) / delta;

		s = chi + eta * rho_N - eta * rho_N * rho_N / ( 2.0 * delta );
	}
	else
	{
		const double lam = 2.0 * ( chi - s_B ) / delta;

		s = chi + lam * rho_N + lam * rho_N * rho_N / ( 2.0 * delta );
	}

	return 1.0 / s;
}

//------------------ Colisao completa de um sitio -------------------------------------------------------------------//
//
//      grad_N   grad( rho^N ), calculado fora ( precisa dos vizinhos )
//      Q_R,Q_B  Q_alpha de cada fluido, eqs. (39)-(41), tambem calculados fora
//      gx,gy,gz aceleracao de corpo ( entra pelo esquema de Guo em espaco de momentos )

#if nvel == 19

#pragma acc routine seq
void coll_LYL_D3Q19 ( double *f_R, double *f_B, double *grad_N, double *Q_R, double *Q_B,
                      double gx, double gy, double gz, LATTICE lattice, PARAMETERS parameters )
{
	const double eps = 1.0e-30;

	//------------- Densidades, velocidade da mistura e rho^N ------------------------------------//

	double rho_R = 0.0, rho_B = 0.0;
	double jx = 0.0, jy = 0.0, jz = 0.0;

	for ( int i = 0; i < nvel; i++ )
	{
		rho_R = rho_R + f_R[i];
		rho_B = rho_B + f_B[i];

		const double s = f_R[i] + f_B[i];

		jx = jx + s * lattice.c_i[ i * dim + 0 ];
		jy = jy + s * lattice.c_i[ i * dim + 1 ];
		jz = jz + s * lattice.c_i[ i * dim + 2 ];
	}

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	//  Velocidade da mistura, com o meio passo de Guo quando ha forca de corpo.

	const double ux = ( jx + 0.5 * gx * rho ) / rho;
	const double uy = ( jy + 0.5 * gy * rho ) / rho;
	const double uz = ( jz + 0.5 * gz * rho ) / rho;

	const double alfa_R = parameters.alfa_R;
	const double alfa_B = parameters.alfa_B;

	const double a_R = rho_R / parameters.rho_ini_R;
	const double a_B = rho_B / parameters.rho_ini_B;

	const double rho_N = ( a_R + a_B > eps ) ? ( a_R - a_B ) / ( a_R + a_B ) : 0.0;

	//------------- Matriz de relaxacao ----------------------------------------------------------//

	const double tau_v = tau_v_LYL ( rho_N, parameters.tau_R, parameters.tau_B, parameters.delta_tau );

	const double tau_e = parameters.tau_m;
	const double tau_q = parameters.tau_m;
	const double tau_p = parameters.tau_m;

	double S[nvel];

	S[0] = S[1] = S[2] = S[3] = 1.0;
	S[4] = 1.0 / tau_e;
	S[5] = S[6] = S[7] = S[8] = S[9] = 1.0 / tau_v;
	for ( int i = 10; i < 16; i++ ) S[i] = 1.0 / tau_q;
	S[16] = S[17] = S[18] = 1.0 / tau_p;

	//------------- Direcao normal e modulo do gradiente de cor ----------------------------------//

	const double mod_grad = sqrt ( grad_N[0] * grad_N[0] + grad_N[1] * grad_N[1] + grad_N[2] * grad_N[2] );

	double n_x = 0.0, n_y = 0.0, n_z = 0.0;

	if ( mod_grad > eps )
	{
		n_x = grad_N[0] / mod_grad;
		n_y = grad_N[1] / mod_grad;
		n_z = grad_N[2] / mod_grad;
	}

	//  Perturbacao em espaco de momentos ( forma fechada; ver o cabecalho ).

	double mP[nvel];

	for ( int i = 0; i < nvel; i++ ) mP[i] = 0.0;

	if ( mod_grad > eps )
	{
		const double um_nono = 1.0 / 9.0;

		mP[ 4] = -4.0 * um_nono;
		mP[ 5] =  2.0 * ( 3.0 * n_x * n_x - 1.0 ) * um_nono;
		mP[ 6] =  2.0 * ( n_y * n_y - n_z * n_z ) * um_nono;
		mP[ 7] =  2.0 * n_x * n_y * um_nono;
		mP[ 8] =  2.0 * n_x * n_z * um_nono;
		mP[ 9] =  2.0 * n_y * n_z * um_nono;
		mP[16] = -n_z * n_z * um_nono;
		mP[17] = -n_y * n_y * um_nono;
		mP[18] = -n_x * n_x * um_nono;
	}

	//------------- Colisao de cada fluido em espaco de momentos ---------------------------------//

	double f_star[nvel];

	for ( int i = 0; i < nvel; i++ ) f_star[i] = 0.0;

	for ( int k = 0; k < 2; k++ )
	{
		double *f_k = ( k == 0 ) ? f_R : f_B;

		const double rho_k = ( k == 0 ) ? rho_R : rho_B;
		const double alfa_k = ( k == 0 ) ? alfa_R : alfa_B;
		const double A_k = ( k == 0 ) ? parameters.A_fact_R : parameters.A_fact_B;
		const double *Q_k = ( k == 0 ) ? Q_R : Q_B;

		//----- m = M f -------------------------------------------------------------------------//

		double m[nvel], m_eq[nvel];

		for ( int a = 0; a < nvel; a++ )
		{
			double s = 0.0;

			for ( int i = 0; i < nvel; i++ ) s = s + MRT_LYL_M[a][i] * f_k[i];

			m[a] = s;
		}

		if ( parameters.lyl_melhorado != 0.0 )
		{
			mom_eq_LYL ( m_eq, ux, uy, uz, rho_k, alfa_k );
		}
		else
		{
			//  Sem o termo de ordem alta a forma fechada nao vale: transforma o equilibrio da
			//  eq. (6) para o espaco de momentos.  So e usado nas comparacoes com o modelo antigo.

			double f_eq[nvel];

			dist_eq_LYL ( f_eq, ux, uy, uz, rho_k, alfa_k, 0.0, lattice );

			for ( int a = 0; a < nvel; a++ )
			{
				double t = 0.0;

				for ( int i = 0; i < nvel; i++ ) t = t + MRT_LYL_M[a][i] * f_eq[i];

				m_eq[a] = t;
			}
		}

		//----- Termo de correcao em espaco de momentos ( eq. 52 ) ------------------------------//

		double C[nvel];

		for ( int a = 0; a < nvel; a++ ) C[a] = 0.0;

		if ( parameters.lyl_melhorado != 0.0 )
		{
			C[4] = Q_k[0] + Q_k[1] + Q_k[2];
			C[5] = 2.0 * Q_k[0] - Q_k[1] - Q_k[2];
			C[6] = Q_k[1] - Q_k[2];
		}

		//----- Forca de corpo, esquema de Guo em espaco de momentos ----------------------------//
		//
		//  F = rho_k g ; os momentos de F_i = w_i [ (e-u)/cs^2 + (e.u) e/cs^4 ] . F valem, para os
		//  momentos que importam: 2 j.F no momento 4, e as componentes do tensor F u + u F nos
		//  momentos 5..9.  Os momentos 1..3 recebem F, mas ali s = 1 e o equilibrio ja usa u com o
		//  meio passo, entao a injecao correta ja esta contida em m_eq.

		double C_F[nvel];

		for ( int a = 0; a < nvel; a++ ) C_F[a] = 0.0;

		if ( gx != 0.0 || gy != 0.0 || gz != 0.0 )
		{
			const double Fx = rho_k * gx, Fy = rho_k * gy, Fz = rho_k * gz;

			C_F[ 4] = 2.0 * ( ux * Fx + uy * Fy + uz * Fz );
			C_F[ 5] = 2.0 * ( 2.0 * ux * Fx - uy * Fy - uz * Fz );
			C_F[ 6] = 2.0 * ( uy * Fy - uz * Fz );
			C_F[ 7] = ux * Fy + uy * Fx;
			C_F[ 8] = ux * Fz + uz * Fx;
			C_F[ 9] = uy * Fz + uz * Fy;
		}

		//----- Omega em espaco de momentos ------------------------------------------------------//
		//
		//      m_bar = - S ( m - m^eq )  +  ( I - S/2 ) ( C + C_F )  +  S ( A_k/2 ) |grad| mP
		//
		//  O ultimo termo e a eq. (60): a perturbacao passa pela matriz de relaxacao, o que torna
		//  a tensao interfacial independente de tau.

		const double pert = 0.5 * A_k * mod_grad;

		double m_bar[nvel];

		for ( int a = 0; a < nvel; a++ )
		{
			m_bar[a] = - S[a] * ( m[a] - m_eq[a] )
			           + ( 1.0 - 0.5 * S[a] ) * ( C[a] + C_F[a] )
			           + S[a] * pert * mP[a];
		}

		//----- Volta ao espaco das distribuicoes e acumula f* -----------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double s = 0.0;

			for ( int a = 0; a < nvel; a++ ) s = s + MRT_LYL_M_INV[i][a] * m_bar[a];

			f_k[i] = f_k[i] + s;

			f_star[i] = f_star[i] + f_k[i];
		}
	}

	//------------- Recoloracao de Latva-Kokko ( eqs. 9 e 10 ) -----------------------------------//

	const double conc_R = rho_R / rho;
	const double conc_B = 1.0 - conc_R;

	const double fator = parameters.recoll * conc_R * conc_B;

	for ( int i = 0; i < nvel; i++ )
	{
		//  soma_k f^eq_i( rho_k , alpha_k , u = 0 ) = soma_k rho_k phi^k_i

		double phi_R, phi_B;

		if      ( i == 0 ) { phi_R = alfa_R;                  phi_B = alfa_B;                  }
		else if ( i <  7 ) { phi_R = ( 1.0 - alfa_R ) / 12.0; phi_B = ( 1.0 - alfa_B ) / 12.0; }
		else               { phi_R = ( 1.0 - alfa_R ) / 24.0; phi_B = ( 1.0 - alfa_B ) / 24.0; }

		const double soma_eq = rho_R * phi_R + rho_B * phi_B;

		double cos_phi = 0.0;

		if ( i > 0 && mod_grad > eps )
		{
			const double cx = lattice.c_i[ i * dim + 0 ];
			const double cy = lattice.c_i[ i * dim + 1 ];
			const double cz = ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0;

			const double mod_c = sqrt ( cx * cx + cy * cy + cz * cz );

			cos_phi = ( cx * n_x + cy * n_y + cz * n_z ) / mod_c;
		}

		const double desvio = fator * cos_phi * soma_eq;

		f_R[i] = conc_R * f_star[i] + desvio;
		f_B[i] = conc_B * f_star[i] - desvio;
	}
}



//====================================================================================================================//

//------ Verificacao do modelo na partida ---------------------------------------------------------------------------//
//
//  Confere, com numeros, tres coisas que sao faceis de errar ao transcrever o artigo:
//
//    1) M^-1 M = I;
//    2) mom_eq_LYL() ( forma fechada ) == M . dist_eq_LYL() ;
//    3) as eqs. (20) e (23): segundo momento = rho u u + p_k I, e terceiro momento fora da diagonal
//       carregando p_k -- que e justamente o que o termo de ordem alta existe para garantir.

bool verifica_LYL ( LATTICE lattice, double alfa_k )
{
	bool ok = true;

	//----- 1) M^-1 M = I ------------------------------------------------------------------------//

	double err_inv = 0.0;

	for ( int a = 0; a < nvel; a++ )
	{
		for ( int b = 0; b < nvel; b++ )
		{
			double s = 0.0;

			for ( int i = 0; i < nvel; i++ ) s = s + MRT_LYL_M_INV[a][i] * MRT_LYL_M[i][b];

			err_inv = max ( err_inv, fabs ( s - ( a == b ? 1.0 : 0.0 ) ) );
		}
	}

	cout << "   max | M^-1 M - I |               = " << err_inv << endl;

	if ( err_inv > 1.0e-12 ) ok = false;

	//----- 2) e 3) sobre um estado arbitrario ---------------------------------------------------//

	const double rho_k = 1.37, ux = 0.031, uy = -0.017, uz = 0.023;

	double f_eq[nvel], m_eq[nvel];

	dist_eq_LYL ( f_eq, ux, uy, uz, rho_k, alfa_k, 1.0, lattice );

	mom_eq_LYL  ( m_eq, ux, uy, uz, rho_k, alfa_k );

	double err_m = 0.0;

	for ( int a = 0; a < nvel; a++ )
	{
		double s = 0.0;

		for ( int i = 0; i < nvel; i++ ) s = s + MRT_LYL_M[a][i] * f_eq[i];

		err_m = max ( err_m, fabs ( s - m_eq[a] ) );
	}

	cout << "   max | M f^eq - m^eq |            = " << err_m << endl;

	if ( err_m > 1.0e-12 ) ok = false;

	//----- momentos ---------------------------------------------------------------------------//

	const double cs2_k = 0.5 * ( 1.0 - alfa_k );

	const double p_k = rho_k * cs2_k;

	const double u[3] = { ux, uy, uz };

	double err_2 = 0.0, err_3 = 0.0;

	for ( int a = 0; a < dim; a++ )
	{
		for ( int b = 0; b < dim; b++ )
		{
			double s = 0.0;

			for ( int i = 0; i < nvel; i++ )
				s = s + lattice.c_i[ i * dim + a ] * lattice.c_i[ i * dim + b ] * f_eq[i];

			err_2 = max ( err_2, fabs ( s - ( rho_k * u[a] * u[b] + ( a == b ? p_k : 0.0 ) ) ) );

			for ( int g = 0; g < dim; g++ )
			{
				double t = 0.0;

				for ( int i = 0; i < nvel; i++ )
					t = t + lattice.c_i[ i * dim + a ] * lattice.c_i[ i * dim + b ]
					      * lattice.c_i[ i * dim + g ] * f_eq[i];

				//  eq. (23): fora da diagonal a pressao e p_k;  na diagonal continua rho_k c^2/3

				const double fat = ( a == b && b == g ) ? rho_k / 3.0 : p_k;

				const double alvo = fat * ( u[a] * ( b == g ) + u[b] * ( a == g ) + u[g] * ( a == b ) );

				err_3 = max ( err_3, fabs ( t - alvo ) );
			}
		}
	}

	cout << "   max erro no 2o momento, eq.(20)  = " << err_2 << endl;
	cout << "   max erro no 3o momento, eq.(23)  = " << err_3 << endl;

	if ( err_2 > 1.0e-12 || err_3 > 1.0e-12 ) ok = false;

	return ok;
}

#endif   // nvel == 19

//====================================================================================================================//


//=============================== Collision step using the Spencer, Halliday & Care model (immicisble) ===============//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

//=============================== Recoloracao de Halliday, Hollis & Care ( 2007 ) ====================================//
//
//      A regra que Spencer, Halliday & Care usam de fato no artigo de 2010.  Na Sec. II eles nao a
//      escrevem -- dizem so "we therefore employ identical segregation rules to those analyzed
//      previously [14]", e [14] e  Halliday, Hollis & Care, Phys. Rev. E 76, 026708 (2007).  O unico
//      lugar em que ela aparece explicita e a eq. (69), a generalizacao para N componentes:
//
//          f_i^{n+} = ( rho_n / rho ) f_i  +  soma_{m != n}  beta_nm  t_i ( rho_n rho_m / rho^2 ) c_i . n_nm
//
//      que com dois fluidos da
//
//          f_i^{R+} = conc_R f_i  +  beta conc_R conc_B rho w_i ( c_i . n )
//
//      A DIFERENCA para recolloring(), a de Latva-Kokko, e UMA so: aqui e  c_i . n , sem dividir por
//      | c_i | .  Latva-Kokko e D'Ortona usam  cos( phi_i ) = ( c_i . n ) / | c_i | ; SHC aparecem
//      apenas como linhagem no artigo ( "derived from that first introduced by D'Ortona et al. [24]
//      and, later, Latva-Kokko and Rothman [25]" ), nao como a regra empregada.  E a mesma regra das
//      eqs. (18)-(19) de Saito et al. (2023), que a atribuem a Halliday et al.
//
//      A eq. (69) impressa nao traz o fator rho ( a amplitude e  t_i rho_n rho_m / rho^2 ).  Aqui ele
//      esta, para casar com a forma de Saito e para que a UNICA diferenca em relacao a recolloring()
//      seja o  1 / | c_i | ;  com rho ~ 1, que e o caso destes ensaios, a distincao e de ordem do
//      salto de Laplace ( < 1 % ).
//
//      QUANTO MUDA.  O fluxo de segregacao  soma_i d_i c_i  vale  beta conc_R conc_B rho cs^2  nesta
//      forma e  beta conc_R conc_B rho S  na de Latva-Kokko, com
//
//          S = soma_i w_i ( c_i . n )^2 / | c_i |
//
//      As duas sao isotropicas ( conferido em 64 direcoes, 2e-16 ), mas S < cs^2 :
//
//          D2Q9   S / cs^2 = 0.90237        D3Q19  S / cs^2 = 0.80474        D3Q27  0.82286
//
//      Ou seja: no mesmo beta a recoloracao de Latva-Kokko segrega menos e deixa a interface mais
//      grossa -- 24 % mais grossa no D3Q19, que e a rede dos programas de SHC desta arvore.
//
//      E ha um efeito colateral que nao e obvio: esta forma degenera EXATAMENTE de  nz = 1  para o
//      D2Q9 ( o fluxo vale cs^2 nas tres redes ), e a de Latva-Kokko NAO ( 0.2682 e 0.2743 contra
//      0.3008 ).  Um caso "2D" rodado como nz = 1 com Latva-Kokko nao e o modelo D2Q9 verdadeiro.
//
//      Assinatura identica a de recolloring(), e a mesma convencao de sinal: 'grad' e o que os
//      mediadores devolvem, ou seja  - grad( rho^N ) .
//
//      Input: f, grad, |grad|, fracoes de massa, rho, beta
//      Output: f_R, f_B

#pragma acc routine seq
void recolloring_HHC ( double* f, double* f_R, double* f_B, double* grad, double mod_mM,
						double conc_R, double conc_B, double rho, double beta, LATTICE lattice )
{
	const double fator = beta * conc_R * conc_B;

	f_R[0] = conc_R * f[0];
	f_B[0] = conc_B * f[0];

	for ( int i = 1; i < nvel; i++ )
	{
		const double *c_i = lattice.c_i + i * dim;

		const double prod = c_i[0] * grad[0] + c_i[1] * grad[1]
		                    + ( ( dim == 3 ) ? c_i[2] * grad[2] : 0.0 );

		//  c_i . n , com  n = grad / | grad | .  Sem o  1 / | c_i |  de Latva-Kokko.

		const double ci_n = ( mod_mM > 0.0 ) ? prod / mod_mM : 0.0;

		const double d = fator * lattice.w[i] * rho * ci_n;

		f_R[i] = conc_R * f[i] - d;
		f_B[i] = conc_B * f[i] + d;
	}
}

//====================================================================================================================//



//=============================== Seletor da recoloracao dos modelos de SHC =========================================//
//
//      parameters.recoll_halliday = 0  ( padrao )  Latva-Kokko, recolloring()  -- eq. (1.3) do
//                                                  relatorio LMPT, e o que esta arvore sempre usou.
//                                 = 1              Halliday, Hollis & Care    -- eq. (69) do artigo
//                                                  de 2010, recolloring_HHC().

#pragma acc routine seq
inline void recolloring_SHC ( double* f, double* f_R, double* f_B, double* grad, double mod_mM,
								double conc_R, double conc_B, double rho, double beta,
								LATTICE lattice, PARAMETERS parameters )
{
	if ( parameters.recoll_halliday )
		recolloring_HHC ( f, f_R, f_B, grad, mod_mM, conc_R, conc_B, rho, beta, lattice );
	else
		recolloring     ( f, f_R, f_B, grad, mod_mM, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//



void coll_SHC_2010 ( double* f_R, double* f_B, double* f_m, double gx, double gy, double gz, LATTICE lattice, 
					PARAMETERS parameters )
{
	//------------- Calcula densidades e concentrações -------------------------------------------//
	
	double rho_R = density ( f_R );
	double rho_B = density ( f_B );
	
	double rho = rho_R + rho_B;
	
	double conc_R = rho_R / rho; 
	double conc_B = 1. - conc_R;
			
	//------------- Calcula a direção do gradiente -----------------------------------------------//

    double grad_x, grad_y, grad_z;
            
    momentum ( f_m, grad_x, grad_y, grad_z, lattice );
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
        
    double n[3] = {grad_uni_x, grad_uni_y, grad_uni_z};
   
    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );
   
    //------------------- Colisão monofásica -----------------------------------------------------//
    
    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    
    else parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
    coll_BGK ( f, gx, gy, gz, lattice, parameters );

	//------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//
	
	double fat_int = parameters.A_fact;
	
	double tau = parameters.tau;

    double fator = rho * ( 1. / tau ) * fat_int * mod_grad * ( 3.0 );

    interf_tension_SHC ( f, fator, n, lattice );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {grad_x, grad_y, grad_z};    

	double beta = parameters.recoll;
	
    recolloring_SHC ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice, parameters );
}

//====================================================================================================================//





//=============================== Collision step using the Reis-Phillips/Leclaire model (immiscible) ================//
//
//      Input: red and blue distribution functions, mediator distribution, acceleration, lattice and parameters
//      Output: post-collisional red and blue distribution functions
//
//      This version mirrors the mediator-based coll_SHC_2010 implementation. The interface normal/gradient is obtained
//      from the mediator population f_m through momentum(f_m,...), and the RPL perturbation is applied after the BGK
//      collision of the mixture and before Latva-Kokko recoloring.
//
//      parameters.A_fact is interpreted as the RPL perturbation amplitude A.
//      With the usual D3Q19 normalization, sigma_lu ~= (2/9) A.
//
//====================================================================================================================//

#pragma acc routine seq
void coll_RPL_D3Q19 ( double* f_R, double* f_B, double* f_m, double gx, double gy, double gz, LATTICE lattice,
                      PARAMETERS parameters )
{
    constexpr double eps = 1.0e-30;

    //------------- Calcula densidades e concentrações -------------------------------------------//

    double rho_R = density ( f_R );
    double rho_B = density ( f_B );

    double rho = rho_R + rho_B;

    if ( rho <= eps ) return;

    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;

    //------------- Calcula o gradiente a partir dos mediadores ----------------------------------//

    double grad_x, grad_y, grad_z;

    momentum ( f_m, grad_x, grad_y, grad_z, lattice );

    double grad[3] = { grad_x, grad_y, grad_z };

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    //------------------- Colisão monofásica da mistura ------------------------------------------//

    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    else                   parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

    coll_BGK ( f, gx, gy, gz, lattice, parameters );

    //------------------- Imposição da tensão interfacial RPL ------------------------------------//

    double A = parameters.A_fact;

    if ( mod_grad > eps && A != 0.0 )
    {
        interf_tension_RPL ( f, A, mod_grad, grad, lattice );
    }

    //------------------- Etapa de recoloração (Latva-Kokko) -------------------------------------//

    double beta = parameters.recoll;

    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//

//=============================== Collision step using the Spencer, Halliday & Care model (immicisble) ===============//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

//=============================== Modelo de Latva-Kokko & Rothman ( 2005 ) ===========================================//
//
//      M. Latva-Kokko & D. H. Rothman,
//        "Diffusion properties of gradient-based lattice Boltzmann models of immiscible fluids",
//        Phys. Rev. E 71, 056702 (2005);
//        "Static contact angle in lattice Boltzmann models of immiscible fluids",
//        Phys. Rev. E 72, 046701 (2005).
//
//      Este e o modelo de REFERENCIA desta arvore -- todos os outros color-gradient daqui descendem
//      dele --, e vale ser explicito sobre o que ele e e sobre o que ja estava implementado.
//
//      O passo bifasico tem tres partes ( secao II do primeiro artigo ):
//
//      1) BGK monofasico da populacao total,  eq. (5):  N'_i = ( 1 + lambda ) N_i - lambda N_i^eq .
//         Com  lambda = -1/tau  e o BGK de sempre, e o equilibrio e o classico ( c_s^2 = 1/3,
//         p = rho / 3 ).  Nada de velocidade do som variavel: este modelo e de densidades iguais;
//
//      2) PERTURBACAO de Gunstensen, eq. (8) do artigo do angulo de contato:
//
//             N''_i = N'_i + A | f | [ ( c_i . f )^2 / ( f . f ) - 1/2 ]
//
//         com o gradiente de cor da eq. (9),  f = soma_i c_i [ rho_R - rho_B ]( x + c_i ) ;
//
//      3) RECOLORACAO, eqs. (9)-(10) do artigo da difusao -- a contribuicao central do trabalho:
//
//             R_i = ( R / rho ) N''_i + beta ( R B / rho^2 ) N_i^eq( rho, 0 ) cos( phi_i )
//             B_i = ( B / rho ) N''_i - beta ( R B / rho^2 ) N_i^eq( rho, 0 ) cos( phi_i )
//
//         que substitui a recoloracao "de maximo trabalho" de Gunstensen -- aquela que enche as
//         direcoes mais proximas ao gradiente em ordem, e que PRENDE a interface na rede.  Com esta,
//         as cores se misturam moderadamente e a distribuicao fica simetrica em torno do gradiente.
//
//      DUAS COISAS JA ESTAVAM AQUI, com outro nome:
//
//        * a etapa (3) e exatamente  recolloring()  de LBM_functions.cpp.  E o que os modelos de
//          Spencer-Halliday-Care, Reis-Phillips e Liu-Valocchi-Kang desta arvore ja usavam;
//        * a etapa (2), escrita para uma rede COM PESOS, e  interf_tension_LVK()  com chi = 1/2 --
//          isto e, o mesmo  interf_tension_RPL().  A forma  [ ( c_i . f )^2 / | f |^2 - 1/2 ]  do
//          artigo e a versao FCHC, onde todas as 24 direcoes tem o mesmo peso; num D3Q19 com pesos
//          o operador conservativo e  A | f | [ w_i ( c_i . n )^2 - B_i ] , e os B_i formam a
//          familia de um parametro da eq. (29) de Liu-Valocchi-Kang.  Que a tensao interfacial nao
//          depende de chi ja foi conferido numericamente ( ver LIU_VALOCCHI_KANG.md ).
//
//      TENSAO INTERFACIAL.  A eq. (22) do segundo artigo,  sigma = - 192 A rho / lambda , vale para
//      a rede FCHC ( b = 24, c^2 = 2, D = 4 ) e NAO se transporta para o D3Q19.  Aqui, com o
//      gradiente isotropico dos mediadores, o mesmo calculo tensorial da
//
//          soma_i Omega_i c_ia c_ib = ( 2 A / 9 ) | f | ( n_a n_b - delta_ab )
//
//          =>   sigma = ( 4 / 9 ) A tau | f | / | grad rho^N |  =  ( 4 / 9 ) A tau rho
//
//      porque aqui  f = grad( rho_R - rho_B ) = rho grad( rho^N )  para rho uniforme.  Com rho = 1
//      isso e a mesma constante do modelo de Liu-Valocchi-Kang, como tem de ser: o operador e o
//      mesmo, so a normalizacao do gradiente muda.
//
//      MOLHABILIDADE.  Pela eq. (12) do segundo artigo o solido carrega  N_red,wall = rho p_red , e a
//      eq. (28) da  cos( theta ) = p_red .  Nesta arvore o valor de parede e  parameters.wett_R ,
//      injetado por propag_site_med() no lugar do vizinho solido, e a fase levada pelos mediadores e
//      rho^N -- entao  wett_R = p_red  e  cos( theta ) = wett_R  diretamente.  wett_R = 0 e
//      molhamento neutro; +-1, molhamento total por uma das cores.
//
//      A CURVA INTERFACIAL.  A eq. (15) do artigo da difusao,  d phi / ds = K phi ( 1 - phi ) , e a
//      previsao mais util deste modelo para o resto da arvore: ela FIXA a espessura da interface a
//      partir de beta, e a espessura e o que domina a corrente espuria em todos os modelos medidos
//      ate aqui.  calc_perfil_logistico(), em Other_functions.cpp, mede K e o residuo do ajuste.
//
//      Input: f_R, f_B, mediadores, aceleracao de corpo, lattice, parameters
//      Output: f_R, f_B pos-colisao
//
//====================================================================================================================//

#pragma acc routine seq
void coll_LKR_D3Q19 ( double* f_R, double* f_B, double* f_m, double gx, double gy, double gz,
						LATTICE lattice, PARAMETERS parameters )
{
	constexpr double eps = 1.0e-30;

	const double rho_R = density ( f_R );
	const double rho_B = density ( f_B );

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	const double conc_R = rho_R / rho;
	const double conc_B = 1.0 - conc_R;

	//------------- Gradiente de cor, eq. (9), pelos mediadores ----------------------------------//
	//
	//  Os mediadores entregam  - grad( fase ) .  A fase deste modelo e rho^N, entao o gradiente de
	//  cor do artigo,  f = grad( rho_R - rho_B ) , e  rho  vezes este, para rho uniforme.  Como o
	//  operador de perturbacao usa so a DIRECAO de f e o modulo entra linearmente em A | f |, a
	//  diferenca e absorvida na calibracao de A -- que o cabecalho declara: sigma = 4 A tau rho / 9.

	double grad_x, grad_y, grad_z;

	momentum ( f_m, grad_x, grad_y, grad_z, lattice );

	double grad[3] = { grad_x, grad_y, grad_z };

	const double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

	//------------- (1) BGK monofasico, eq. (5) --------------------------------------------------//
	//
	//  tau interpolado linearmente na cor, que e o que a eq. (21) de modelos posteriores faz.  Com
	//  tau_R = tau_B -- o caso do artigo, em que os dois fluidos so diferem pela cor -- nao muda nada.

	const double rho_N = ( rho_R - rho_B ) / rho;

	parameters.tau = 0.5 * ( 1.0 + rho_N ) * parameters.tau_R
	               + 0.5 * ( 1.0 - rho_N ) * parameters.tau_B;

	double f[nvel];

	for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

	coll_BGK ( f, gx, gy, gz, lattice, parameters );

	//------------- (2) Perturbacao de Gunstensen, eq. (8) ---------------------------------------//
	//
	//  chi = 1/2 e o membro da familia que corresponde ao  - 1/2  do artigo na rede FCHC.

	if ( mod_grad > eps && parameters.A_fact != 0.0 )
	{
		interf_tension_LVK ( f, parameters.A_fact * rho, mod_grad, grad, 0.5, lattice );
	}

	//------------- (3) Recoloracao de Latva-Kokko, eqs. (9)-(10) --------------------------------//

	recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, parameters.recoll, lattice );
}

//====================================================================================================================//


//=============================== Modelo de Saito et al. ( 2023 ), D3Q27 =============================================//
//
//      Saito, Takada, Baba, Someya & Ito, Phys. Rev. E 108, 065305 (2023),
//      "Generalized equilibria for color-gradient lattice Boltzmann model based on higher-order
//       Hermite polynomials: A simplified implementation with central moments".
//
//      Este bloco implementa o modelo inteiro.  Ele difere dos outros color-gradient desta arvore
//      em tres pontos, e vale separar o que e novidade do artigo do que e so heranca:
//
//      1) A TENSAO INTERFACIAL NAO VEM DE OPERADOR DE PERTURBACAO.  Vem de uma forca de corpo do
//         tipo continuum surface force, eq. (15):
//
//             F_s = ( 1/2 ) sigma kappa grad( phi ) ,      kappa = - div( n ) ,   n = grad(phi)/|grad(phi)|
//
//         Isto e: sigma e ENTRADA, nao resultado de calibracao.  Nao ha "sigma = C A" para medir --
//         o teste da gota estatica serve para CONFERIR o modelo, nao para calibra-lo.  A curvatura
//         nao e local ( precisa dos vizinhos de n ), entao quem a calcula e o programa.
//
//      2) O EQUILIBRIO GENERALIZADO, eq. (56), que e a contribuicao do artigo.  Em espaco de
//         populacoes ele e um polinomio de Hermite ate SEXTA ordem em u, com dezenas de termos.  Em
//         espaco de MOMENTOS CENTRAIS ele e isto, e nada mais ( eqs. 57-63 ):
//
//             k000 = rho ;   k200 = k020 = k002 = p ;
//             k220 = k202 = k022 = p cs^2 ;   k222 = p cs^4 ;   todos os outros = 0
//
//         INDEPENDENTE DA VELOCIDADE.  E dai que sai a invariancia galileana: nos equilibrios
//         anteriores ( modelos A, B e C do artigo ) sobram termos proporcionais a ( p - rho cs^2 )
//         vezes potencias de u nos momentos de terceira ordem e acima, e sao eles que deformam uma
//         gota que se move.
//
//         Por isso esta implementacao NAO escreve a eq. (56): escreve os sete numeros acima e
//         volta ao espaco de populacoes pela transformada inversa.  E o mesmo equilibrio, exato, e
//         cabe em cinco linhas em vez de trinta termos.
//
//      3) A RECOLORACAO e a de Halliday, eqs. (18)-(19):
//
//             f_i^{r,*} = ( rho_r / rho ) f_i^*  +  ( w_i / cs^2 ) c_i . R ,    R = beta ( rho_r rho_b / rho^2 ) p n
//
//         Repare que e  c_i . R , e nao  |c_i| cos(phi_i) : em relacao a recoloracao de
//         Latva-Kokko usada nos outros modelos desta arvore, sobra um fator | c_i | nas direcoes
//         diagonais.  E a amplitude e p, a pressao, e nao rho w_i.
//
//      A pressao e  p = soma_k rho_k ( cs^k )^2 , eq. (12), com uma velocidade do som por fluido --
//      e dela que vem a razao de densidades, eq. (13):  gamma = ( cs^b )^2 / ( cs^r )^2 .
//
//      A rede e D3Q27 ( eq. 1 ).  Numa caixa com nz = 1 ela degenera EXATAMENTE no D2Q9:
//      axial 2/27 + 2 x 1/54 = 1/9 , diagonal 1/54 + 2 x 1/216 = 1/36 , repouso 8/27 + 2 x 2/27 = 4/9.
//      O teste da gota estatica do artigo ( secao IV A, Tabela II ) e justamente bidimensional.
//
//====================================================================================================================//

//------------------ Indice da rede a partir do trio ( a, b, c ) ----------------------------------------------------//
//
//  a, b, c em { 0, 1, 2 }  correspondem a  c_alpha em { 0, +1, -1 }.  As 27 velocidades do D3Q27 sao
//  o produto tensorial dos tres eixos, e e isso que faz a transformada de momentos centrais ser
//  SEPARAVEL: em vez de uma matriz 27 x 27, tres passagens de tamanho 3, uma por eixo.
//
//  SAITO_MAP[ a * 9 + b * 3 + c ] = indice i da velocidade ( c_x, c_y, c_z ) nesta biblioteca.

#define SAITO_NM 27

static const int SAITO_MAP[27] =
{
	 0,  5,  6,
	 3, 16, 18,
	 4, 17, 15,
	 1, 11, 13,
	 7, 19, 24,
	 9, 25, 22,
	 2, 14, 12,
	10, 21, 26,
	 8, 23, 20
};

#pragma acc declare copyin( SAITO_MAP )

//------------------ Transformada direta: populacoes -> momentos centrais -------------------------------------------//
//
//      k_abc = soma_i f_i ( c_ix - u_x )^a ( c_iy - u_y )^b ( c_iz - u_z )^c ,   a, b, c em { 0, 1, 2 }
//
//  Feita eixo a eixo.  Num eixo, com f0, f+ e f- as tres populacoes daquela fibra:
//
//      m0 = f0 + f+ + f-        m1 = f+ - f-        m2 = f+ + f-
//      k0 = m0                  k1 = m1 - u m0      k2 = m2 - 2 u m1 + u^2 m0
//
//  Saida em k[ a * 9 + b * 3 + c ].

#pragma acc routine seq
void cm_saito ( const double *f, double ux, double uy, double uz, double *k )
{
	double t[27], s[27];

	//------ eixo x: a fibra e ( a, b, c ) com b, c fixos ------------------------------------------//

	for ( int b = 0; b < 3; b++ )
	{
		for ( int c = 0; c < 3; c++ )
		{
			const double f0 = f[ SAITO_MAP[ 0 * 9 + b * 3 + c ] ];
			const double fp = f[ SAITO_MAP[ 1 * 9 + b * 3 + c ] ];
			const double fm = f[ SAITO_MAP[ 2 * 9 + b * 3 + c ] ];

			const double m0 = f0 + fp + fm;
			const double m1 = fp - fm;
			const double m2 = fp + fm;

			t[ 0 * 9 + b * 3 + c ] = m0;
			t[ 1 * 9 + b * 3 + c ] = m1 - ux * m0;
			t[ 2 * 9 + b * 3 + c ] = m2 - 2.0 * ux * m1 + ux * ux * m0;
		}
	}

	//------ eixo y ---------------------------------------------------------------------------------//

	for ( int a = 0; a < 3; a++ )
	{
		for ( int c = 0; c < 3; c++ )
		{
			const double f0 = t[ a * 9 + 0 * 3 + c ];
			const double fp = t[ a * 9 + 1 * 3 + c ];
			const double fm = t[ a * 9 + 2 * 3 + c ];

			const double m0 = f0 + fp + fm;
			const double m1 = fp - fm;
			const double m2 = fp + fm;

			s[ a * 9 + 0 * 3 + c ] = m0;
			s[ a * 9 + 1 * 3 + c ] = m1 - uy * m0;
			s[ a * 9 + 2 * 3 + c ] = m2 - 2.0 * uy * m1 + uy * uy * m0;
		}
	}

	//------ eixo z ---------------------------------------------------------------------------------//

	for ( int a = 0; a < 3; a++ )
	{
		for ( int b = 0; b < 3; b++ )
		{
			const double f0 = s[ a * 9 + b * 3 + 0 ];
			const double fp = s[ a * 9 + b * 3 + 1 ];
			const double fm = s[ a * 9 + b * 3 + 2 ];

			const double m0 = f0 + fp + fm;
			const double m1 = fp - fm;
			const double m2 = fp + fm;

			k[ a * 9 + b * 3 + 0 ] = m0;
			k[ a * 9 + b * 3 + 1 ] = m1 - uz * m0;
			k[ a * 9 + b * 3 + 2 ] = m2 - 2.0 * uz * m1 + uz * uz * m0;
		}
	}
}

//------------------ Transformada inversa: momentos centrais -> populacoes ------------------------------------------//
//
//  Num eixo:   m0 = k0 ,  m1 = k1 + u k0 ,  m2 = k2 + 2 u k1 + u^2 k0
//              f0 = m0 - m2 ,  f+ = ( m2 + m1 ) / 2 ,  f- = ( m2 - m1 ) / 2

#pragma acc routine seq
void cm_saito_inv ( const double *k, double ux, double uy, double uz, double *f )
{
	double t[27], s[27];

	//------ eixo z ---------------------------------------------------------------------------------//

	for ( int a = 0; a < 3; a++ )
	{
		for ( int b = 0; b < 3; b++ )
		{
			const double k0 = k[ a * 9 + b * 3 + 0 ];
			const double k1 = k[ a * 9 + b * 3 + 1 ];
			const double k2 = k[ a * 9 + b * 3 + 2 ];

			const double m0 = k0;
			const double m1 = k1 + uz * k0;
			const double m2 = k2 + 2.0 * uz * k1 + uz * uz * k0;

			t[ a * 9 + b * 3 + 0 ] = m0 - m2;
			t[ a * 9 + b * 3 + 1 ] = 0.5 * ( m2 + m1 );
			t[ a * 9 + b * 3 + 2 ] = 0.5 * ( m2 - m1 );
		}
	}

	//------ eixo y ---------------------------------------------------------------------------------//

	for ( int a = 0; a < 3; a++ )
	{
		for ( int c = 0; c < 3; c++ )
		{
			const double k0 = t[ a * 9 + 0 * 3 + c ];
			const double k1 = t[ a * 9 + 1 * 3 + c ];
			const double k2 = t[ a * 9 + 2 * 3 + c ];

			const double m0 = k0;
			const double m1 = k1 + uy * k0;
			const double m2 = k2 + 2.0 * uy * k1 + uy * uy * k0;

			s[ a * 9 + 0 * 3 + c ] = m0 - m2;
			s[ a * 9 + 1 * 3 + c ] = 0.5 * ( m2 + m1 );
			s[ a * 9 + 2 * 3 + c ] = 0.5 * ( m2 - m1 );
		}
	}

	//------ eixo x ---------------------------------------------------------------------------------//

	for ( int b = 0; b < 3; b++ )
	{
		for ( int c = 0; c < 3; c++ )
		{
			const double k0 = s[ 0 * 9 + b * 3 + c ];
			const double k1 = s[ 1 * 9 + b * 3 + c ];
			const double k2 = s[ 2 * 9 + b * 3 + c ];

			const double m0 = k0;
			const double m1 = k1 + ux * k0;
			const double m2 = k2 + 2.0 * ux * k1 + ux * ux * k0;

			f[ SAITO_MAP[ 0 * 9 + b * 3 + c ] ] = m0 - m2;
			f[ SAITO_MAP[ 1 * 9 + b * 3 + c ] ] = 0.5 * ( m2 + m1 );
			f[ SAITO_MAP[ 2 * 9 + b * 3 + c ] ] = 0.5 * ( m2 - m1 );
		}
	}
}

//------------------ Equilibrio generalizado, eq. (56) --------------------------------------------------------------//
//
//  Escrito pelos momentos centrais das eqs. (57)-(63), que sao independentes da velocidade, e
//  trazido de volta ao espaco de populacoes.  E EXATAMENTE a eq. (56) -- a transformada e bijetiva.
//
//      k000 = rho ;  k200 = k020 = k002 = p ;  k220 = k202 = k022 = p cs^2 ;  k222 = p cs^4
//
//  Conferido em runtime por check_saito(): soma f^eq = rho, soma f^eq c = rho u, e os momentos de
//  terceira ordem da eq. (64).
//
//      Input: velocidade, densidade total, pressao, lattice
//      Output: f_eq[27]

#pragma acc routine seq
void dist_eq_saito ( double *f_eq, double ux, double uy, double uz, double rho, double p, LATTICE lattice )
{
	const double cs2 = lattice.c_s2;

	double k[27];

	for ( int n = 0; n < 27; n++ ) k[n] = 0.0;

	k[ 0 * 9 + 0 * 3 + 0 ] = rho;						// k000

	k[ 2 * 9 + 0 * 3 + 0 ] = p;							// k200
	k[ 0 * 9 + 2 * 3 + 0 ] = p;							// k020
	k[ 0 * 9 + 0 * 3 + 2 ] = p;							// k002

	k[ 2 * 9 + 2 * 3 + 0 ] = p * cs2;					// k220
	k[ 2 * 9 + 0 * 3 + 2 ] = p * cs2;					// k202
	k[ 0 * 9 + 2 * 3 + 2 ] = p * cs2;					// k022

	k[ 2 * 9 + 2 * 3 + 2 ] = p * cs2 * cs2;				// k222

	cm_saito_inv ( k, ux, uy, uz, f_eq );
}

//====================================================================================================================//



//=============================== Campo auxiliar do modelo de Saito =================================================//
//
//      A forca de tensao interfacial e a correcao Q envolvem DIVERGENCIAS, e divergencia le os
//      vizinhos.  Esta funcao prepara, sitio a sitio, o que o programa precisa juntar depois:
//
//          C[0..2] = n = grad(phi) / |grad(phi)|          ( eq. 17, aponta para o vermelho )
//          C[3..5] = ( p - rho cs^2 ) u                    ( o que entra em Q, eq. 65 )
//
//      O gradiente vem dos mediadores, que entregam  - grad(phi) : dai o sinal trocado em n.
//
//      CORTE: longe da interface  grad(phi)  e so ruido de arredondamento e n viraria um campo de
//      modulo 1 apontando para todo lado, cuja divergencia e enorme.  Abaixo de 'corte' zeramos n.
//      Como a forca e proporcional a grad(phi) LOCAL -- e nao a n --, o corte nao deixa degrau: ele
//      so evita que a curvatura de um vizinho de fora da interface contamine o gather.
//
//      Input: f_R, f_B, mediadores, u, cs2_R, cs2_B, corte, lattice
//      Output: C[6]

#pragma acc routine seq
void campo_saito ( double *f_R, double *f_B, double *f_m, double ux, double uy, double uz,
					double *C, LATTICE lattice, PARAMETERS parameters, double corte )
{
	constexpr double eps = 1.0e-30;

	for ( int n = 0; n < 6; n++ ) C[n] = 0.0;

	const double rho_R = density ( f_R );
	const double rho_B = density ( f_B );

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	double gx, gy, gz;

	momentum ( f_m, gx, gy, gz, lattice );				// = - grad( phi )

	const double mod_grad = sqrt ( gx * gx + gy * gy + gz * gz );

	if ( mod_grad >= corte )
	{
		C[0] = - gx / mod_grad;							// n = grad(phi) / |grad(phi)|
		C[1] = - gy / mod_grad;
		C[2] = - gz / mod_grad;
	}

	const double p = rho_R * parameters.cs2_R + rho_B * parameters.cs2_B;

	const double dp = p - rho * lattice.c_s2;

	C[3] = dp * ux;
	C[4] = dp * uy;
	C[5] = dp * uz;
}

//====================================================================================================================//



//=============================== Recoloracao de Halliday, eqs. (18)-(19) ===========================================//
//
//          f_i^{r,*} = ( rho_r / rho ) f_i^*  +  ( w_i / cs^2 ) c_i . R
//          f_i^{b,*} = ( rho_b / rho ) f_i^*  -  ( w_i / cs^2 ) c_i . R
//
//          R = beta ( rho_r rho_b / rho^2 ) p n
//
//      Conserva f_i^r + f_i^b = f_i^* exatamente, e a soma_i das duas correcoes e zero por paridade:
//      massa e momento TOTAIS intactos.  Cada cor ganha momento -- e o fluxo de segregacao.
//
//      Duas diferencas para recolloring(), a de Latva-Kokko usada nos outros modelos:
//        * aqui e  c_i . n , sem dividir por | c_i | -- as diagonais pesam | c_i | vezes mais;
//        * a amplitude e a PRESSAO p, e nao rho w_i .
//
//      n vem na convencao da eq. (17), apontando para o vermelho.
//
//      -------------------------------------------------------------------------------------------
//      POSITIVIDADE  ( limita = parameters.recoll_positivo )
//
//      f_i^{b,*} >= 0  exige      f_i / ( rho w_i )  >=  beta ( rho_r / rho ) ( p / rho cs^2 ) c_i.n
//
//      Em repouso o lado esquerdo vale 1 ( o equilibrio deste modelo e  f_i = rho w_i  quando
//      p = rho cs^2 ), o pior elo tem  c_i . n = sqrt(2)  no plano e  sqrt(3)  em 3D, e a condicao
//      da  beta <= 1/sqrt(2) = 0.707  em 2D -- beta = 0.7 passa raspando.  Numa gota que se move
//      com velocidade u o lado esquerdo cai, nos elos contra o movimento, para
//
//          g( u ) = 1 - 3 u + 3 u^2      ( 0.86 em u = 0.05 ,  0.62 em u = 0.15 ,  0.44 em u = 0.25 )
//
//      e a folga acaba: a partir de u ~ 0.1 aparecem populacoes de cor negativas na cauda da
//      interface, phi sai de [ -1, 1 ], | grad phi | vira ruido no seio do fluido, n passa a ser um
//      versor aleatorio em TODO o dominio e a curvatura -- e com ela a forca CSF -- fica espuria.
//      A gota freia e se deforma.  Nao e o equilibrio generalizado que falha ( com cs2_R = cs2_B os
//      momentos centrais de equilibrio sao exatamente independentes de u e Q e identicamente nulo ):
//      e a recoloracao.
//
//      O limitador multiplica todos os d_i do sitio pelo MESMO lambda <= 1, o maior que mantem as
//      duas cores nao negativas.  Sendo um fator global do sitio, soma_i d_i continua zero: a massa
//      de cada cor segue exatamente conservada, e o momento total tambem.  Onde nao ha ameaca a
//      positividade lambda = 1 e nada muda -- em repouso o resultado e bit a bit o de limita = 0.

#pragma acc routine seq
void recolloring_saito ( double *f, double *f_R, double *f_B, double *n, double rho_R, double rho_B,
							double p, double beta, int limita, int forma, LATTICE lattice )
{
	constexpr double eps = 1.0e-30;

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	const double conc_R = rho_R / rho;
	const double conc_B = 1.0 - conc_R;

	//------------- Amplitude, conforme a forma escolhida -----------------------------------------//
	//
	//  forma = 0 ( RECOLL_SAITO_HHC )  Halliday, eqs. (18)-(19):   d_i = fat w_i ( c_i . n )
	//  forma = 1 ( RECOLL_SAITO_LKR )  Latva-Kokko:                d_i = fat w_i ( c_i . n ) / | c_i |
	//
	//  A amplitude e a mesma nas duas quando p = rho cs^2 :  beta c_R c_B p / cs^2 = beta c_R c_B rho .
	//  Escrevemos a de Halliday, que e a do artigo; a de Latva-Kokko so divide por | c_i |.

	double fat = beta * conc_R * conc_B * p * lattice.one_over_c_s2;

	//------------- Limitador de positividade -----------------------------------------------------//
	//
	//  Vale para as duas formas: um unico lambda <= 1 por sitio, o maior que mantem as duas cores
	//  nao negativas.  Sendo fator global do sitio, soma_i d_i continua zero.

	if ( limita )
	{
		double lambda = 1.0;

		for ( int i = 0; i < nvel; i++ )
		{
			const double *c_i = lattice.c_i + i * dim;

			double ci_n = c_i[0] * n[0] + c_i[1] * n[1] + ( ( dim == 3 ) ? c_i[2] * n[2] : 0.0 );

			if ( forma )
			{
				const double c2 = c_i[0] * c_i[0] + c_i[1] * c_i[1]
				                  + ( ( dim == 3 ) ? c_i[2] * c_i[2] : 0.0 );

				ci_n = ( c2 > 0.0 ) ? ci_n / sqrt ( c2 ) : 0.0;
			}

			const double d = fat * lattice.w[i] * ci_n;

			const double lim = ( d > 0.0 ) ? conc_B * f[i] : conc_R * f[i];

			const double ad = ( d > 0.0 ) ? d : - d;

			if ( ad > lim )
			{
				const double r = ( ad > 0.0 ) ? lim / ad : 0.0;

				if ( r < lambda ) lambda = r;
			}
		}

		if ( lambda < 0.0 ) lambda = 0.0;

		fat = fat * lambda;
	}

	//------------- Recoloracao -------------------------------------------------------------------//

	for ( int i = 0; i < nvel; i++ )
	{
		const double *c_i = lattice.c_i + i * dim;

		double ci_n = c_i[0] * n[0] + c_i[1] * n[1] + ( ( dim == 3 ) ? c_i[2] * n[2] : 0.0 );

		if ( forma )
		{
			const double c2 = c_i[0] * c_i[0] + c_i[1] * c_i[1]
			                  + ( ( dim == 3 ) ? c_i[2] * c_i[2] : 0.0 );

			ci_n = ( c2 > 0.0 ) ? ci_n / sqrt ( c2 ) : 0.0;
		}

		const double d = fat * lattice.w[i] * ci_n;

		f_R[i] = conc_R * f[i] + d;
		f_B[i] = conc_B * f[i] - d;
	}
}

//====================================================================================================================//



//=============================== Colisao do modelo de Saito, em momentos centrais ==================================//
//
//      Eq. (67) do artigo.  Com todos os coeficientes de relaxacao iguais a 1 menos omega_1 -- o da
//      viscosidade --, os momentos centrais pos-colisao sao quase todos CONSTANTES, e so seis
//      precisam ser calculados a partir de f:  k110, k011, k101, k200, k020, k002.
//
//          k100* = Fx/2                       k110* = ( 1 - w ) k110
//          k200* - k020* = ( 1 - w )( k200 - k020 ) + ( 1 - w/2 )( Qx - Qy )
//          k200* - k002* = ( 1 - w )( k200 - k002 ) + ( 1 - w/2 )( Qx - Qz )
//          k200* + k020* + k002* = 3 p + ( Qx + Qy + Qz ) / 2
//          k120* = k102* = cs^2 Fx / 2        k111* = 0
//          k220* = k202* = k022* = p cs^2     k211* = k121* = k112* = 0
//          k122* = cs^4 Fx / 2                k222* = p cs^4
//
//      VELOCIDADE.  Pela eq. (11),  rho u = soma_i f_i c_i + F / 2  ANTES da colisao; e como
//      k100* = Fx/2, depois dela  soma_i f_i* c_i = rho u + F / 2 .  Ou seja: a velocidade fisica
//      pos-colisao e  ( soma_i f_i* c_i - F / 2 ) / rho , e nao o momento cru.  As medidas desta
//      biblioteca fazem essa subtracao quando lattice.ini_force esta preenchido -- ver
//      calc_drop_radius(), em Other_functions.cpp.  Sem isso o que aparece como "corrente espuria"
//      e a metade da forca.
//
//      Q e a correcao da eq. (65),  Q = - 3 div[ ( p - rho cs^2 ) u ] , calculada pelo programa.
//      Em gota estatica ela e nula ( u = 0 ), como toda a maquinaria de invariancia galileana.
//
//      Input: f_R, f_B, forca F, correcao Q, normal n, omega, lattice, parameters
//      Output: f_R, f_B pos-colisao
//
//====================================================================================================================//

#pragma acc routine seq
void coll_saito ( double *f_R, double *f_B, double Fx, double Fy, double Fz,
					double Qx, double Qy, double Qz, double *n, double omega,
					LATTICE lattice, PARAMETERS parameters )
{
	constexpr double eps = 1.0e-30;

	const double cs2 = lattice.c_s2;

	const double rho_R = density ( f_R );
	const double rho_B = density ( f_B );

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	const double p = rho_R * parameters.cs2_R + rho_B * parameters.cs2_B;

	double f[nvel];

	for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

	//------------- Velocidade, eq. (11): inclui a meia-forca ------------------------------------//

	double mx, my, mz;

	momentum ( f, mx, my, mz, lattice );

	const double ux = ( mx + 0.5 * Fx ) / rho;
	const double uy = ( my + 0.5 * Fy ) / rho;
	const double uz = ( mz + 0.5 * Fz ) / rho;

	//------------- Os seis momentos centrais que a colisao precisa ------------------------------//

	double k[27];

	cm_saito ( f, ux, uy, uz, k );

	const double k110 = k[ 1 * 9 + 1 * 3 + 0 ];
	const double k101 = k[ 1 * 9 + 0 * 3 + 1 ];
	const double k011 = k[ 0 * 9 + 1 * 3 + 1 ];

	const double k200 = k[ 2 * 9 + 0 * 3 + 0 ];
	const double k020 = k[ 0 * 9 + 2 * 3 + 0 ];
	const double k002 = k[ 0 * 9 + 0 * 3 + 2 ];

	//------------- Momentos pos-colisao, eq. (67) ------------------------------------------------//

	for ( int m = 0; m < 27; m++ ) k[m] = 0.0;

	k[ 0 * 9 + 0 * 3 + 0 ] = rho;

	k[ 1 * 9 + 0 * 3 + 0 ] = 0.5 * Fx;
	k[ 0 * 9 + 1 * 3 + 0 ] = 0.5 * Fy;
	k[ 0 * 9 + 0 * 3 + 1 ] = 0.5 * Fz;

	k[ 1 * 9 + 1 * 3 + 0 ] = ( 1.0 - omega ) * k110;
	k[ 1 * 9 + 0 * 3 + 1 ] = ( 1.0 - omega ) * k101;
	k[ 0 * 9 + 1 * 3 + 1 ] = ( 1.0 - omega ) * k011;

	//  As tres diagonais saem de duas diferencas e de um traco.

	const double D1 = ( 1.0 - omega ) * ( k200 - k020 ) + ( 1.0 - 0.5 * omega ) * ( Qx - Qy );
	const double D2 = ( 1.0 - omega ) * ( k200 - k002 ) + ( 1.0 - 0.5 * omega ) * ( Qx - Qz );

	const double S  = 3.0 * p + 0.5 * ( Qx + Qy + Qz );

	const double k200n = ( S + D1 + D2 ) / 3.0;

	k[ 2 * 9 + 0 * 3 + 0 ] = k200n;
	k[ 0 * 9 + 2 * 3 + 0 ] = k200n - D1;
	k[ 0 * 9 + 0 * 3 + 2 ] = k200n - D2;

	k[ 1 * 9 + 2 * 3 + 0 ] = 0.5 * cs2 * Fx;			// k120
	k[ 1 * 9 + 0 * 3 + 2 ] = 0.5 * cs2 * Fx;			// k102
	k[ 2 * 9 + 1 * 3 + 0 ] = 0.5 * cs2 * Fy;			// k210
	k[ 0 * 9 + 1 * 3 + 2 ] = 0.5 * cs2 * Fy;			// k012
	k[ 2 * 9 + 0 * 3 + 1 ] = 0.5 * cs2 * Fz;			// k201
	k[ 0 * 9 + 2 * 3 + 1 ] = 0.5 * cs2 * Fz;			// k021

	k[ 1 * 9 + 1 * 3 + 1 ] = 0.0;						// k111

	k[ 2 * 9 + 2 * 3 + 0 ] = p * cs2;					// k220
	k[ 2 * 9 + 0 * 3 + 2 ] = p * cs2;					// k202
	k[ 0 * 9 + 2 * 3 + 2 ] = p * cs2;					// k022

	k[ 1 * 9 + 2 * 3 + 2 ] = 0.5 * cs2 * cs2 * Fx;		// k122
	k[ 2 * 9 + 1 * 3 + 2 ] = 0.5 * cs2 * cs2 * Fy;		// k212
	k[ 2 * 9 + 2 * 3 + 1 ] = 0.5 * cs2 * cs2 * Fz;		// k221

	k[ 2 * 9 + 2 * 3 + 2 ] = p * cs2 * cs2;				// k222

	cm_saito_inv ( k, ux, uy, uz, f );

	//------------- Recoloracao, eqs. (18)-(19) ---------------------------------------------------//

	recolloring_saito ( f, f_R, f_B, n, rho_R, rho_B, p, parameters.recoll,
							parameters.recoll_positivo, parameters.recoll_saito, lattice );
}

//====================================================================================================================//


//=============================== Operador de perturbacao de Liu, Valocchi & Kang ( 2012 ) ===========================//
//
//      Liu, Valocchi & Kang, Phys. Rev. E 85, 046309 (2012), eq. (26):
//
//          Omega_i^(2) = ( A / 2 ) | grad rho^N |  [ w_i ( e_i . grad rho^N )^2 / | grad rho^N |^2  -  B_i ]
//
//      somado sobre as duas cores ( A_R = A_B = A ), o que da o fator A e nao A/2.  Os B_i vem da
//      eq. (29) e formam uma familia de UM parametro livre chi:
//
//          B_0 = - ( 2 + 2 chi ) / ( 3 chi + 12 ) ,   B_1..6 = chi / ( 6 chi + 24 ) ,
//          B_7..18 = 1 / ( 6 chi + 24 )
//
//      Qualquer chi satisfaz as tres condicoes da eq. (28) --  soma B_i = 1/3 ,  soma B_i e_i = 0  e
//      soma B_i e_ia e_ib = delta_ab / 3 -- e sao ELAS, e so elas, que fixam o tensor S.  Por isso a
//      tensao interfacial NAO depende de chi:
//
//          soma_i Omega_i^(2) e_ia e_ib = ( 2 A / 9 ) | grad rho^N | ( n_a n_b - delta_ab )
//
//          S = - tau soma_i Omega_i^(2) e_i e_i = ( sigma / 2 ) | grad rho^N | ( I - n n )
//
//          =>   sigma = ( 4 / 9 ) A tau                                             eq. (30)
//
//      que e a eq. (25) do artigo -- a forma da tensao de superficie do continuum surface force.  E
//      justamente isso que o artigo mostra que o operador de Tolke, eq. (31), NAO consegue: aquele
//      da um S ( eq. 32 ) que nao se escreve na forma da eq. (25).
//
//      O artigo usa chi = 2, que da  B_0 = -1/3 ,  B_1..6 = 1/18 ,  B_7..18 = 1/36  -- isto e,
//      B_i = w_i para i >= 1, e B_0 = -w_0.  interf_tension_RPL(), nesta mesma biblioteca, e o mesmo
//      operador com chi = 1/2 ( B = -2/9, 1/54, 1/27 ): as duas dao a MESMA sigma.
//
//      Conservacao: soma_i Omega_i = A |grad| ( cs^2 - 1/3 ) = 0 e soma_i Omega_i e_i = 0 por
//      paridade -- eqs. (15) e (16).  Massa e momento intactos, para qualquer chi.
//
//      Input: f ( populacao total, pos-BGK ), A, | grad rho^N |, grad rho^N, chi, lattice
//      Output: f com a perturbacao somada
//
//====================================================================================================================//

#pragma acc routine seq
void interf_tension_LVK ( double *f, double A, double mod_grad, double *grad, double chi, LATTICE lattice )
{
	constexpr double eps = 1.0e-30;

	if ( mod_grad <= eps || A == 0.0 ) return;

	const double den = 6.0 * chi + 24.0;

	const double B_0 = - ( 2.0 + 2.0 * chi ) / ( 3.0 * chi + 12.0 );
	const double B_1 =           chi / den;
	const double B_7 =           1.0 / den;

	const double inv_mod2 = 1.0 / ( mod_grad * mod_grad );

	for ( int i = 0; i < nvel; i++ )
	{
		const double *c_i = lattice.c_i + i * dim;

		const double ci_grad = c_i[0] * grad[0] + c_i[1] * grad[1]
		                       + ( ( dim == 3 ) ? c_i[2] * grad[2] : 0.0 );

		const double B_i = ( i == 0 ) ? B_0 : ( ( i < 7 ) ? B_1 : B_7 );

		f[i] = f[i] + A * mod_grad * ( lattice.w[i] * ci_grad * ci_grad * inv_mod2 - B_i );
	}
}

//====================================================================================================================//



//=============================== Recoloracao de Latva-Kokko com equilibrio de alpha variavel ========================//
//
//      Liu, Valocchi & Kang, eqs. (35)-(37):
//
//          f_R,i = ( rho_R / rho ) f_i*  +  beta ( rho_R rho_B / rho^2 ) cos(phi_i) f_i^eq |_{u=0}
//          f_B,i = ( rho_B / rho ) f_i*  -  beta ( rho_R rho_B / rho^2 ) cos(phi_i) f_i^eq |_{u=0}
//
//          cos(phi_i) = e_i . grad rho^N / ( | e_i | | grad rho^N | )
//
//      A diferenca para recolloring(), que ja existe nesta biblioteca, e UMA so: aqui o equilibrio
//      em repouso e
//
//          f_i^eq |_{u=0} = rho_R phi_i^R + rho_B phi_i^B ,       phi_i^k da eq. (9)
//
//      e nao  rho w_i .  As duas coincidem exatamente quando alpha_R = alpha_B = 1/3, que e o caso
//      de densidades iguais; para  alpha != 1/3  -- que e como o modelo alcanca razao de densidades
//      -- a distribuicao em repouso deixa de ser w_i e a recoloracao tem de acompanhar, senao ela
//      injeta massa onde o equilibrio nao a tem.
//
//      SINAL.  grad aqui vem dos mediadores, isto e, vale  - grad rho^N .  Por isso os sinais estao
//      trocados em relacao as eqs. (35)-(36), exatamente como em recolloring().
//
//      Input: f ( total, pos-perturbacao ), grad, |grad|, rho_R, rho_B, beta, alpha_R, alpha_B
//      Output: f_R, f_B
//
//====================================================================================================================//

#pragma acc routine seq
void recolloring_LVK ( double *f, double *f_R, double *f_B, double *grad, double mod_grad,
						double rho_R, double rho_B, double beta, double alfa_R, double alfa_B,
						LATTICE lattice )
{
	constexpr double eps = 1.0e-30;

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	const double conc_R = rho_R / rho;
	const double conc_B = 1.0 - conc_R;

	const double fator = beta * conc_R * conc_B;

	//  f_i^eq |_{u=0} = rho_R phi_i^R + rho_B phi_i^B , nos tres grupos do D3Q19

	const double eq_0 = rho_R * alfa_R + rho_B * alfa_B;
	const double eq_1 = ( rho_R * ( 1.0 - alfa_R ) + rho_B * ( 1.0 - alfa_B ) ) / 12.0;
	const double eq_7 = ( rho_R * ( 1.0 - alfa_R ) + rho_B * ( 1.0 - alfa_B ) ) / 24.0;

	f_R[0] = conc_R * f[0];
	f_B[0] = conc_B * f[0];

	for ( int i = 1; i < nvel; i++ )
	{
		const double *c_i = lattice.c_i + i * dim;

		const double c_2 = c_i[0] * c_i[0] + c_i[1] * c_i[1]
		                   + ( ( dim == 3 ) ? c_i[2] * c_i[2] : 0.0 );

		const double ci_grad = c_i[0] * grad[0] + c_i[1] * grad[1]
		                       + ( ( dim == 3 ) ? c_i[2] * grad[2] : 0.0 );

		const double cos_phi = ( mod_grad > eps ) ? ci_grad / ( mod_grad * sqrt ( c_2 ) ) : 0.0;

		const double f_eq_0 = ( i < 7 ) ? eq_1 : eq_7;

		f_R[i] = conc_R * f[i] - fator * f_eq_0 * cos_phi;
		f_B[i] = conc_B * f[i] + fator * f_eq_0 * cos_phi;
	}

	//  eq_0 nao entra: o repouso nao tem direcao.  Fica declarado so para deixar claro de onde
	//  vem a normalizacao -- soma_i ( phi^R_i rho_R + phi^B_i rho_B ) = rho .

	( void ) eq_0;
}

//====================================================================================================================//



//=============================== Modelo de Liu, Valocchi & Kang ( 2012 ), D3Q19 =====================================//
//
//      Phys. Rev. E 85, 046309 (2012).  Tres subetapas, eq. (2):
//
//          Omega_i^k = ( Omega^k )^(3) [ ( Omega^k )^(1) + ( Omega^k )^(2) ]
//
//          (1) BGK em direcao ao equilibrio da eq. (8), com o phi_i^k da eq. (9);
//          (2) perturbacao da eq. (26), que gera a tensao interfacial;
//          (3) recoloracao de Latva-Kokko, eqs. (35)-(36).
//
//      Tres coisas distinguem este modelo dos outros color-gradient desta arvore:
//
//      1) O EQUILIBRIO.  f_i^{k,eq} = rho_k [ phi_i^k + w_i ( 3 e.u + 4.5 (e.u)^2 - 1.5 u^2 ) ] ,
//         com phi_i^k = alpha_k , (1-alpha_k)/12 , (1-alpha_k)/24 .  Cada fluido ganha sua propria
//         velocidade do som, cs_k^2 = ( 1 - alpha_k ) / 2 , e portanto sua propria pressao,
//         p_k = rho_k ( 1 - alpha_k ) / 2 .  Continuidade da pressao numa interface plana exige
//
//             lambda = rho_R / rho_B = ( 1 - alpha_B ) / ( 1 - alpha_R )                eq. (11)
//
//         que e como o modelo chega a razao de densidades ate 1000.  Com alpha = 1/3 o equilibrio
//         recai no de sempre ( phi_i = w_i , cs^2 = 1/3 ).  Reutiliza dist_eq_LYL com hi_order = 0.
//
//      2) A PERTURBACAO vem do continuum surface force, nao de um ansatz.  Ver interf_tension_LVK.
//         O ganho pratico anunciado pelo artigo e sobre o operador de Tolke, eq. (31), que nao
//         reproduz a forma correta do tensor.
//
//      3) O TEMPO DE RELAXACAO e interpolado LINEARMENTE em rho^N atraves da interface, eq. (13):
//
//             tau = ( 1 + rho^N ) / 2 tau_R  +  ( 1 - rho^N ) / 2 tau_B
//
//         em vez do degrau  ( conc_R > conc_B ? tau_R : tau_B )  usado em coll_RPL_D3Q19.  Com
//         tau_R = tau_B nao faz diferenca nenhuma; com viscosidades diferentes, faz.
//
//      A colisao BGK e feita na populacao TOTAL, com o tau interpolado e o equilibrio total
//      f_i^eq = soma_k f_i^{k,eq} .  E o que as eqs. (13) e (35) pressupoem: a eq. (35) recolore
//      f_i*, que ja e a distribuicao total pos-perturbacao.
//
//      A aceleracao de corpo entra pelo esquema de Guo, como no resto da biblioteca.  ATENCAO: se
//      voce usar forca de corpo, a velocidade fisica pos-colisao e  soma_i f_i e_i / rho - F / 2 ,
//      nao o momento cru -- ver a correcao da meia-forca em calc_drop_radius().  Na lei de Laplace
//      a aceleracao e nula e a questao nao aparece.
//
//      O gradiente do campo de fase vem dos mediadores, que entregam  - grad rho^N  com o estencil
//      isotropico de 4a ordem  ( 3 / c^2 ) soma_i w_i rho^N( x + e_i ) e_i  -- exatamente a eq. (33)
//      do artigo.  A fase emitida tem de ser rho^N ( simetrica, em [-1,1] ), nao phi.
//
//      sigma = ( 4 / 9 ) A tau ,  com  A = parameters.A_fact .
//
//      Input: f_R, f_B, mediadores, aceleracao de corpo, lattice, parameters
//      Output: f_R, f_B pos-colisao
//
//====================================================================================================================//

#pragma acc routine seq
void coll_LVK_D3Q19 ( double *f_R, double *f_B, double *f_m, double gx, double gy, double gz,
						LATTICE lattice, PARAMETERS parameters )
{
	constexpr double eps = 1.0e-30;

	const double rho_R = density ( f_R );
	const double rho_B = density ( f_B );

	const double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	const double rho_N = ( rho_R - rho_B ) / rho;

	//------------- Gradiente do campo de fase, pelos mediadores ---------------------------------//

	double grad_x, grad_y, grad_z;

	momentum ( f_m, grad_x, grad_y, grad_z, lattice );

	double grad[3] = { grad_x, grad_y, grad_z };

	const double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

	//------------- Populacao total e velocidade -------------------------------------------------//

	double f[nvel];

	for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

	double vx, vy, vz, rho_f;

	calcula ( f, vx, vy, vz, rho_f, lattice );

	//------------- (1) BGK, eq. (3), com o tau interpolado da eq. (13) ---------------------------//

	const double tau = 0.5 * ( 1.0 + rho_N ) * parameters.tau_R
	                 + 0.5 * ( 1.0 - rho_N ) * parameters.tau_B;

	//  Guo: desloca a velocidade de a/2 no equilibrio e adiciona o termo fonte.

	const double vx_alt = vx + 0.5 * gx;
	const double vy_alt = vy + 0.5 * gy;
	const double vz_alt = vz + 0.5 * gz;

	double S[nvel];

	source ( gx * rho, gy * rho, gz * rho, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );

	//  Equilibrio da eq. (8): soma das duas cores.  dist_eq_LYL com hi_order = 0 e exatamente ele.

	double f_eq_R[nvel], f_eq_B[nvel];

	dist_eq_LYL ( f_eq_R, vx_alt, vy_alt, vz_alt, rho_R, parameters.alfa_R, 0.0, lattice );
	dist_eq_LYL ( f_eq_B, vx_alt, vy_alt, vz_alt, rho_B, parameters.alfa_B, 0.0, lattice );

	const double one_over_tau = 1.0 / tau;

	for ( int i = 0; i < nvel; i++ )
	{
		f[i] = f[i] + ( f_eq_R[i] + f_eq_B[i] - f[i] ) * one_over_tau + S[i];
	}

	//------------- (2) Perturbacao, eq. (26) ----------------------------------------------------//

	interf_tension_LVK ( f, parameters.A_fact, mod_grad, grad, parameters.chi_LVK, lattice );

	//------------- (3) Recoloracao, eqs. (35)-(36) ----------------------------------------------//

	recolloring_LVK ( f, f_R, f_B, grad, mod_grad, rho_R, rho_B, parameters.recoll,
						parameters.alfa_R, parameters.alfa_B, lattice );
}

//====================================================================================================================//


//=============================== Spencer-Halliday-Care 2010 -- versao H do relatorio LMPT ===========================//
//
//      Mattila, Bresolin, Philippi & Emerich, relatorio LMPT/UFSC ( 2013 ), eqs. (1.2) e (1.6):
//
//          Pi^srf_ab = ( rho / tau ) A_H | grad w^rb | ( n_a n_b - cs^2 delta_ab )
//          Gamma_i   = Pi^srf_ab K_i,ab ,      K_i,ab = w_i ( c_ia c_ib - cs^2 delta_ab ) / ( 2 cs^4 )
//
//      que da, em forma fechada,
//
//          Gamma_i = coef w_i [ ( n.c_i )^2 - cs^2 c_i^2 - cs^2 + D cs^4 ] / ( 2 cs^4 )
//
//      Essa base K_i e a que faz  soma_i Gamma_i c_ia c_ib = Pi^srf_ab  exatamente, e ao mesmo tempo
//      soma_i Gamma_i = 0  e  soma_i Gamma_i c_i = 0 .
//
//      Diferenca em relacao a coll_SHC_2010: aquela impoe  ( 2/3 )( n_a n_b - delta_ab )  em vez de
//      ( n_a n_b - cs^2 delta_ab ), isto e, amplitude anisotropica 2/3 da correta e parte isotropica
//      tres vezes maior.  Ver MODELO_SHC.md, na pasta Spencer-Halliday-Care.
//
//      A analise de Chapman-Enskog do relatorio ( secao 1.4.3 ) mostra que esta versao -- com
//      | grad w^rb | -- e livre dos termos que quebram a invariancia galileana na versao G, que usa
//      w^r w^b no lugar do modulo do gradiente.
//
//      Input: f_R, f_B, mediadores, aceleracao de corpo, lattice, parameters
//      Output: f_R, f_B pos-colisao
//
//====================================================================================================================//

#pragma acc routine seq
void coll_SHC_H ( double* f_R, double* f_B, double* f_m, double gx, double gy, double gz,
					LATTICE lattice, PARAMETERS parameters )
{
	constexpr double eps = 1.0e-30;

	double rho_R = density ( f_R );
	double rho_B = density ( f_B );

	double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	double conc_R = rho_R / rho;
	double conc_B = 1.0 - conc_R;

	//------------- Gradiente da cor, pelos mediadores -------------------------------------------//

	double grad_x, grad_y, grad_z;

	momentum ( f_m, grad_x, grad_y, grad_z, lattice );

	double grad[3] = { grad_x, grad_y, grad_z };

	double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

	double nx, ny, nz;

	unit_vector ( grad_x, grad_y, grad_z, nx, ny, nz );

	//------------------- Colisao monofasica da mistura ------------------------------------------//

	if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
	else                   parameters.tau = parameters.tau_B;

	double f[nvel];

	for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

	coll_BGK ( f, gx, gy, gz, lattice, parameters );

	//------------------- Termo de interface ------------------------------------------------------//

	double A = parameters.A_fact;

	if ( mod_grad > eps && A != 0.0 )
	{
		const double cs2 = lattice.c_s2;

		const double coef = ( rho / parameters.tau ) * A * mod_grad / ( 2.0 * cs2 * cs2 );

		const double termo_const = - cs2 + ( double ) dim * cs2 * cs2;

		for ( int i = 0; i < nvel; i++ )
		{
			const double *c_i = lattice.c_i + i * dim;

			const double n_ci = c_i[0] * nx + c_i[1] * ny + ( dim == 3 ? c_i[2] * nz : 0.0 );

			const double c_2 = dot_product ( ( double* ) c_i, ( double* ) c_i );

			f[i] = f[i] + coef * lattice.w[i] * ( n_ci * n_ci - cs2 * c_2 + termo_const );
		}
	}

	//------------------- Recoloracao de Latva-Kokko ----------------------------------------------//

	recolloring_SHC ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, parameters.recoll, lattice, parameters );
}

//====================================================================================================================//



//=============================== Tensao interfacial Dsigma'_ab ( Spencer-Halliday-Care 2010 ) =======================//
//
//      Eq. (5) do artigo:      Dsigma'_ab = alpha0 | grad rho^N |  n_a n_b
//
//      Guarda as seis componentes independentes em T, na ordem  xx, yy, zz, xy, xz, yz .
//      alpha0 = parameters.A_fact.
//
//      Note que o artigo NAO poe rho nem 1/tau aqui: Dsigma' e uma tensao, e a forca sai da sua
//      divergencia ( eq. 22 ).  Quem faz a divergencia e o programa, porque ela nao e local.
//
//====================================================================================================================//

//  Nucleo: recebe o gradiente de cor JA calculado, na convencao dos mediadores
//  ( grad = - grad rho^N , que aponta para fora do vermelho, como n^ ).  Existe para que o
//  programa possa alimenta-lo com um estencil mais isotropico que o dos mediadores.

#pragma acc routine seq
void stress_SHC_g ( double rho_R, double rho_B, double grad_x, double grad_y, double grad_z,
					double* T, PARAMETERS parameters, bool peso_analitico )
{
	constexpr double eps = 1.0e-30;

	for ( int k = 0; k < 6; k++ ) T[k] = 0.0;

	double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

	if ( mod_grad <= eps ) return;

	double nx, ny, nz;

	unit_vector ( grad_x, grad_y, grad_z, nx, ny, nz );

	//  PESO.  O artigo da duas formas equivalentes no equilibrio, eqs. (8) e (9):
	//
	//      (8)   w = - alpha0 ( n^ . grad rho^N ) = - alpha0 | grad rho^N |
	//      (9)   w = - alpha0 beta [ 1 - ( rho^N )^2 ] = - 4 alpha0 beta rho_R rho_B / rho^2
	//
	//  As duas coincidem para o perfil  rho^N = tanh( beta s )  que a recoloracao de Latva-Kokko
	//  produz.  Usamos a (9): ela e LOCAL e ANALITICA -- e o que o resumo do artigo chama de "an
	//  interface-inducing force distribution which is analytic".  Isso importa porque a forca e a
	//  DIVERGENCIA deste tensor: com a forma (8), w ja e um gradiente numerico, e a divergencia
	//  vira uma segunda derivada discreta de um perfil de tres sitios -- ruidosa, e a corrente
	//  espuria sobe.  Medido em R = 40, sigma = 4.3e-03, tau = 1, caixa 121x121:
	//
	//      forma (8), | grad rho^N |     |u|max = 2.20e-04
	//      forma (9), analitica          |u|max = ver LAPLACE_SHC.md
	//
	//  SINAL: pela eq. (4), n^ = - grad(rho^N)/|grad(rho^N)| aponta para FORA do vermelho, e as duas
	//  formas dao w < 0.  Importa: pela eq. (56), F = ( n^ . grad w ) n^ + K w n^ , e numa gota
	//  K = +1/R com n^ para fora -- so com w < 0 a forca de curvatura aponta para DENTRO.
	//
	//  beta e parameters.recoll, o mesmo da recoloracao: as duas coisas descrevem a mesma interface.

	const double w = peso_analitico
	               ? - 4.0 * parameters.A_fact * parameters.recoll * rho_R * rho_B / ( rho * rho )
	               : - parameters.A_fact * mod_grad;

	T[0] = w * nx * nx;			// xx
	T[1] = w * ny * ny;			// yy
	T[2] = w * nz * nz;			// zz
	T[3] = w * nx * ny;			// xy
	T[4] = w * nx * nz;			// xz
	T[5] = w * ny * nz;			// yz
}

//  Rota padrao: o gradiente vem dos mediadores.

#pragma acc routine seq
void stress_SHC ( double* f_R, double* f_B, double* f_m, double* T, LATTICE lattice, PARAMETERS parameters,
					bool peso_analitico )
{
	double grad_x, grad_y, grad_z;

	momentum ( f_m, grad_x, grad_y, grad_z, lattice );

	stress_SHC_g ( density ( f_R ), density ( f_B ), grad_x, grad_y, grad_z, T, parameters, peso_analitico );
}

//====================================================================================================================//



//=============================== Campo do normal e do peso ( rota da eq. 56 do SHC 2010 ) ===========================//
//
//      A eq. (22),  F_a = d_b ( Dsigma' )_ab , com  Dsigma'_ab = w n^_a n^_b , exige DUAS derivadas
//      discretas encadeadas: o gradiente ( pelos mediadores ) que da n^ e w, e depois a divergencia
//      do produto.  Compor dois estenceis  w_i c_i  nao e uma segunda derivada bem construida, e o
//      erro aparece como corrente espuria.
//
//      A eq. (56) do artigo abre a divergencia:
//
//          F = ( n^ . grad w ) n^  +  K w n^  +  w ( n^ . grad ) n^ ,      K = div n^        (56)
//
//      O ultimo termo e a parte TANGENCIAL de grad|grad rho^N| dividida por |grad rho^N|; ele se
//      anula identicamente para uma interface de simetria radial ( gota, cilindro ), que e o caso
//      da lei de Laplace.  O artigo o omite.  Aqui tambem.
//
//      O ganho e este: com o peso ANALITICO da eq. (9),
//
//          w = - alpha0 beta [ 1 - ( rho^N )^2 ]      =>      grad w = 2 alpha0 beta rho^N grad rho^N
//
//      e, como  n^ = - grad rho^N / | grad rho^N | ,
//
//          n^ . grad w = - 2 alpha0 beta rho^N | grad rho^N |                              LOCAL
//
//      isto e: o primeiro termo sai FECHADO dos mediadores, sem estencil nenhum.  Sobra so a
//      curvatura  K = div n^ , UMA derivada primeira de um campo de modulo 1 -- que e liso e O(1),
//      ao contrario de  Dsigma' , que e um pico de tres sitios.
//
//      Esta funcao so prepara o campo; quem faz o estencil de  div n^  e o programa, porque le os
//      vizinhos.  Grava, por sitio:
//
//          C[0..2] = n^          C[3] = w          C[4] = n^ . grad w        C[5] = | grad rho^N |
//
//      CORTE.  unit_vector() normaliza qualquer vetor nao nulo: longe da interface  grad rho^N  e
//      so ruido de arredondamento, e n^ viraria um campo de modulo 1 apontando para todo lado --
//      cuja divergencia e enorme.  Por isso zeramos n^ ( e com ele a forca, que e toda proporcional
//      a n^ ) onde  | grad rho^N | < corte .  Como F ~ n^ local, o corte nao deixa residuo: fora da
//      interface a forca e exatamente zero, e nao "quase zero".
//
//      Input: f_R, f_B, mediadores, lattice, parameters, peso_analitico, corte
//      Output: C[6]
//
//====================================================================================================================//

#pragma acc routine seq
void campo_normal_SHC_g ( double rho_R, double rho_B, double grad_x, double grad_y, double grad_z,
						double* C, PARAMETERS parameters, bool peso_analitico, double corte )
{
	constexpr double eps = 1.0e-30;

	for ( int k = 0; k < 6; k++ ) C[k] = 0.0;

	double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

	C[5] = mod_grad;

	if ( mod_grad < corte ) return;

	double nx, ny, nz;

	unit_vector ( grad_x, grad_y, grad_z, nx, ny, nz );

	C[0] = nx;
	C[1] = ny;
	C[2] = nz;

	//  Peso: eqs. (8) e (9) do artigo, as mesmas duas formas de stress_SHC.  w < 0 sempre.

	const double rho_N = ( rho_R - rho_B ) / rho;

	C[3] = peso_analitico
	     ? - parameters.A_fact * parameters.recoll * ( 1.0 - rho_N * rho_N )
	     : - parameters.A_fact * mod_grad;

	//  n^ . grad w .  Com o peso analitico e fechado; com o peso da eq. (8) nao e ( grad|grad rho^N|
	//  nao e local ), e o programa tem de calcula-lo por estencil sobre C[3].  Sinal: os mediadores
	//  entregam  soma_i M_i c_i = - grad rho^N , entao  ( grad_x, grad_y, grad_z )  JA e  -grad rho^N,
	//  o mesmo sentido de n^, e  n^ . grad rho^N = - mod_grad .

	C[4] = peso_analitico ? ( - 2.0 * parameters.A_fact * parameters.recoll * rho_N * mod_grad ) : 0.0;
}

//  Rota padrao: o gradiente vem dos mediadores.

#pragma acc routine seq
void campo_normal_SHC ( double* f_R, double* f_B, double* f_m, double* C, LATTICE lattice,
						PARAMETERS parameters, bool peso_analitico, double corte )
{
	double grad_x, grad_y, grad_z;

	momentum ( f_m, grad_x, grad_y, grad_z, lattice );

	campo_normal_SHC_g ( density ( f_R ), density ( f_B ), grad_x, grad_y, grad_z, C, parameters,
						peso_analitico, corte );
}

//====================================================================================================================//



//=============================== Spencer-Halliday-Care 2010, como publicado =========================================//
//
//      O artigo ( secao II C, logo apos a eq. 21 ) testa impor Dsigma' como termo fonte direto e
//      DESCARTA essa rota, porque "the level of the interfacial microcurrent increases
//      significantly".  A eq. (59) mostra que aquela rota reproduz o operador de Gunstensen (1991).
//      O que o artigo usa e a DIVERGENCIA do tensor, como forca de corpo:
//
//          F_a = d_b ( Dsigma' )_ab                                             (22)
//
//      entrando pelo forcamento de Guo, com  B = ( 1 - 1/2tau ) F  (15)  e a correcao da velocidade
//      rho u = soma_i f_i c_i + dt F / 2  (17).  E exatamente o que coll_BGK() desta biblioteca ja
//      faz quando recebe uma aceleracao.
//
//      A divergencia NAO e local: quem a calcula e o programa, lendo Dsigma' dos vizinhos, e passa
//      o resultado aqui em ( Fx, Fy, Fz ).  ( gx, gy, gz ) e a aceleracao de corpo, se houver.
//
//      Input: f_R, f_B, mediadores, forca interfacial F, aceleracao de corpo, lattice, parameters
//      Output: f_R, f_B pos-colisao
//
//====================================================================================================================//

//  Nucleo: o gradiente de cor da recoloracao entra pronto, para que o programa possa usar um
//  estencil mais isotropico que o dos mediadores.  Mesma convencao de sinal ( grad = - grad rho^N ).

#pragma acc routine seq
void coll_SHC_2010_guo_g ( double* f_R, double* f_B, double grad_x, double grad_y, double grad_z,
						double Fx, double Fy, double Fz,
						double gx, double gy, double gz, LATTICE lattice, PARAMETERS parameters )
{
	constexpr double eps = 1.0e-30;

	double rho_R = density ( f_R );
	double rho_B = density ( f_B );

	double rho = rho_R + rho_B;

	if ( rho <= eps ) return;

	double conc_R = rho_R / rho;
	double conc_B = 1.0 - conc_R;

	double grad[3] = { grad_x, grad_y, grad_z };

	double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

	//------------------- Colisao monofasica com forcamento de Guo --------------------------------//

	if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
	else                   parameters.tau = parameters.tau_B;

	double f[nvel];

	for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

	//  coll_BGK recebe ACELERACAO: desloca a velocidade de a/2 no equilibrio e na fonte, que e a
	//  receita de Guo das eqs. (15) e (17) do artigo.

	coll_BGK ( f, gx + Fx / rho, gy + Fy / rho, gz + Fz / rho, lattice, parameters );

	//------------------- Recoloracao de Latva-Kokko ----------------------------------------------//

	recolloring_SHC ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, parameters.recoll, lattice, parameters );
}

//  Rota padrao: o gradiente vem dos mediadores.

#pragma acc routine seq
void coll_SHC_2010_guo ( double* f_R, double* f_B, double* f_m, double Fx, double Fy, double Fz,
						double gx, double gy, double gz, LATTICE lattice, PARAMETERS parameters )
{
	double grad_x, grad_y, grad_z;

	momentum ( f_m, grad_x, grad_y, grad_z, lattice );

	coll_SHC_2010_guo_g ( f_R, f_B, grad_x, grad_y, grad_z, Fx, Fy, Fz, gx, gy, gz, lattice, parameters );
}

//====================================================================================================================//



#pragma acc routine seq
void coll_SHC_2010 ( double* f_R, double* f_B, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry,
						LATTICE lattice, PARAMETERS parameters )
{
	
	//------------- Calcula densidades e concentrações -------------------------------------------//

	double rho_R = density ( f_R );
	double rho_B = density ( f_B );
	
	double rho = rho_R + rho_B;
	
	double conc_R = rho_R / rho; 
	double conc_B = 1. - conc_R;
	
	//------------- Calcula a direção do gradiente ------------------------------------------------/
    
	double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
        
    double n[3] = {grad_uni_x, grad_uni_y, grad_uni_z};
   
    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );
   
    //------------------- Colisão monofásica -----------------------------------------------------//
    
    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    
    else parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
    coll_BGK ( f, gx, gy, gz, lattice, parameters );

	//------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//
	
	double fat_int = parameters.A_fact;
	
	double tau = parameters.tau;

    double fator = rho * ( 1. / tau ) * fat_int * mod_grad * ( 3.0 );

    interf_tension_SHC ( f, fator, n, lattice );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {grad_x, grad_y, grad_z};    

	double beta = parameters.recoll;
	
    recolloring_SHC ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice, parameters );
}

//====================================================================================================================//




//=============================== Collision step using the Reis-Phillips/Leclaire model (immiscible) ================//
//
//      Input: red and blue distribution functions, site coordinates, external acceleration,
//             geometry, lattice and parameters
//      Output: post-collisional red and blue distribution functions
//
//      Notes:
//      1) The gradient and recoloring conventions are kept identical to coll_SHC_2010.
//      2) parameters.A_fact is interpreted as the RPL perturbation amplitude A.
//         With the usual D3Q19 normalization, sigma_lu ~= (2/9) A.
//      3) The perturbation is applied after the hydrodynamic BGK collision and before recoloring.
//
//====================================================================================================================//

#pragma acc routine seq
void coll_RPL_D3Q19 ( double* f_R, double* f_B, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry,
                        LATTICE lattice, PARAMETERS parameters )
{
    constexpr double eps = 1.0e-30;

    //------------- Calcula densidades e concentrações -------------------------------------------//

    double rho_R = density ( f_R );
    double rho_B = density ( f_B );

    double rho = rho_R + rho_B;

    if ( rho <= eps ) return;

    double conc_R = rho_R / rho;
    double conc_B = 1.0 - conc_R;

    //------------- Calcula o gradiente do campo de fase -----------------------------------------//

    double grad_x, grad_y, grad_z;

    gradient ( lattice.ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );

    // Mantém a mesma convenção de sinal usada em coll_SHC_2010 e na recoloração existente.
    grad_x = -grad_x;
    grad_y = -grad_y;
    grad_z = -grad_z;

    double grad[3] = { grad_x, grad_y, grad_z };

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    //------------------- Colisão monofásica da mistura ------------------------------------------//

    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    else                   parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

    coll_BGK ( f, gx, gy, gz, lattice, parameters );

    //------------------- Imposição da tensão interfacial RPL ------------------------------------//

    double A = parameters.A_fact;

    if ( mod_grad > eps && A != 0.0 )
    {
        interf_tension_RPL ( f, A, mod_grad, grad, lattice );
    }

    //------------------- Etapa de recoloração (Latva-Kokko) -------------------------------------//

    double beta = parameters.recoll;

    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//



//=============================== Collision step using the model H (immicisble) ======================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_recoll_BGK ( double* f_R, double* f_B, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry,
						LATTICE lattice, PARAMETERS parameters )
{
	
	double rho_R = density ( f_R );
	double rho_B = density ( f_B );
	
	double rho = rho_R + rho_B;
	
	double conc_R = rho_R / rho; 
	double conc_B = 1. - conc_R;
	
	//------------- Calcula a direção do gradiente ------------------------------------------------/
    
	double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
        
    double n[3] = {grad_uni_x, grad_uni_y, grad_uni_z};
   
    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );
    
    //------------------- Colisão monofásica -----------------------------------------------------//
    
    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    
    else parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
    coll_BGK ( f, gx, gy, gz, lattice, parameters );

	//------------------- Imposição da tensão interfacial (model H) ------------------------------//
	
	interf_tension ( f, parameters.A_fact, parameters.tau, n, mod_grad, lattice );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {grad_x, grad_y, grad_z};    

	double beta = parameters.recoll;
	
    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//



//=============================== Collision step using the model H (immicisble) ======================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_recoll_ExtReg ( double* f_R, double* f_B, int x, int y, int z, double gx, double gy, double gz, 
							GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	
	double rho_R = density ( f_R );
	double rho_B = density ( f_B );
	
	double rho = rho_R + rho_B;
	
	double conc_R = rho_R / rho; 
	double conc_B = 1. - conc_R;
	
	//------------- Calcula a direção do gradiente ------------------------------------------------/
    
	double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
        
    double n[3] = {grad_uni_x, grad_uni_y, grad_uni_z};
   
    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );
    
    //------------------- Colisão monofásica -----------------------------------------------------//
    
    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    
    else parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
    coll_Ext_Reg ( f, gx, gy, gz, lattice, parameters );

	//------------------- Imposição da tensão interfacial (model H) ------------------------------//
	
	interf_tension ( f, parameters.A_fact, parameters.tau, n, mod_grad, lattice );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {grad_x, grad_y, grad_z};    

	double beta = parameters.recoll;
	
    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//



//=============================== Collision step using the model H (immicisble) ======================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_recoll_GR ( double* f_R, double* f_B, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry,
						LATTICE lattice, PARAMETERS parameters )
{
	
	double rho_R = density ( f_R );
	double rho_B = density ( f_B );
	
	double rho = rho_R + rho_B;
	
	double conc_R = rho_R / rho; 
	double conc_B = 1. - conc_R;
	
	//------------- Calcula a direção do gradiente ------------------------------------------------/
    
	double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
        
    double n[3] = {grad_uni_x, grad_uni_y, grad_uni_z};
   
    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );
    
    //------------------- Colisão monofásica -----------------------------------------------------//
    
    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    
    else parameters.tau = parameters.tau_B;

    double f[nvel];

    for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
    coll_BGK ( f, gx, gy, gz, lattice, parameters );

	//------------------- Imposição da tensão interfacial (model H) ------------------------------//
	
	interf_tension_GR ( f, parameters.A_fact, parameters.tau, n, mod_grad, lattice );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {grad_x, grad_y, grad_z};    

	double beta = parameters.recoll;
	
    recolloring ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice );
}

//====================================================================================================================//



//=============================== Collision step using the Spencer, Halliday & Care model (immicisble) ===============//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_SHC_diff ( double* f_R, double* f_B, int x, int y, int z, double gx, double gy, double gz, GEOMETRY geometry,
						LATTICE lattice, PARAMETERS parameters )
{
	
	double rho_R = density ( f_R );
	double rho_B = density ( f_B );
	
	double rho = rho_R + rho_B;
	
	double conc_R = rho_R / rho; 
	double conc_B = 1. - conc_R;
	
	//------------- Calcula a direção do gradiente ------------------------------------------------/
    
	double* ini_psi = lattice.ini_psi;
    
    double grad_x, grad_y, grad_z;
    
    gradient ( ini_psi, x, y, z, grad_x, grad_y, grad_z, lattice, geometry );
    
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;
    
    double grad_uni_x, grad_uni_y, grad_uni_z;

    unit_vector ( grad_x, grad_y, grad_z, grad_uni_x, grad_uni_y, grad_uni_z );    
        
    double n[3] = {grad_uni_x, grad_uni_y, grad_uni_z};
   
    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );
    
    //------------------- Atração R-R ------------------------------------------------------------//
    
    double rho_std = parameters.rho_ini_R;
    
    double alpha_2 = ( ( rho - rho_std ) / rho_std ) * ( ( rho - rho_std ) / rho_std );
    
    double fact_dens = exp( - 10 * alpha_2 );
    
    double fat_int_R = parameters.A_fact_R;
   
    double acc_x = gx + fat_int_R * fact_dens * grad_uni_x;
    double acc_y = gy + fat_int_R * fact_dens * grad_uni_y;
    double acc_z = gz + fat_int_R * fact_dens * grad_uni_z;
    
    //------------------- Colisão monofásica -----------------------------------------------------//
    
    if ( conc_R > conc_B ) parameters.tau = parameters.tau_R;
    
    else parameters.tau = parameters.tau_B;
    
    double f[nvel];
    
    if ( conc_R < conc_B )
    {		
		for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
		coll_BGK ( f, gx, gy, gz, lattice, parameters );
	}
	else 
	{
		for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];
    
		coll_BGK ( f, acc_x, acc_y, acc_z, lattice, parameters );
	}

	//------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//
	
	double fat_int = parameters.A_fact;
	
	double tau = parameters.tau;

    double fator = rho * ( 1. / tau ) * fat_int * mod_grad * ( 3.0 );

    interf_tension_SHC ( f, fator, n, lattice );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//
    
    double grad[3] = {grad_x, grad_y, grad_z};    

	double beta = parameters.recoll;
	
    recolloring_SHC ( f, f_R, f_B, grad, mod_grad, conc_R, conc_B, rho, beta, lattice, parameters );
}

//====================================================================================================================//


//====================================================================================================================//
//====================================================================================================================//
//
//                    MODELO MONTESSORI, HEGELE JR. & LAURICELLA (2024) -- FLUIDOS IMISCIVEIS
//
//  Referencias:
//      [MHL]     Montessori, Hegele Jr. & Lauricella, "A high-performance lattice Boltzmann model for
//                multicomponent turbulent jet simulations", arXiv:2403.15773 (2024).
//      [Talita]  Silva, T. L., "Modelagem e simulacao de jato multicomponente pelo metodo lattice
//                Boltzmann", dissertacao de mestrado, UDESC (2026), secoes 3.1.4 e 3.1.5.
//
//  O modelo substitui a estrutura de duas cores (f_R, f_B + perturbacao interfacial + recoloracao)
//  por:
//
//      1) uma unica populacao hidrodinamica f, com colisao regularizada de segunda ordem
//         (equilibrio truncado em 2a ordem e nao-equilibrio reconstruido do tensor de fluxo de
//         momento) -- eqs. 17 a 19 de [Talita];
//
//      2) um campo de fase phi transportado por uma segunda populacao g, com relaxacao unitaria,
//         equilibrio de primeira ordem e um termo anti-difusivo Lambda_i que mantem a interface
//         fina -- eqs. 24 a 26 de [Talita], eqs. 23 e 24 de [MHL];
//
//      3) a tensao superficial entrando como forca de corpo F = -sigma * kappa * grad(phi),
//         incorporada pelo esquema de forcante de He/Guo -- eqs. 22 e 23 de [Talita].
//
//  Correcoes em relacao ao texto das referencias (erros de digitacao evidentes):
//
//      * [Talita] eq. 23 traz "(c_i - u)/c_s^2 + (c_i - u)/c_s^4"; o segundo termo do esquema de
//        He/Guo e (c_i.u) c_i / c_s^4. A forma correta e exatamente a que a rotina source() desta
//        biblioteca ja calcula, e e ela que se usa aqui.
//
//      * [Talita] eq. 25 traz "(c_i - u)"; como o proprio texto diz que o equilibrio e truncado em
//        primeira ordem, o correto e (c_i . u).
//
//      * [Talita] eq. 22 escreve F_alpha = -sigma*kappa*delta_int, com o lado esquerdo tendo indice
//        livre e o direito sendo escalar: falta a direcao. A forma consistente, que e a que
//        reproduz a lei de Laplace, e F_alpha = -sigma * kappa * d_alpha(phi), ou seja, usando
//        n_alpha |grad phi| = d_alpha(phi) no lugar de delta_int.
//
//  Mapeamento dos parametros de entrada (data_in.txt inalterado):
//
//      fat_R_B  ->  parameters.A_fact    ->  gamma, constante anti-difusiva da eq. 26
//      fat_R_R  ->  parameters.A_fact_R  ->  sigma, coeficiente de tensao superficial
//      tau_R    ->  tempo de relaxacao onde phi = 1
//      tau_B    ->  tempo de relaxacao onde phi = 0
//      wett_R   ->  valor de phi imposto na parede (molhabilidade)
//      recoll   ->  sem uso neste modelo (nao ha recoloracao)
//
//====================================================================================================================//
//====================================================================================================================//




//=============================== Gradiente do campo de fase =========================================================//
//
//      Input: indice do sitio, campo de fase (um valor por sitio), rede e parametros
//      Output: grad_x, grad_y, grad_z
//
//      Usa a estimativa isotropica de rede
//
//              d_alpha(phi) = (1/c_s^2) * soma_i w_i c_i_alpha phi(x + c_i)
//
//      Os vizinhos sao alcancados pelo mapa de propagacao ini_stream, que ja existe: para um elo
//      fluido, ini_stream[pto*nvel + i] = sitio_vizinho*nvel + i, entao a divisao inteira por nvel
//      devolve o indice do vizinho. E leitura apenas, portanto segura em paralelo.
//
//      Em um elo que aponta para solido, ini_solid marca o elo e o valor lido passa a ser wett_R:
//      e assim que a molhabilidade entra no modelo, do mesmo modo que na versao com mediadores.
//
//====================================================================================================================//

#pragma acc routine seq
void grad_phase_MHL ( int pto, const double* ini_phi, LATTICE lattice, PARAMETERS parameters,
                      double& grad_x, double& grad_y, double& grad_z )
{
    grad_x = 0.0;
    grad_y = 0.0;
    grad_z = 0.0;

    const int base = pto * nvel;

    //  i = 0 nao contribui: c_0 = 0.

    for ( int i = 1; i < nvel; i++ )
    {
        const int index = base + i;

        double phi_viz;

        if ( lattice.ini_solid[index] ) phi_viz = parameters.wett_R;

        else phi_viz = ini_phi[ lattice.ini_stream[index] / nvel ];

        const double w_phi = lattice.w[i] * phi_viz;

        grad_x = grad_x + w_phi * lattice.c_i[ i * dim + 0 ];
        grad_y = grad_y + w_phi * lattice.c_i[ i * dim + 1 ];
        grad_z = grad_z + w_phi * ( ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0 );
    }

    const double inv_c_s2 = lattice.one_over_c_s2;

    grad_x = grad_x * inv_c_s2;
    grad_y = grad_y * inv_c_s2;
    grad_z = grad_z * inv_c_s2;
}

//====================================================================================================================//




//=============================== Curvatura da interface =============================================================//
//
//      Input: indice do sitio, gradiente do campo de fase (dim valores por sitio), rede
//      Output: kappa = div(n), com n = grad(phi)/|grad(phi)|
//
//      Mesma estimativa isotropica do gradiente, aplicada componente a componente e somada:
//
//              kappa = (1/c_s^2) * soma_i w_i c_i_alpha n_alpha(x + c_i)
//
//      A normal e formada no vizinho, a partir do gradiente ja armazenado. Onde |grad phi| e
//      desprezivel (dentro das fases) o vizinho nao contribui, o que evita normalizar ruido.
//
//      Em elos solidos ini_stream devolve o proprio sitio, o que equivale a derivada nula de n na
//      parede.
//
//====================================================================================================================//

#pragma acc routine seq
double curvature_MHL ( int pto, const double* ini_grad, LATTICE lattice )
{
    const double eps = 1.0e-12;

    double div = 0.0;

    const int base = pto * nvel;

    for ( int i = 1; i < nvel; i++ )
    {
        const int viz = lattice.ini_stream[ base + i ] / nvel;

        const double gx = ini_grad[ viz * dim + 0 ];
        const double gy = ini_grad[ viz * dim + 1 ];
        const double gz = ( ( dim == 3 ) ? ini_grad[ viz * dim + 2 ] : 0.0 );

        const double mod = sqrt ( gx * gx + gy * gy + gz * gz );

        if ( mod < eps ) continue;

        const double inv_mod = 1.0 / mod;

        div = div + lattice.w[i] * inv_mod * ( lattice.c_i[ i * dim + 0 ] * gx
                                             + lattice.c_i[ i * dim + 1 ] * gy
                                             + ( ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0 ) * gz );
    }

    return div * lattice.one_over_c_s2;
}

//====================================================================================================================//




//=============================== Colisao do modelo Montessori-Hegele-Lauricella =====================================//
//
//      Input: populacao hidrodinamica do sitio, campo de fase e seu gradiente, curvatura,
//             aceleracao externa, rede e parametros
//      Output: populacao pos-colisional; densidade e velocidade do sitio (por referencia)
//
//      Colisao regularizada de segunda ordem, eqs. 17 a 19 de [Talita]:
//
//              f_i^pos = f_i^eq + ( 1 - 1/tau ) f_i^neq + S_i
//
//      com f_i^neq reconstruido do tensor de fluxo de momento fora do equilibrio. A estrutura e a
//      mesma de coll_Reg(); esta reescrita nao chama aquela rotina para nao repetir o calculo de
//      densidade e velocidade, que aqui ja e necessario para o campo de fase.
//
//      A tensao superficial entra como forca de corpo
//
//              F_alpha = - sigma * kappa * d_alpha(phi)
//
//      incorporada pelo forcante de He/Guo, que e o que source() calcula. A velocidade devolvida
//      ja traz a correcao de meia forca, u = ( soma_i c_i f_i + F/2 ) / rho, que e a velocidade
//      fisica e a que deve alimentar o equilibrio do campo de fase.
//
//      O tempo de relaxacao e interpolado pela fase em frequencia,
//
//              omega = phi * omega_R + ( 1 - phi ) * omega_B,
//
//      que atravessa a interface suavemente. Com tau_R = tau_B -- o caso das duas referencias, em
//      que os fluidos tem a mesma viscosidade -- a interpolacao e inocua.
//
//====================================================================================================================//

#pragma acc routine seq
void coll_MHL ( double* f, double phi, double grad_x, double grad_y, double grad_z, double kappa,
                double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters,
                double& rho_out, double& ux_out, double& uy_out, double& uz_out )
{
    const double eps = 1.0e-12;

    //------------------- Tempo de relaxacao interpolado pela fase -----------------------------------//

    double phi_lim = phi;

    if ( phi_lim < 0.0 ) phi_lim = 0.0;
    if ( phi_lim > 1.0 ) phi_lim = 1.0;

    const double omega_R = ( parameters.tau_R > eps ) ? 1.0 / parameters.tau_R : 1.0;
    const double omega_B = ( parameters.tau_B > eps ) ? 1.0 / parameters.tau_B : 1.0;

    const double omega = phi_lim * omega_R + ( 1.0 - phi_lim ) * omega_B;

    const double tau = 1.0 / omega;

    //------------------- Densidade e velocidade -----------------------------------------------------//

    double vx, vy, vz, rho;

    calcula ( f, vx, vy, vz, rho, lattice );

    rho_out = rho;

    if ( rho <= eps )
    {
        ux_out = 0.0;
        uy_out = 0.0;
        uz_out = 0.0;

        return;
    }

    //------------------- Forca de tensao superficial ------------------------------------------------//
    //
    //  sigma = parameters.A_fact_R (campo fat_R_R do data_in.txt).

    const double sigma = parameters.A_fact_R;

    const double Fx = - sigma * kappa * grad_x;
    const double Fy = - sigma * kappa * grad_y;
    const double Fz = - sigma * kappa * grad_z;

    //  Aceleracao total: a externa mais a de origem interfacial.

    const double inv_rho = 1.0 / rho;

    const double a_x = acc_x + Fx * inv_rho;
    const double a_y = acc_y + Fy * inv_rho;
    const double a_z = acc_z + Fz * inv_rho;

    //  Velocidade fisica, com a correcao de meia forca.

    const double ux = vx + 0.5 * a_x;
    const double uy = vy + 0.5 * a_y;
    const double uz = vz + 0.5 * a_z;

    ux_out = ux;
    uy_out = uy;
    uz_out = uz;

    //------------------- Termo forcante (He/Guo) ----------------------------------------------------//

    double S[nvel];

    source ( a_x * rho, a_y * rho, a_z * rho, ux, uy, uz, rho, tau, S, lattice );

    //------------------- Equilibrio de segunda ordem ------------------------------------------------//

    double f_eq[nvel];

    dist_eq ( f_eq, ux, uy, uz, rho, lattice );

    //------------------- Tensor de fluxo de momento fora do equilibrio ------------------------------//

    double Pi_xx = 0.0;
    double Pi_yy = 0.0;
    double Pi_zz = 0.0;

    double Pi_xy = 0.0;
    double Pi_xz = 0.0;
    double Pi_yz = 0.0;

    for ( int i = 0; i < nvel; i++ )
    {
        const double f_neq = f[i] - f_eq[i];

        const double cx = lattice.c_i[ i * dim + 0 ];
        const double cy = lattice.c_i[ i * dim + 1 ];
        const double cz = ( ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0 );

        Pi_xx = Pi_xx + f_neq * cx * cx;
        Pi_yy = Pi_yy + f_neq * cy * cy;
        Pi_zz = Pi_zz + f_neq * cz * cz;

        Pi_xy = Pi_xy + f_neq * cx * cy;
        Pi_xz = Pi_xz + f_neq * cx * cz;
        Pi_yz = Pi_yz + f_neq * cy * cz;
    }

    //------------------- Reconstrucao regularizada e colisao ----------------------------------------//

    const double one_over_2cs4 = 0.5 * lattice.one_over_c_s2 * lattice.one_over_c_s2;

    const double fator_neq = 1.0 - 1.0 / tau;

    for ( int i = 0; i < nvel; i++ )
    {
        const int ind = i * dim * dim;

        const double f_1 = one_over_2cs4 * lattice.w[i]
                         * ( lattice.Q_i[ ind     ] * Pi_xx
                           + lattice.Q_i[ ind + 4 ] * Pi_yy
                           + lattice.Q_i[ ind + 8 ] * Pi_zz
                           + 2.0 * lattice.Q_i[ ind + 1 ] * Pi_xy
                           + 2.0 * lattice.Q_i[ ind + 2 ] * Pi_xz
                           + 2.0 * lattice.Q_i[ ind + 5 ] * Pi_yz );

        f[i] = f_eq[i] + fator_neq * f_1 + S[i];
    }
}

//====================================================================================================================//




//=============================== Evolucao do campo de fase (MHL) ====================================================//
//
//      Input: indice do sitio, campo de fase, velocidade do sitio, gradiente do campo de fase,
//             rede e parametros
//      Output: escreve g nos vizinhos (lattice.inif_m_new)
//
//      Eqs. 24 a 26 de [Talita] / 23 e 24 de [MHL]:
//
//              g_i( x + c_i, t + 1 ) = g_i^eq( phi, u ) + Lambda_i
//
//              g_i^eq   = w_i phi ( 1 + c_i.u / c_s^2 )
//              Lambda_i = w_i gamma phi ( 1 - phi ) ( c_i . n )
//
//      Como a relaxacao e unitaria, o valor pos-colisional e o proprio equilibrio: a rotina escreve
//      direto no vizinho, exatamente como a etapa de emissao de mediadores ja fazia. A difusividade
//      da interface fica fixada pela rede, D = c_s^2 ( 1 - 1/2 ) = c_s^2 / 2 = 1/6 em D3Q19; quem
//      contrabalanca essa difusao e o termo Lambda_i, cuja intensidade e gamma.
//
//      Conservacao: soma_i g_i^eq = phi e soma_i Lambda_i = 0, entao a emissao preserva o campo de
//      fase; a propagacao e uma bijecao e tambem o preserva.
//
//      Em elos que apontam para solido escreve-se w_i * wett_R, que impoe o valor de phi visto a
//      partir da parede -- e a mesma convencao de molhabilidade da versao com mediadores.
//
//====================================================================================================================//

#pragma acc routine seq
void emit_phase_MHL ( int pto, double phi, double ux, double uy, double uz,
                      double grad_x, double grad_y, double grad_z, LATTICE lattice,
                      PARAMETERS parameters )
{
    const double eps = 1.0e-12;

    //  gamma = parameters.A_fact (campo fat_R_B do data_in.txt).

    const double gamma = parameters.A_fact;

    const double inv_c_s2 = lattice.one_over_c_s2;

    //------------------- Normal unitaria a interface ------------------------------------------------//

    double n_x = 0.0;
    double n_y = 0.0;
    double n_z = 0.0;

    const double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    if ( mod_grad > eps )
    {
        const double inv_mod = 1.0 / mod_grad;

        n_x = grad_x * inv_mod;
        n_y = grad_y * inv_mod;
        n_z = grad_z * inv_mod;
    }

    const double anti_dif = gamma * phi * ( 1.0 - phi );

    //------------------------------------------------------------------------------------------------//

    const int base = pto * nvel;

    for ( int i = 0; i < nvel; i++ )
    {
        const int index = base + i;

        const double w_i = lattice.w[i];

        const double cx = lattice.c_i[ i * dim + 0 ];
        const double cy = lattice.c_i[ i * dim + 1 ];
        const double cz = ( ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0 );

        const double g_eq = w_i * phi * ( 1.0 + ( cx * ux + cy * uy + cz * uz ) * inv_c_s2 );

        const double lambda_i = w_i * anti_dif * ( cx * n_x + cy * n_y + cz * n_z );

        lattice.inif_m_new[ lattice.ini_stream[index] ] = lattice.ini_solid[index]
                                                        ? w_i * parameters.wett_R
                                                        : ( g_eq + lambda_i );
    }
}

//====================================================================================================================//



//=============================== Collision step using the Shan-Chen model (immiscible, two components) ==============//
//
//      Input: red and blue distribution functions, force densities acting on each component,
//             lattice and parameters
//      Output: post-collisional distribution functions
//
//      Modelo de Shan & Chen (1993/1994) para dois componentes imisciveis. O que caracteriza o
//      modelo, e o que uma colisao BGK aplicada separadamente a cada especie nao reproduz, e que
//      as duas relaxam para equilibrios construidos sobre uma unica velocidade composta:
//
//              u' = soma_k ( rho_k u_k / tau_k ) / soma_k ( rho_k / tau_k )
//
//      e a forca entra deslocando essa velocidade em cada equilibrio:
//
//              u_k^eq = u' + tau_k F_k / rho_k
//
//      Sem a velocidade comum os dois fluidos so trocam momento pela forca de interacao: nao ha
//      acoplamento viscoso entre eles e um pode deslizar sobre o outro livremente. Com ela, a
//      mistura se comporta como um unico fluido de velocidade baricentrica.
//
//      As forcas entram como DENSIDADE de forca (forca por volume), que e o que force_SC devolve
//      e o que a equacao de estado do modelo pressupoe -- p = rho c_s^2 + (G c_s^2/2) psi^2 vem de
//      F = -G psi grad(psi) por volume. Atencao: coll_BGK() e op_bgk() desta mesma biblioteca
//      recebem ACELERACAO e fazem F = acc * rho internamente; passar uma no lugar da outra
//      multiplica a interacao por rho.
//
//      O momento total se conserva porque a soma das forcas sobre todo o dominio se anula, ainda
//      que localmente F_R + F_B nao seja zero.
//
//====================================================================================================================//

#pragma acc routine seq
void coll_SC ( double* f_R, double* f_B, double Fx_R, double Fy_R, double Fz_R,
               double Fx_B, double Fy_B, double Fz_B, LATTICE lattice, PARAMETERS parameters )
{
    const double eps = 1.0e-30;

    //------------------- Momentos de cada especie ---------------------------------------------------//

    double vx_R, vy_R, vz_R, rho_R;

    double vx_B, vy_B, vz_B, rho_B;

    calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );

    calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

    //------------------- Velocidade composta --------------------------------------------------------//

    double tau_R = parameters.tau_R;
    double tau_B = parameters.tau_B;

    if ( tau_R <= eps ) tau_R = parameters.tau;
    if ( tau_B <= eps ) tau_B = parameters.tau;

    double peso_R = rho_R / tau_R;
    double peso_B = rho_B / tau_B;

    double peso = peso_R + peso_B;

    double ux = 0.0, uy = 0.0, uz = 0.0;

    if ( peso > eps )
    {
        ux = ( peso_R * vx_R + peso_B * vx_B ) / peso;
        uy = ( peso_R * vy_R + peso_B * vy_B ) / peso;
        uz = ( peso_R * vz_R + peso_B * vz_B ) / peso;
    }

    //------------------- Vermelho: equilibrio deslocado pela forca ----------------------------------//

    if ( rho_R > eps )
    {
        double desl = tau_R / rho_R;

        double f_eq[nvel];

        dist_eq ( f_eq, ux + desl * Fx_R, uy + desl * Fy_R, uz + desl * Fz_R, rho_R, lattice );

        double one_over_tau = 1.0 / tau_R;

        for ( int i = 0; i < nvel; i++ ) f_R[i] = f_R[i] + ( f_eq[i] - f_R[i] ) * one_over_tau;
    }

    //------------------- Azul ------------------------------------------------------------------------//

    if ( rho_B > eps )
    {
        double desl = tau_B / rho_B;

        double f_eq[nvel];

        dist_eq ( f_eq, ux + desl * Fx_B, uy + desl * Fy_B, uz + desl * Fz_B, rho_B, lattice );

        double one_over_tau = 1.0 / tau_B;

        for ( int i = 0; i < nvel; i++ ) f_B[i] = f_B[i] + ( f_eq[i] - f_B[i] ) * one_over_tau;
    }
}

//====================================================================================================================//




//=============================== Shan-Chen pseudopotential ==========================================================//
//
//      Input: local density, reference density
//      Output: pseudopotential psi = rho_0 [ 1 - exp( - rho / rho_0 ) ]
//
//      Forma classica de Shan & Chen: cresce como rho quando rho << rho_0 e satura em rho_0, o que
//      limita a forca em densidades altas e da estabilidade ao modelo.
//
//====================================================================================================================//

#pragma acc routine seq
double psi_SC ( double rho, double rho_0 )
{
    if ( rho_0 <= 0.0 ) return rho;

    return rho_0 * ( 1.0 - exp ( - rho / rho_0 ) );
}

//====================================================================================================================//

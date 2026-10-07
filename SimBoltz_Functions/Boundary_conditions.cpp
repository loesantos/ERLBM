


//====================================================================================================================//
//
//   A GEOMETRIA NO DISPOSITIVO
//
//   As fronteiras deste arquivo sao as unicas rotinas da biblioteca que saem do indice de fluido
//   e voltam ao indice de voxel: elas leem  geometry.ini  para achar, em cada coluna, o sitio e o
//   seu vizinho.  Os kernels do laco principal nunca precisam disso, entao um programa tipico
//   NAO mapeia  geometry.ini  -- e ai uma clausula  present  sobre ele falha em execucao, com
//
//       FATAL ERROR: data in PRESENT clause was not found on device
//
//   A struct GEOMETRY inteira nao pode ir para o dispositivo, porque contem um std::string; so o
//   ponteiro do mapa de indices pode.  Como aqui  geometry.ini  e apenas LIDO, a solucao segura e
//   declara-lo com  copyin :  do OpenACC 2.5 em diante  copyin  e "presente-ou-copia", isto e,
//   reaproveita a copia do dispositivo quando ela existe e faz uma temporaria quando nao existe.
//
//   As populacoes continuam em  present , e isso e proposital: elas sao ESCRITAS, e uma copia
//   temporaria receberia as escritas e as jogaria fora ao fim da regiao.
//
//   Copiar o mapa a cada chamada funciona mas e desperdicio.  O aviso abaixo sai uma unica vez e
//   diz como evitar: mapear a geometria uma vez, no programa.
//
//====================================================================================================================//

#ifdef _OPENACC

static void confere_geometria_no_dispositivo ( const int *geometry_ini, int geometry_size )
{
	static bool avisado = false;

	if ( avisado ) return;

	if ( acc_is_present ( ( void * ) geometry_ini,
	                      ( size_t ) geometry_size * sizeof ( int ) ) ) return;

	avisado = true;

	cerr << "\n   AVISO: geometry.ini nao esta no dispositivo.  As condicoes de fronteira vao"
	        "\n          copia-lo a cada chamada; o resultado esta certo, mas custa uma"
	        "\n          transferencia por chamada.  Para evitar, acrescente ao programa, junto"
	        "\n          dos outros  enter data :"
	        "\n"
	        "\n              int *geo_ini = geometry.ini;"
	        "\n              #pragma acc enter data copyin( geo_ini[ 0 : ( long long ) nx * ny * nz ] )"
	        "\n"
	        "\n          e o  exit data delete  correspondente ao fim.\n" << endl;
}

#define CONFERE_GEOMETRIA( p, n )   confere_geometria_no_dispositivo ( ( p ), ( n ) )

#else

#define CONFERE_GEOMETRIA( p, n )   ( ( void ) 0 )

#endif


//=============================== Boundary condition - null derivative of the velocity in x direction ================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void devnull_rho_x ( int pos_x, double rho, const GEOMETRY &geometry, LATTICE &lattice )
{
    int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;

    // Somente dados triviais e ponteiros registrados entram no kernel. GEOMETRY
    // contem std::string e nao pode ser capturada pelo device; LATTICE tambem nao
    // deve ser desreferenciada indiretamente a partir de uma referencia do host.
    int *geometry_ini = geometry.ini;
    double *populations = lattice.inif;
    int geometry_size = nx * ny * nz;
    int distribution_size = geometry.fluid * nvel;

    //  Copia LOCAL da LATTICE, por valor.  LATTICE e POD e os pesos sao o array fixo w[nvel]
    //  DENTRO da struct, entao eles viajam junto na copia -- nao ha ponteiro a anexar e nao e
    //  preciso mapear nada.  Copiar aqui, no hospedeiro, evita desreferenciar a referencia do
    //  host dentro do kernel, que e o que o comentario acima adverte.  Os ponteiros que vierem
    //  na copia sao enderecos de hospedeiro e NAO podem ser usados no kernel; dist_eq() so le
    //  lattice.w, entao o uso abaixo e legitimo.
    LATTICE lat = lattice;
    
    int infx = -1;
    
    if ( pos_x < nx / 2 ) infx = 1;

	CONFERE_GEOMETRIA ( geometry_ini, geometry_size );

#ifdef _OPENACC
	//  geometry_ini so e lido: copyin serve, e reaproveita a copia do programa quando existe.
	//  populations e escrito: tem de estar mesmo mapeado, entao continua em present.
	#pragma acc parallel loop collapse(2) firstprivate(lat) \
		copyin(geometry_ini[0:geometry_size]) present(populations[0:distribution_size])
#else
	#pragma omp parallel for collapse(2)
#endif

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int* meio_pto = geometry_ini + pos_x + y * nx + z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = populations + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = geometry_ini + ( pos_x + infx ) + y * nx + z * ny * nx;

                if ( *meio_adj )
                {
                    double* f_adj = populations + ( *meio_adj - 1 ) * nvel;
                    
                    double rho_adj = density( f_adj );
                    
                    double fator = rho / rho_adj;
                
					for ( int i = 0; i < nvel; i++ ) f_pto[i] = f_adj[i] * fator;
                }
                else
                {
                    //  Equilibrio em repouso:  f_eq_i = w_i rho .  Antes estava escrito a mao com
                    //  os pesos de D3Q19 ( 1/3, 1/18, 1/36 ), o que dava valores ERRADOS, sem
                    //  aviso, em qualquer outra rede -- em D2Q9 os pesos sao 4/9, 1/9 e 1/36.
                    //  dist_eq() e  acc routine seq  e faz exatamente  feq[i] = lattice.w[i] * T ,
                    //  com os pesos que a montagem da rede gravou.  Vale para D2Q9, D3Q15, D3Q19
                    //  e D3Q27 sem nenhum caso especial, e sem depender da ordem das velocidades:
                    //  w e indexado como f.
                    dist_eq ( f_pto, rho, lat );
                }
            }
        }
    }
}
//====================================================================================================================//




//=============================== Boundary condition - null derivative of the velocity in z direction ================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void devnull_rho_z ( int pos_z, double rho, GEOMETRY geometry, LATTICE lattice )
{
    int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;
    
    int infz = -1;
    
    if ( pos_z < nz / 2 ) infz = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int x = 0; x < nx; x++ )
        {
            int* meio_pto = geometry.ini + x + y * nx + pos_z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = lattice.inif + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = geometry.ini + x + y * nx + ( pos_z + infz ) * ny * nx;

                if ( *meio_adj )
                {
                    double* f_adj = lattice.inif + ( *meio_adj - 1 ) * nvel;
                    
                    double rho_adj = density( f_adj );
                    
                    double fator = rho / rho_adj;
                
					for ( int i = 0; i < nvel; i++ ) f_pto[i] = f_adj[i] * fator;
                }
                else
                {
                    dist_eq (f_pto, 0., 0., 0., rho, lattice );
                }
            }
        }
    }
}
//====================================================================================================================//



//=============================== Boundary condition - null derivative of the velocity in z direction ================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void imp_eq_z ( int pos_z, double rho, double vx, double vy, double vz, double* ini_f, 
				GEOMETRY geometry, LATTICE lattice )
{
    int nx = geometry.nx;
    int ny = geometry.ny;
    
    int z = pos_z;
    
	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int x = 0; x < nx; x++ )
        {
            int* meio_pto = geometry.ini + x + y * nx + z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = ini_f + ( *meio_pto - 1 ) * nvel;
				
                dist_eq (f_pto, vx, vy, vz, rho, lattice );
               
            }
        }
    }
}
//====================================================================================================================//



//=============================== Boundary condition - null derivative of the velocity in z direction ================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void devnull_rho_Red_z ( int pos_z, double rho, GEOMETRY geometry, LATTICE lattice )
{
    int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;
    
    int infz = -1;
    
    if ( pos_z < nz / 2 ) infz = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int x = 0; x < nx; x++ )
        {
            int* meio_pto = geometry.ini + x + y * nx + pos_z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = lattice.inif_R + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = geometry.ini + x + y * nx + ( pos_z + infz ) * ny * nx;

                if ( *meio_adj )
                {
                    double* f_adj = lattice.inif_R + ( *meio_adj - 1 ) * nvel;
                    
                    double rho_adj = density( f_adj );
                    
                    double fator = rho / rho_adj;
                
					for ( int i = 0; i < nvel; i++ ) f_pto[i] = f_adj[i] * fator;
                }
                else
                {
                    dist_eq ( f_pto, 0., 0., 0., rho, lattice );
                }
            }
        }
    }
}
//====================================================================================================================//



//=============================== Boundary condition - null derivative of the velocity in z direction ================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void devnull_rho_Blue_z ( int pos_z, double rho, GEOMETRY geometry, LATTICE lattice )
{
    int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;
    
    int infz = -1;
    
    if ( pos_z < nz / 2 ) infz = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int x = 0; x < nx; x++ )
        {
            int* meio_pto = geometry.ini + x + y * nx + pos_z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = lattice.inif_B + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = geometry.ini + x + y * nx + ( pos_z + infz ) * ny * nx;

                if ( *meio_adj )
                {
                    double* f_adj = lattice.inif_B + ( *meio_adj - 1 ) * nvel;
                    
                    double rho_adj = density( f_adj );
                    
                    double fator = rho / rho_adj;
                
					for ( int i = 0; i < nvel; i++ ) f_pto[i] = f_adj[i] * fator;
                }
                else
                {
                    dist_eq (f_pto, 0., 0., 0., rho, lattice );
                }
            }
        }
    }
}
//====================================================================================================================//



//============================ Zou & He Boundary - vx imposed - x direction ==========================================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void ZouHeD3Q19_vx_x ( int pos_x, double vx, GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;
    
    int x = pos_x;

    int infx = 0;

    if ( pos_x < nx / 2 ) infx = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = geometry.ini + x + y * nx + z * ny * nx;

            if ( *meio && infx )
            {
                double *f = lattice.inif + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xminus = f[1] + f[7] + f[9] + f[11] + f[13];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
                                     + f[15] + f[16] + f[17] + f[18];

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double rho = ( 2.0 * sum_f_xminus + sum_f_xzero ) / ( 1.0 - vx );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[1]  = f[1] + rho * vx / 3.0;

                f[7]  = f[7] + rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[9]  = f[9] + rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[11] = f[11] + rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[13] = f[13] + rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
            else if ( *meio )
            {
                double *f = lattice.inif + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xplus = f[2] + f[8] + f[10] + f[12] + f[14];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
                                     + f[15] + f[16] + f[17] + f[18];

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double rho = ( 2.0 * sum_f_xplus + sum_f_xzero ) / ( 1.0 - vx );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[2]  = f[2] - rho * vx / 3.0;

                f[10] = f[10] - rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[8]  = f[8] - rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[14] = f[14] - rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[12] = f[12] - rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
        }
    }
}

//====================================================================================================================//




//===================== Boundary (Zou & He) imposing a parabollic velocity field in the inlet ========================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void parabolic_D3Q19_x ( int pos_x, double vx_max, int L, GEOMETRY geometry, LATTICE lattice  )
                       
{	
	int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;
    
    int x = pos_x;
    
    double *vel = new double[ ny ];

    for ( int h = 0; h < ( L + 1 ); h++ )
    {
        double d_h = ( double ) h - 0.5;

        vel[ h ] = -4.0 * vx_max / ( L * L ) * ( d_h ) * ( d_h - L );       
    }

    for ( int z = 0; z < nz; z++ )
    {
        int h = 0;

        for ( int y = 0; y < ny; y++ )
        {
            int *meio = geometry.ini + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //--------------------------------------------------------------------------------//

                h++;
                                
                double vx = vel[ h ];

                //--------------------------------------------------------------------------------//
                
                double *f = lattice.inif + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xminus = f[2] + f[8] + f[10] + f[12] + f[14];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
										+ f[15] + f[16] + f[17] + f[18];

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double rho = ( 2.0 * sum_f_xminus + sum_f_xzero ) / ( 1.0 - vx );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[1]  = f[2] + rho * vx / 3.0;

                f[7]  = f[8] + rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[9]  = f[10] + rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[11] = f[12] + rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[13] = f[14] + rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
        }
    }
		
    delete[] vel;
}
//====================================================================================================================//




//===================== Open Boundary- direção y =========-----=======================================================//
//
//      Input: position, geometry, lattice
//      Output:
//
//====================================================================================================================//

void OpenBoundaryD2Q9_y ( int pos_y, GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
    int ny = geometry.ny;

	int delta_y = -1;

    if ( pos_y < ny / 2 ) delta_y = 1;

    for ( int x = 0; x < nx; x++ )
    {   
		int* meio = geometry.ini + x + pos_y * nx;
		int* meio_adj = geometry.ini + x + ( pos_y + delta_y ) * nx;

		double* f = lattice.inif + ( *meio - 1 ) * nvel;
		double* f_adj = lattice.inif + ( *meio_adj - 1 ) * nvel;
		
		if ( delta_y > 0 ) 
		{
			f[3] = f_adj[3];
			f[5] = f_adj[5];
			f[8] = f_adj[8];
 		}
 		else
 		{
			f[4] = f_adj[4];
			f[6] = f_adj[6];
			f[7] = f_adj[7];
 		}
	}
}
//====================================================================================================================//



//===================== Open Boundary- direção y =========-----=======================================================//
//
//      Input: position, geometry, lattice
//      Output:
//
//====================================================================================================================//

void OpenBoundaryD2Q9_x ( int pos_x, GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
    int ny = geometry.ny;

	int delta_x = -1;

    if ( pos_x < nx / 2 ) delta_x = 1;

    for ( int y = 0; y < ny; y++ )
    {   
		int* meio = geometry.ini + pos_x + y * nx;
		int* meio_adj = geometry.ini + ( pos_x + delta_x ) + y * nx;

		double* f = lattice.inif + ( *meio - 1 ) * nvel;
		double* f_adj = lattice.inif + ( *meio_adj - 1 ) * nvel;
		
		if ( delta_x > 0 ) 
		{
			f[1] = f_adj[1];
			f[5] = f_adj[5];
			f[7] = f_adj[7];
 		}
 		else
 		{
			f[2] = f_adj[2];
			f[6] = f_adj[6];
			f[8] = f_adj[8];
 		}
	}
}
//====================================================================================================================//




//===================== Condição de contorno de derivada nula da velocidade - direção x ==============================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//=====================================================================================================================//

void convective_x ( int posx, GEOMETRY geometry, LATTICE lattice )
{
	//double Mach = 0.1;
		
	//double U_max = Mach * sqrt ( lattice.c_s2 );
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

	#pragma omp parallel for

	for ( int y = 0; y < ny; y++ )
	{
		for ( int z = 0; z < nz; z++ )
		{
			int *meio = geometry.ini + x + y * nx + z * ny * nx;

			if ( *meio )
			{
				//int alea = rand() % 100;
				
				//----------- Calcula a velocidade -----------------------------------------------//
				
				int* meio_A = geometry.ini + ( x + infx ) + y * nx + z * ny * nx;

				double *f_A = lattice.inif + ( *meio_A - 1 ) * nvel;
				
				double U, vy, vz, rho;
								
				calcula ( f_A, U, vy, vz, rho, lattice );		
				
				//if ( alea == 0 )	cout << "\nU = " << U << endl;
												
				//--------------------------------------------------------------------------------//
								
				meio = geometry.ini + x + y * nx + z * ny * nx;

				double *f = lattice.inif + ( *meio - 1 ) * nvel;
								
				//double U_old, vy_old, vz_old, rho_old;
				/*
				if ( alea == 0 )
				{
					for ( int i = 0; i < nvel; i ++ ) cout << "\nf[" << i << "] = " << f[i];
					
					calcula ( f, U_old, vy_old, vz_old, rho_old, lattice );		
				
					cout << "\nU_old = " << U_old << "		rho_old = " << rho_old << endl;
				
					getchar();
				}
				*/
				double factor = 1.0 / ( 1.0 + U );
				
				for ( int i = 0 ; i < nvel; i++ )
				{
					f[i] = factor * ( f[i] + U * f_A[i] );
				}

				//--------------------------------------------------------------------------------//
			}
		}
	}
}
//====================================================================================================================//




//============================ Zou & He Boundary - vx imposed - x direction ==========================================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void ZouHeD2Q9_vx_x ( int pos_x, double vx, GEOMETRY geometry, LATTICE lattice )
{
    int x = pos_x;
    
    int nx = geometry.nx;
    int ny = geometry.ny;
    
    for ( int y = 0; y < ny; y++ )
    {
        int *meio = geometry.ini + x + y * nx;

        if ( *meio )
        {
            double *f = lattice.inif + ( *meio - 1 ) * nvel;

			double rho = ( f[0] + f[3] + f[4] + 2 * ( f[2] + f[8] + f[6] ) ) / ( 1. - vx );
			
			f[1] = ( rho / 9.0 ) * ( 1.0 + 3.0 * vx + 3.0 * vx * vx );
			
			f[5] = ( rho / 36.0 ) * ( 1.0 + 3.0 * vx + 3.0 * vx * vx );
			
			f[7] = ( rho / 36.0 ) * ( 1.0 + 3.0 * vx + 3.0 * vx * vx );
		
        }
    }
}
//====================================================================================================================//




//===================== Boundary (Zou & He) imposing a parabollic velocity field in the inlet ========================//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void parabolic_D2Q9_x ( int pos_x, double vx_max, int L, GEOMETRY geometry, LATTICE lattice  )
                       
{	
	int nx = geometry.nx;
    int ny = geometry.ny;
    
    int x = pos_x;
    
    double *vel = new double[ ny ];

    for ( int h = 0; h < ( L + 1 ); h++ )
    {
        double d_h = ( double ) h - 0.5;

        vel[ h ] = -4.0 * vx_max / ( L * L ) * ( d_h ) * ( d_h - L );       
    }

	//----------------------------------------------------------------------------------------//
	
	int h = 0;

	for ( int y = 0; y < ny; y++ )
	{
		int *meio = geometry.ini + x + y * nx;

		if ( *meio )
		{
			//--------------------------------------------------------------------------------//

			h++;
							
			double vx = vel[ h ];

			//--------------------------------------------------------------------------------//
			
			double *f = lattice.inif + ( *meio - 1 ) * nvel;

			double rho = ( f[0] + f[3] + f[4] + 2 * ( f[2] + f[8] + f[6] ) ) / ( 1. - vx );
			
			f[1] = ( rho / 9.0 ) * ( 1.0 + 3.0 * vx + 3.0 * vx * vx );
			
			f[5] = ( rho / 36.0 ) * ( 1.0 + 3.0 * vx + 3.0 * vx * vx );
			
			f[7] = ( rho / 36.0 ) * ( 1.0 + 3.0 * vx + 3.0 * vx * vx );
		}
	}

    delete[] vel;
}


//=============================== Boundary condition - imposing temperature for heat conduction - x direction ========//
//
//      Input: position, density, geometry, lattice
//      Output: distribution function
//
//====================================================================================================================//

void imp_T_x ( int pos_x, double T_imp, GEOMETRY geometry, LATTICE lattice )
{
    int nx = geometry.nx;
    int ny = geometry.ny;
    int nz = geometry.nz;
    
    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int* meio_pto = geometry.ini + pos_x + y * nx + z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = lattice.inif + ( *meio_pto - 1 ) * nvel;
				
				dist_eq ( f_pto, T_imp, lattice );                
            }
        }
    }
}
//====================================================================================================================//


//=============================== Boundary condition - null derivative on an x face ==================================//
//
//      Derivada nula em x para as DUAS cores e para os mediadores: a coluna de fronteira recebe,
//      populacao por populacao, uma copia da coluna interior vizinha.
//
//      Serve ao ensaio de onda capilar, em que as duas faces x ficam no seio de um fluido puro e em
//      repouso.  'rho_R_alvo' e 'rho_B_alvo' reescalam a copia para uma densidade dada, o que ancora
//      a pressao no campo distante; NEGATIVO desliga a reescala daquela cor e deixa a copia crua.
//      O uso normal ancora so a cor presente em cada face:
//
//          devnull_RB_x ( 0,      rho_ini_R, -1.0, geometry, lattice );
//          devnull_RB_x ( nx - 1, -1.0, rho_ini_B, geometry, lattice );
//
//      Sem ancoragem a massa total escorre devagar -- medido em 0.3 % a cada 3000 passos no ensaio
//      de calibracao.  Com ela, a deriva some.
//
//      Chame no ALTO do passo de tempo, antes da emissao dos mediadores.
//
//      Input : pos_x, rho_R_alvo, rho_B_alvo, geometry, lattice
//      Output: inif_R, inif_B e inif_m na coluna pos_x
//
//====================================================================================================================//

void devnull_RB_x ( int pos_x, double rho_R_alvo, double rho_B_alvo,
                    GEOMETRY geometry, LATTICE lattice )
{
	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	//  So dados triviais e ponteiros registrados entram no kernel ( ver devnull_rho_x ).

	int *geometry_ini = geometry.ini;

	double *pop_R = lattice.inif_R;
	double *pop_B = lattice.inif_B;
	double *pop_m = lattice.inif_m;

	const int geometry_size = nx * ny * nz;
	const int distribution_size = geometry.fluid * nvel;

	const bool tem_med = ( pop_m != nullptr );

	const int infx = ( pos_x < nx / 2 ) ? 1 : -1;

	CONFERE_GEOMETRIA ( geometry_ini, geometry_size );

#ifdef _OPENACC
	#pragma acc parallel loop collapse(2) \
		copyin(geometry_ini[0:geometry_size]) \
		present(pop_R[0:distribution_size], pop_B[0:distribution_size])
#else
	#pragma omp parallel for collapse(2)
#endif

	for ( int y = 0; y < ny; y++ )
	{
		for ( int z = 0; z < nz; z++ )
		{
			const int pos = pos_x + y * nx + z * nx * ny;

			const int meio_pto = geometry_ini[pos];

			if ( meio_pto == 0 ) continue;

			const int meio_adj = geometry_ini[ ( pos_x + infx ) + y * nx + z * nx * ny ];

			if ( meio_adj == 0 ) continue;

			const int pto = meio_pto - 1;
			const int adj = meio_adj - 1;

			//  Reescala so quando a cor esta presente: nas faces a cor minoritaria tem
			//  densidade ~ 0 e dividir por ela nao faria sentido.

			double fat_R = 1.0, fat_B = 1.0;

			if ( rho_R_alvo >= 0.0 )
			{
				const double rho_adj = density ( pop_R + adj * nvel );

				if ( rho_adj > 1.e-12 ) fat_R = rho_R_alvo / rho_adj;
			}

			if ( rho_B_alvo >= 0.0 )
			{
				const double rho_adj = density ( pop_B + adj * nvel );

				if ( rho_adj > 1.e-12 ) fat_B = rho_B_alvo / rho_adj;
			}

			for ( int i = 0; i < nvel; i++ )
			{
				pop_R[ pto * nvel + i ] = pop_R[ adj * nvel + i ] * fat_R;
				pop_B[ pto * nvel + i ] = pop_B[ adj * nvel + i ] * fat_B;
			}

		}
	}

	//  Os mediadores vao num laco proprio.  Eles NAO podem entrar na clausula  present  do laco
	//  acima porque  pop_m  e nullptr nos modelos sem mediador, e  present  sobre ponteiro nulo
	//  e erro.  Deixar  pop_m  fora da clausula tambem nao serve: o nvc++ nao descobre o tamanho
	//  da copia no dispositivo e aborta a regiao
	//  ( "Could not find allocated-variable index for symbol - pop_m" ).
	//  Com o teste  tem_med  no hospedeiro, a clausula so e avaliada quando o ponteiro existe.

	if ( tem_med )
	{
#ifdef _OPENACC
		#pragma acc parallel loop collapse(2) \
			copyin(geometry_ini[0:geometry_size]) present(pop_m[0:distribution_size])
#else
		#pragma omp parallel for collapse(2)
#endif

		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				const int pos = pos_x + y * nx + z * nx * ny;

				const int meio_pto = geometry_ini[pos];

				if ( meio_pto == 0 ) continue;

				const int meio_adj = geometry_ini[ ( pos_x + infx ) + y * nx + z * nx * ny ];

				if ( meio_adj == 0 ) continue;

				const int pto = meio_pto - 1;
				const int adj = meio_adj - 1;

				for ( int i = 0; i < nvel; i++ )
				{
					pop_m[ pto * nvel + i ] = pop_m[ adj * nvel + i ];
				}
			}
		}
	}
}

//====================================================================================================================//

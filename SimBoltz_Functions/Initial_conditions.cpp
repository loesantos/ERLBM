

//=============================== Initial Condition of a Sphere (immiscible fluids) ==================================//
//
//      Input: radius, position, geometry, lattice, parameters
//      Output: distribution function
//
//====================================================================================================================//

void initial_conditions_sphere ( double radius, int x0, int y0, int z0, double vx, double vy, double vz, GEOMETRY geometry, LATTICE lattice, 
							PARAMETERS parameters )
{
	double R2 = radius * radius;
	
	int *geo = geometry.ini;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{	
				int pos = x + y * nx + z * nx * ny;			
				
				//  CORRECAO: nos sitios solidos geo[pos] vale zero, e o indice geo[pos] - 1 fica
				//  negativo -- a escrita caia antes do inicio dos vetores. Passava despercebido em
				//  dominios sem solidos, mas corrompia a memoria em qualquer geometria porosa.
				//
				//  No solido nao ha populacao, mas o campo escalar existe em todos os voxels e e
				//  lido pelos vizinhos em gradient(): recebe a molhabilidade, como faz
				//  initial_conditions_random.
				
				if ( geo[pos] == 0 )
				{
					if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = parameters.wett_R;
					
					continue;
				}
				
				double rho_imp_R = 0.0;		
				double rho_imp_B = 0.0;				
					
				double in_sphere = ( x - x0 ) * ( x - x0 ) + ( y - y0 ) * ( y - y0 ) + ( z - z0 ) * ( z - z0 );
				
				if ( in_sphere < R2 ) rho_imp_R = parameters.rho_ini_R;
				
				else rho_imp_B = parameters.rho_ini_B;	
				
				double *f_R = lattice.inif_R + ( geo[pos] - 1 ) * nvel;
				double *f_B = lattice.inif_B + ( geo[pos] - 1 ) * nvel;	
				
				dist_eq ( f_R, vx, vy, vz, rho_imp_R, lattice );
				dist_eq ( f_B, vx, vy, vz, rho_imp_B, lattice );

				//  CORRECAO: os campos de fase nao eram semeados aqui, ao contrario do que fazem
				//  initial_conditions_random e initial_conditions_flat_x. Um programa que use os
				//  mediadores ou o campo escalar partindo de uma esfera lia memoria nao
				//  inicializada. As convencoes seguidas sao as do resto do arquivo:
				//  inif_m guarda rho_R/(rho_R+rho_B), em [0,1], e ini_psi guarda
				//  (rho_R-rho_B)/(rho_R+rho_B), em [-1,1].
				
				double rho_tot = rho_imp_R + rho_imp_B;
				
				if ( lattice.inif_m != nullptr )
				{
					double phase = ( rho_tot > 0.0 ) ? rho_imp_R / rho_tot : 0.0;
					
					double *f_m = lattice.inif_m + ( geo[pos] - 1 ) * nvel;
					
					emite_mediadores ( f_m, phase, lattice );
				}
				
				if ( lattice.ini_psi != nullptr )
				{
					double phase = ( rho_tot > 0.0 ) ? ( rho_imp_R - rho_imp_B ) / rho_tot : 0.0;
					
					lattice.ini_psi[pos] = phase;
				}
			}
		}
	}
}

//====================================================================================================================//



//=============================== Initial Condition of a Sphere (immiscible fluids) ==================================//
//
//      Input: radius, position, geometry, lattice, parameters
//      Output: distribution function
//
//====================================================================================================================//

void initial_conditions_sphereSC ( double radius, int x0, int y0, int z0, GEOMETRY geometry, LATTICE lattice, 
							PARAMETERS parameters )
{
	double R2 = radius * radius;
	
	int *geo = geometry.ini;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{	
				int pos = x + y * nx + z * nx * ny;			
				
				//  CORRECAO: nos sitios solidos geo[pos] vale zero, e o indice geo[pos] - 1 fica
				//  negativo -- a escrita caia antes do inicio dos vetores. Passava despercebido em
				//  dominios sem solidos, mas corrompia a memoria em qualquer geometria porosa.
				//
				//  No solido nao ha populacao, mas o campo escalar existe em todos os voxels e e
				//  lido pelos vizinhos em gradient(): recebe a molhabilidade, como faz
				//  initial_conditions_random.
				
				if ( geo[pos] == 0 )
				{
					if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = parameters.wett_R;
					
					continue;
				}
				
				double rho_imp_R = 0.15 * parameters.rho_ini_R;
				double rho_imp_B = 0.15 * parameters.rho_ini_B;			
					
				double in_sphere = ( x - x0 ) * ( x - x0 ) + ( y - y0 ) * ( y - y0 ) + ( z - z0 ) * ( z - z0 );
				
				if ( in_sphere < R2 ) rho_imp_R = 1.01 * parameters.rho_ini_R;
				
				else rho_imp_B = 0.9 * parameters.rho_ini_B;	
				
				double *f_R = lattice.inif_R + ( geo[pos] - 1 ) * nvel;
				double *f_B = lattice.inif_B + ( geo[pos] - 1 ) * nvel;	
				
				dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
				dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );

				//  CORRECAO: idem a initial_conditions_sphere -- os campos de fase passam a ser
				//  semeados, com as convencoes do resto do arquivo.
				
				double rho_tot = rho_imp_R + rho_imp_B;
				
				if ( lattice.inif_m != nullptr )
				{
					double phase = ( rho_tot > 0.0 ) ? rho_imp_R / rho_tot : 0.0;
					
					double *f_m = lattice.inif_m + ( geo[pos] - 1 ) * nvel;
					
					emite_mediadores ( f_m, phase, lattice );
				}
				
				if ( lattice.ini_psi != nullptr )
				{
					double phase = ( rho_tot > 0.0 ) ? ( rho_imp_R - rho_imp_B ) / rho_tot : 0.0;
					
					lattice.ini_psi[pos] = phase;
				}
			}
		}
	}
}

//====================================================================================================================//



//=============================== Initial Condition from a file (immiscible fluids) ==================================//
//
//      Input: geometry, lattice, parameters
//      Output: distribution function
//
//====================================================================================================================//

void initial_conditions_file ( GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	string nome_geo = geometry.file;
	
	ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();
	
	int *geo = geometry.ini;
	
	int meio = 0;
	
	for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;
                								
				double rho_imp_R = 0.0;		
				double rho_imp_B = 0.0;		
				
				fmatriz >> meio;
								
				if ( meio == 1 ) rho_imp_B = parameters.rho_ini_B;
				if ( meio == 2 ) rho_imp_R = parameters.rho_ini_R;
								
				if ( geo[pos] )
				{
					double *f_R = lattice.inif_R + ( geo[pos] - 1 ) * nvel;
					double *f_B = lattice.inif_B + ( geo[pos] - 1 ) * nvel;	
					
					dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
					dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );
					
					double rho_R = density ( f_R );
					double rho_B = density ( f_B );
            
					double phase = rho_R / ( rho_R + rho_B ); 
					
					if ( lattice.inif_m != nullptr ) 
					{		
						double *f_m = lattice.inif_m + ( geo[pos] - 1  ) * nvel;

						emite_mediadores ( f_m, phase, lattice );
					}
				}
				else if ( lattice.ini_psi != nullptr ) 
				{
					double *psi = lattice.ini_psi;
					
					psi[pos] = parameters.wett_R;
				}	
			}
		}
	}
	fmatriz.close();
}

//====================================================================================================================//



//=============================== Initial Condition from a file (immiscible fluids) ==================================//
//
//      Input: geometry, lattice, parameters
//      Output: distribution function
//
//====================================================================================================================//

void initial_conditions_file ( string name_R, string name_B, GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	ifstream file_rhoR ( name_R );
	ifstream file_rhoB ( name_B );
    
	string line, dump;

	for ( int i = 0; i < 10; i++ ) getline( file_rhoR, dump );
	for ( int i = 0; i < 10; i++ ) getline( file_rhoB, dump );
	
	int *geo = geometry.ini;
	
	for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;
                											
				if ( geo[pos] )
				{
					double rho_imp_R, rho_imp_B;
					
					file_rhoR >> rho_imp_R;			
					file_rhoB >> rho_imp_B;		
					
					double *f_R = lattice.inif_R + ( geo[pos] - 1 ) * nvel;
					double *f_B = lattice.inif_B + ( geo[pos] - 1 ) * nvel;		
					
					dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
					dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );
            
					double phase = rho_imp_R / ( rho_imp_R + rho_imp_B ); 
					
					if ( lattice.inif_m != nullptr ) 
					{		
						double *f_m = lattice.inif_m + ( geo[pos] - 1  ) * nvel;

						emite_mediadores ( f_m, phase, lattice );
					}
				}
				else if ( lattice.ini_psi != nullptr ) 
				{
					double *psi = lattice.ini_psi;
					
					psi[pos] = parameters.wett_R;
				}	
			}
		}
	}
	file_rhoR.close();
	file_rhoB.close();
}

//====================================================================================================================//



//=============================== Randomic Intial condition (immiscible fluids) ======================================//
//
//      Input: concentration of one of the fluids, geometry, lattice, parameters
//      Output: distribution function
//
//====================================================================================================================//

void initial_conditions_random ( double conc_R, GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	default_random_engine generator;
	
	uniform_real_distribution<double> distribution( 0.0, 1.0 ); // Limites da distribuição
	
	if ( conc_R ) 
	{
		if ( lattice.inif_m != nullptr )
		{
			for ( int pto = 0; pto < geometry.fluid; pto++ )
			{		
				double *f_R = lattice.inif_R + ( pto ) * nvel;
				double *f_B = lattice.inif_B + ( pto ) * nvel;
						
				double rho_imp_R = 0.1;		
				double rho_imp_B = 0.1;
				
				if ( distribution( generator ) < conc_R ) rho_imp_R = parameters.rho_ini_R * 0.9;
				
				else rho_imp_B = parameters.rho_ini_B * 0.9;
				
				dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
				dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );
				
				double rho_R = density ( f_R );
				double rho_B = density ( f_B );
					
				double phase = rho_R / ( rho_R + rho_B ); 
				
				double *f_m = lattice.inif_m + ( pto ) * nvel;

				emite_mediadores ( f_m, phase, lattice );
			}
		}
		
		else if ( lattice.ini_psi != nullptr )
		{
			int nx = geometry.nx;
			int ny = geometry.ny;
			int nz = geometry.nz;
		
			int* geo = geometry.ini;
		
			for ( int x = 0; x < nx; x++ )
			{
				for ( int y = 0; y < ny; y++ )
				{
					for ( int z = 0; z < nz; z++ )
					{
						int pos = x + y * nx + z * nx * ny;
						
						int pto = geo[ pos ] - 1;
						
						if ( pto >= 0 )
						{		
							double *f_R = lattice.inif_R + ( pto ) * nvel;
							double *f_B = lattice.inif_B + ( pto ) * nvel;
							
							double rho_imp_R = 0.0;		
							double rho_imp_B = 0.0;
					
							if ( distribution( generator ) < conc_R ) rho_imp_R = parameters.rho_ini_R;
					
							else rho_imp_B = parameters.rho_ini_B;
					
							dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
							dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );
						}
						else
						{
							double *psi = lattice.ini_psi;
					
							psi[pos] = parameters.wett_R;							
						}
					}
				}
			}
		}
		
		else
		{
			int nx = geometry.nx;
			int ny = geometry.ny;
			int nz = geometry.nz;
		
			int* geo = geometry.ini;
		
			for ( int x = 0; x < nx; x++ )
			{
				for ( int y = 0; y < ny; y++ )
				{
					for ( int z = 0; z < nz; z++ )
					{
						int pos = x + y * nx + z * nx * ny;
						
						int pto = geo[ pos ] - 1;
						
						if ( pto >= 0 )
						{		
							double *f_R = lattice.inif_R + ( pto ) * nvel;
							double *f_B = lattice.inif_B + ( pto ) * nvel;
							
							double rho_imp_R = parameters.rho_ini_R;	
							double rho_imp_B = parameters.rho_ini_B;
					
							if ( distribution( generator ) < conc_R )
							{
								rho_imp_R = 0.9 * parameters.rho_ini_R;
								rho_imp_B = 0.1 * parameters.rho_ini_B;
							}
					
							else 
							{
								rho_imp_R = 0.1 * parameters.rho_ini_R;
								rho_imp_B = 0.9 * parameters.rho_ini_B;
							}
					
							dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
							dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );
						}
					}
				}
			}
		}
		
	}
	else 
	{
		int nx = geometry.nx;
		int ny = geometry.ny;
		int nz = geometry.nz;
		
		int* geo = geometry.ini;
		
		for ( int x = 0; x < nx; x++ )
		{
			for ( int y = 0; y < ny; y++ )
			{
				for ( int z = 0; z < nz; z++ )
				{
					int pto = geo[ x + y * nx + z * nx * ny ] - 1;
					
					if ( pto >= 0 )
					{
						double *f = lattice.inif + ( pto ) * nvel;
								
						double rho_imp = parameters.rho_ini + 0.01 * distribution( generator );		
						
						dist_eq ( f, 0., 0., 0., rho_imp, lattice );
						
						double rho = density ( f );
						
						if ( lattice.inif_m != nullptr )
						{				
							double *f_m = lattice.inif_m + ( pto ) * nvel;

							emite_mediadores ( f_m, rho, lattice );
						}
						if ( lattice.ini_psi != nullptr ) 
						{
							double *psi = lattice.ini_psi;
							
							psi[ x + y * nx + z * nx * ny  ] = rho;
						}
					}
					else if ( lattice.ini_psi != nullptr ) 
					{
						double *psi = lattice.ini_psi;
						
						psi[ x + y * nx + z * nx * ny  ] = parameters.wett_R;
					}					
				}
			}
		}
	}
}
//====================================================================================================================//



//================================ Initial condition using Poisson's Equation ========================================//
//
//      Input: geometry, lattice, rho_ini
//      Output: rho_o
//
//====================================================================================================================//

void initial_rho_poisson ( GEOMETRY geometry, LATTICE lattice , double* rho0, double rho_ini )
{	
	double* vx_0 = new double [ geometry.fluid ];
	double* vy_0 = new double [ geometry.fluid ];
	double* vz_0 = new double [ geometry.fluid ];
		
	for ( int pto = 0; pto < geometry.fluid; pto++ )
	{
		double *f = lattice.inif + ( pto ) * nvel;

		double vx, vy, vz, rho;

		calcula ( f, vx, vy, vz, rho, lattice );

		vx_0[ pto ] = vx;
		vy_0[ pto ] = vy;
		vz_0[ pto ] = vz;
			  
		double delta_rho = 0.0;

		dist_eq_Poisson ( f, vx_0[ pto ], vy_0[ pto ], vz_0[ pto ], delta_rho, rho_ini, lattice );
		
		double mx, my, mz;
		
		propag_site ( lattice, pto, mx, my, mz );
	}
		
	double sum_delta_old = 0.0;
	double sum_delta_new = 0.0;
	double dif_sum_delta = 0.0;
	
	cout << endl;
	
	do
	{				
		double* temp = lattice.inif;
		lattice.inif = lattice.inif_new;
		lattice.inif_new = temp;
		
		#pragma omp parallel for reduction(+:sum_delta_new)

		for ( int pto = 0; pto < geometry.fluid; pto++ )
		{
			double *f = lattice.inif + ( pto ) * nvel;
			
			double delta_rho = density ( f );
						
			dist_eq_Poisson ( f, vx_0[ pto ], vy_0[ pto ], vz_0[ pto ], delta_rho, rho_ini, lattice );
						 
			double mx, my, mz;
						 
			propag_site ( lattice, pto, mx, my, mz );
			
			sum_delta_new = sum_delta_new + delta_rho;
		}
		
		dif_sum_delta = fabs( sum_delta_new - sum_delta_old );
		
		sum_delta_old = sum_delta_new;
		
		sum_delta_new = 0.0;
		
		if ( rand() % 100 == 0 ) cout << "\rdif_sum_delta = " << dif_sum_delta;
		
	} while ( dif_sum_delta > 1.e-18 );
	
	cout << endl;
	
	for ( int pto = 0; pto < geometry.fluid; pto++ )
	{
		double *f = lattice.inif + ( pto ) * nvel;
		
		double delta_rho = density ( f );
						
		dist_eq_Poisson ( f, vx_0[ pto ], vy_0[ pto ], vz_0[ pto ], delta_rho, rho_ini, lattice );
		
		rho0[ pto ] = rho_ini + density( f );		
	}			
	
	rec_scalar_field ( geometry, "rho_0.vtk", rho0 );
	
	delete [] vx_0;
	delete [] vy_0;
	delete [] vz_0;
}
//====================================================================================================================//




//================================ Initial condition flat in x direction (immiscible) ================================//
//
//      Input: geometry, lattice, parameters
//      Output: 
//
//====================================================================================================================//

void initial_conditions_flat_x ( GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	int *geo = geometry.ini;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{	
				int pos = x + y * nx + z * nx * ny;	
				
				if ( geo[pos] )
				{
					double *f_R = lattice.inif_R + ( geo[pos] - 1 ) * nvel;
					double *f_B = lattice.inif_B + ( geo[pos] - 1 ) * nvel;	
				
					double rho_imp_R = 0.0;		
					double rho_imp_B = 0.0;				

					if ( x < nx / 2 ) rho_imp_R = parameters.rho_ini_R;
				
					else rho_imp_B = parameters.rho_ini_B;				
				
					dist_eq ( f_R, 0., 0., 0., rho_imp_R, lattice );
					dist_eq ( f_B, 0., 0., 0., rho_imp_B, lattice );
			
					double rho_R = density ( f_R );
					double rho_B = density ( f_B );
					
					if ( lattice.inif_m != nullptr )
					{
						double phase = rho_R / ( rho_R + rho_B ); 
			
						double *f_m = lattice.inif_m + ( geo[pos] - 1  ) * nvel;

						emite_mediadores ( f_m, phase, lattice );
					}
					
					if ( lattice.ini_psi != nullptr ) 
					{
						double phase = ( rho_R - rho_B ) / ( rho_R + rho_B ); 
						
						double *psi = lattice.ini_psi + pos;
  
						psi[0] = phase;
					}
				}
			}
		}
	}
}
//====================================================================================================================//



//=============================== Initial condition - random phases (MHL model) ======================================//
//
//      Input: fraction of sites started in the phase phi = 1, geometry, lattice, parameters
//      Output: hydrodynamic population, phase field population and auxiliary fields
//
//      Modelo Montessori-Hegele-Lauricella: a rede hidrodinamica parte do repouso com a densidade
//      lida em data_in.txt, e o campo de fase recebe 0 ou 1 sorteado sitio a sitio. E a condicao
//      de decomposicao espinodal -- os dominios se formam e coalescem sozinhos, o que serve bem
//      para ver o modelo funcionando.
//
//      Como phi so aparece atraves de g, basta g_i = w_i * phi: a soma devolve phi e o primeiro
//      momento e nulo, que e o equilibrio de ordem zero com velocidade nula.
//
//====================================================================================================================//

void initial_conditions_MHL_random ( double conc_R, GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	default_random_engine generator;

	uniform_real_distribution<double> distribution( 0.0, 1.0 );

	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	const double rho_ini = parameters.rho_ini;

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int pos = x + y * nx + z * ny * nx;

				if ( geo[pos] == 0 ) continue;			// sitio solido

				int pto = geo[pos] - 1;

				double phi_0 = ( distribution( generator ) < conc_R ) ? 1.0 : 0.0;

				double *f = lattice.inif_R + pto * nvel;

				dist_eq ( f, 0.0, 0.0, 0.0, rho_ini, lattice );

				double *g = lattice.inif_m + pto * nvel;

				for ( int i = 0; i < nvel; i++ ) g[i] = lattice.w[i] * phi_0;

				if ( lattice.ini_phi != nullptr ) lattice.ini_phi[pto] = phi_0;

				if ( lattice.ini_grad != nullptr )
				{
					for ( int a = 0; a < dim; a++ ) lattice.ini_grad[ pto * dim + a ] = 0.0;
				}
			}
		}
	}
}

//====================================================================================================================//




//=============================== Initial condition of a sphere (MHL model) ==========================================//
//
//      Input: radius, interface width, position of the centre, geometry, lattice, parameters
//      Output: hydrodynamic population, phase field population and auxiliary fields
//
//      Gota de fase phi = 1 imersa em phi = 0, com o perfil de equilibrio do modelo,
//
//              phi(r) = 0.5 * [ 1 + tanh( 2 ( R - r ) / W ) ]
//
//      A largura de equilibrio e W = 2/gamma, com gamma = parameters.A_fact: partindo ja com essa
//      largura a gota quase nao precisa relaxar. E a condicao inicial do teste da lei de Laplace,
//      em que a diferenca de pressao deve sair 2*sigma/R.
//
//      O centro entra como double para permitir meio sitio: num dominio de lado par o centro
//      geometrico e (n-1)/2, e arredondar quebraria a simetria da gota.
//
//====================================================================================================================//

void initial_conditions_MHL_sphere ( double radius, double width, double x0, double y0, double z0,
                                     GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	const double rho_ini = parameters.rho_ini;

	//  Largura nao informada (zero ou negativa): usa a de equilibrio, W = 2/gamma.

	double W = width;

	if ( W <= 0.0 ) W = ( parameters.A_fact > 0.0 ) ? 2.0 / parameters.A_fact : 4.0;

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int pos = x + y * nx + z * ny * nx;

				if ( geo[pos] == 0 ) continue;			// sitio solido

				int pto = geo[pos] - 1;

				double dx = x - x0;
				double dy = y - y0;
				double dz = z - z0;

				double r = sqrt ( dx * dx + dy * dy + dz * dz );

				double phi_0 = 0.5 * ( 1.0 + tanh ( 2.0 * ( radius - r ) / W ) );

				double *f = lattice.inif_R + pto * nvel;

				dist_eq ( f, 0.0, 0.0, 0.0, rho_ini, lattice );

				double *g = lattice.inif_m + pto * nvel;

				for ( int i = 0; i < nvel; i++ ) g[i] = lattice.w[i] * phi_0;

				if ( lattice.ini_phi != nullptr ) lattice.ini_phi[pto] = phi_0;

				if ( lattice.ini_grad != nullptr )
				{
					for ( int a = 0; a < dim; a++ ) lattice.ini_grad[ pto * dim + a ] = 0.0;
				}
			}
		}
	}
}

//=============================== Condicao inicial: laje plana de vermelho ===========================================//
//
//      Irma unidimensional de initial_conditions_drop().  Numa caixa periodica de nx sitios, o
//      fluido vermelho ocupa a metade  [ 0 , nx/2 )  e o azul a outra metade, o que cria DUAS
//      interfaces planas normais a x: uma em  x = nx/2  e outra em  x = 0 , imposta pela
//      periodicidade.  A medida ( calc_flat_tension ) integra so a central.
//
//      A geometria e a mesma da gota, com o raio trocado por meia laje: chamando
//      x0 = nx/4  o centro do vermelho e  meia = nx/4 ,
//
//          phi( x ) = 1 / ( 1 + exp[ k_perfil ( d - meia ) ] ) ,    d = distancia periodica a x0
//
//      e as duas interfaces saem simetricas, sem emenda.  k_perfil = 0 da o degrau abrupto.
//
//      A espessura de equilibrio nao e esta: quem a fixa e a recoloracao ( beta ).  k_perfil serve
//      so para pular o transiente acustico que um degrau lanca numa caixa periodica.
//
//      fase_simetrica e eq_modo tem o mesmo significado de initial_conditions_drop().
//
//      Input: k_perfil, delta_rho, velocidade inicial, geometry, lattice, parameters
//      Output: inif_R, inif_B, inif_m, ini_psi
//
//====================================================================================================================//

void initial_conditions_flat ( double k_perfil, double delta_rho, double vx, double vy, double vz,
							GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters,
							bool fase_simetrica, int eq_modo )
{
	int *geo = geometry.ini;

	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	const double x0   = 0.25 * ( double ) nx;		// centro da laje vermelha
	const double meia = 0.25 * ( double ) nx;		// meia espessura da laje

	for ( int x = 0; x < nx; x++ )
	{
		//  Distancia periodica ao centro da laje.

		double d = ( double ) x - x0;

		if ( d >  0.5 * nx ) d -= nx;
		if ( d < -0.5 * nx ) d += nx;

		d = fabs ( d );

		double phi;

		if ( k_perfil > 0.0 )
		{
			double arg = k_perfil * ( d - meia );

			if ( arg >  60.0 ) arg =  60.0;
			if ( arg < -60.0 ) arg = -60.0;

			phi = 1.0 / ( 1.0 + exp ( arg ) );
		}

		else phi = ( d < meia ) ? 1.0 : 0.0;

		const double rho_tot = parameters.rho_ini_B + phi * ( parameters.rho_ini_R - parameters.rho_ini_B )
							 + phi * delta_rho;

		const double rho_R = phi * rho_tot;
		const double rho_B = ( 1.0 - phi ) * rho_tot;

		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				const int pos = x + y * nx + z * nx * ny;

				if ( geo[pos] == 0 )
				{
					if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = parameters.wett_R;

					continue;
				}

				const int pto = geo[pos] - 1;

				double *fR = lattice.inif_R + pto * nvel;
				double *fB = lattice.inif_B + pto * nvel;

				//  Mesmos tres modos de equilibrio de initial_conditions_drop().

				if ( eq_modo == EQ_ALFA )
				{
					dist_eq_LYL ( fR, vx, vy, vz, rho_R, parameters.alfa_R, 0.0, lattice );
					dist_eq_LYL ( fB, vx, vy, vz, rho_B, parameters.alfa_B, 0.0, lattice );
				}

				else if ( eq_modo == EQ_SAITO )
				{
					const double p_t = rho_R * parameters.cs2_R + rho_B * parameters.cs2_B;

					double f_eq[nvel];

					dist_eq_saito ( f_eq, vx, vy, vz, rho_tot, p_t, lattice );

					const double cR = ( rho_tot > 0.0 ) ? rho_R / rho_tot : 0.0;

					for ( int i = 0; i < nvel; i++ )
					{
						fR[i] = cR * f_eq[i];
						fB[i] = ( 1.0 - cR ) * f_eq[i];
					}
				}

				else
				{
					dist_eq ( fR, vx, vy, vz, rho_R, lattice );
					dist_eq ( fB, vx, vy, vz, rho_B, lattice );
				}

				const double fase_med = fase_simetrica ? ( 2.0 * phi - 1.0 ) : phi;

				if ( lattice.inif_m != nullptr )
					emite_mediadores ( lattice.inif_m + pto * nvel, fase_med, lattice );

				if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = 2.0 * phi - 1.0;
			}
		}
	}
}

//====================================================================================================================//


//=============================== Initial Condition of a Drop with a diffuse interface ===============================//
//
//      Gota de fluido vermelho ( raio 'radius', centro x0,y0,z0 ) num banho de fluido azul, ja com o
//      perfil de interface do modelo, em vez de um degrau:
//
//          phi(r) = 1 / ( 1 + exp( k_perfil ( r - radius ) ) )         phi = rho_R / ( rho_R + rho_B )
//
//      Para os modelos com mediadores ( Santos imiscivel ), o perfil de equilibrio e a logistica de
//      inclinacao  k = lambda_ef / cs^2 = 3 lambda_ef ( D3Q19 ), com
//
//          lambda_ef = lambda / ( 1 - 1 / ( 2 tau_m ) )        lambda = parameters.A_fact
//
//      Partir dai poupa o transiente acustico do degrau, que numa caixa periodica reverbera por
//      milhares de passos.  k_perfil <= 0 volta a dar o degrau ( equivale a initial_conditions_sphere ).
//
//      'cilindro' = true trata a gota como um cilindro de eixo z ( o caso 2D, nz = 1 ): a distancia
//      e medida so em x e y.  Com false a gota e uma esfera.
//
//      'delta_rho' e um acrescimo de densidade DENTRO da gota, para ja entrar com o salto de Laplace
//      ( delta_rho = delta_p / cs^2 ).  ATENCAO: densidade a mais e massa a mais, e a massa se
//      conserva -- a gota de equilibrio fica MAIOR que 'radius'.  delta_rho = 0 e o caso normal.
//
//      'fase_simetrica' escolhe a convencao dos mediadores: false ( padrao ) semeia
//      phi = rho_R / rho, em [0,1], que e a do modelo de Santos; true semeia
//      rho^N = ( rho_R - rho_B ) / rho, em [-1,1], que e a de Spencer-Halliday-Care.
//
//      Input : radius, k_perfil, delta_rho, x0, y0, z0, vx, vy, vz, cilindro, geometry, lattice,
//              parameters, fase_simetrica
//      Output: distribution functions ( inif_R, inif_B ), mediators ( inif_m ), scalar field ( ini_psi )
//
//====================================================================================================================//

void initial_conditions_drop ( double radius, double k_perfil, double delta_rho, int x0, int y0, int z0,
							double vx, double vy, double vz, bool cilindro, GEOMETRY geometry, LATTICE lattice,
							PARAMETERS parameters, bool fase_simetrica, int eq_modo )
{
	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				int pos = x + y * nx + z * nx * ny;

				//  No solido nao ha populacao, mas o campo escalar existe em todos os voxels e e
				//  lido pelos vizinhos em gradient(): recebe a molhabilidade, como fazem as outras
				//  condicoes iniciais deste arquivo.

				if ( geo[pos] == 0 )
				{
					if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = parameters.wett_R;

					continue;
				}

				double dx = ( double ) ( x - x0 );
				double dy = ( double ) ( y - y0 );
				double dz = cilindro ? 0.0 : ( double ) ( z - z0 );

				double r = sqrt ( dx * dx + dy * dy + dz * dz );

				//  Perfil.  O expoente e limitado para nao estourar exp() longe da interface.

				double phi;

				if ( k_perfil > 0.0 )
				{
					double arg = k_perfil * ( r - radius );

					if ( arg >  60.0 ) arg =  60.0;
					if ( arg < -60.0 ) arg = -60.0;

					phi = 1.0 / ( 1.0 + exp ( arg ) );
				}

				else phi = ( r < radius ) ? 1.0 : 0.0;

				//  A densidade total interpola entre os dois fluidos puros; delta_rho e o salto de
				//  Laplace, que so existe dentro da gota.

				double rho_tot = parameters.rho_ini_B + phi * ( parameters.rho_ini_R - parameters.rho_ini_B )
							   + phi * delta_rho;

				int pto = geo[pos] - 1;

				//  eq_modo escolhe COM QUAL EQUILIBRIO as populacoes sao semeadas.  Importa quando
				//  o modelo nao usa o equilibrio classico: com w_i a pressao inicial sairia rho/3
				//  em vez da pressao do modelo, e a caixa inteira arrancaria com uma onda acustica
				//  que reverbera por milhares de passos numa caixa periodica.
				//
				//      EQ_CLASSICO : dist_eq , o de sempre;
				//      EQ_ALFA     : equilibrio de alpha variavel ( eq. (8) de Liu-Valocchi-Kang,
				//                    que e dist_eq_LYL com hi_order = 0 ) -- so difere se alpha != 1/3;
				//      EQ_SAITO    : equilibrio generalizado de Saito ( eq. 56 ), pelos momentos
				//                    centrais.  Precisa da pressao, p = rho_R cs2_R + rho_B cs2_B.

				double *fR = lattice.inif_R + pto * nvel;
				double *fB = lattice.inif_B + pto * nvel;

				const double rho_R = phi * rho_tot;
				const double rho_B = ( 1.0 - phi ) * rho_tot;

				if ( eq_modo == EQ_ALFA )
				{
					dist_eq_LYL ( fR, vx, vy, vz, rho_R, parameters.alfa_R, 0.0, lattice );
					dist_eq_LYL ( fB, vx, vy, vz, rho_B, parameters.alfa_B, 0.0, lattice );
				}

				else if ( eq_modo == EQ_SAITO )
				{
					//  Em Saito o equilibrio e da MISTURA: um so f^eq, com rho e p totais, repartido
					//  entre as duas cores pela concentracao.  E o que a recoloracao pressupoe.

					const double rho_t = rho_R + rho_B;

					const double p_t = rho_R * parameters.cs2_R + rho_B * parameters.cs2_B;

					double f_eq[nvel];

					dist_eq_saito ( f_eq, vx, vy, vz, rho_t, p_t, lattice );

					const double cR = ( rho_t > 0.0 ) ? rho_R / rho_t : 0.0;

					for ( int i = 0; i < nvel; i++ )
					{
						fR[i] = cR * f_eq[i];
						fB[i] = ( 1.0 - cR ) * f_eq[i];
					}
				}

				else
				{
					dist_eq ( fR, vx, vy, vz, rho_R, lattice );
					dist_eq ( fB, vx, vy, vz, rho_B, lattice );
				}

				//  Convencoes do resto do arquivo: ini_psi guarda ( rho_R - rho_B ) / rho, em [-1,1].
				//  Ja os mediadores mudam de modelo para modelo:
				//
				//      fase_simetrica = false  ( Santos ) :  phi = rho_R / rho , em [0,1]
				//      fase_simetrica = true   ( SHC )    :  rho^N = ( rho_R - rho_B ) / rho , em [-1,1]
				//
				//  Nos dois casos os mediadores sao reemitidos do zero a cada passo, entao esta
				//  semeadura so importa para um diagnostico feito antes da primeira colisao.

				double fase_med = fase_simetrica ? ( 2.0 * phi - 1.0 ) : phi;

				if ( lattice.inif_m != nullptr )
				{
					double *f_m = lattice.inif_m + pto * nvel;

					emite_mediadores ( f_m, fase_med, lattice );
				}

				if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = 2.0 * phi - 1.0;
			}
		}
	}
}

//====================================================================================================================//


//=============================== Initial Condition of a Capillary Wave ==============================================//
//
//      Interface plana normal a x, deslocada por uma cossenoide ao longo de y:
//
//          x_int( y ) = x_int  +  amplitude * cos( 2 pi modo y / ny )
//
//      O vermelho fica em x < x_int( y ) e o azul em x > x_int( y ), com o mesmo perfil logistico
//      de initial_conditions_flat:
//
//          phi( x, y ) = 1 / ( 1 + exp[ k_perfil ( x - x_int( y ) ) ] )
//
//      HA UMA INTERFACE SO.  A caixa NAO e periodica em x -- as duas faces x levam derivada nula
//      ( devnull_RB_x ) -- e e periodica em y, que e a direcao de propagacao.  Por isso a distancia
//      aqui nao e periodica em x, ao contrario de initial_conditions_flat.
//
//      x_int <= 0 poe a interface no meio da caixa.  k_perfil <= 0 da o degrau.
//
//      Input : amplitude, modo, x_int, k_perfil, vx, vy, vz, geometry, lattice, parameters,
//              fase_simetrica, eq_modo
//      Output: distribution functions ( inif_R, inif_B ), mediators ( inif_m ), scalar field ( ini_psi )
//
//====================================================================================================================//

void initial_conditions_wave ( double amplitude, int modo, double x_int, double k_perfil,
							double vx, double vy, double vz,
							GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters,
							bool fase_simetrica, int eq_modo )
{
	int *geo = geometry.ini;

	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	const double PI = 3.14159265358979323846;

	const double x0 = ( x_int > 0.0 ) ? x_int : 0.5 * ( double ) nx;

	const double k_onda = 2.0 * PI * ( double ) modo / ( double ) ny;

	for ( int y = 0; y < ny; y++ )
	{
		//  Posicao da interface nesta coluna.

		const double x_i = x0 + amplitude * cos ( k_onda * ( double ) y );

		for ( int x = 0; x < nx; x++ )
		{
			//  d < 0 do lado vermelho, d > 0 do lado azul.  Sem distancia periodica em x.

			const double d = ( double ) x - x_i;

			double phi;

			if ( k_perfil > 0.0 )
			{
				double arg = k_perfil * d;

				if ( arg >  60.0 ) arg =  60.0;
				if ( arg < -60.0 ) arg = -60.0;

				phi = 1.0 / ( 1.0 + exp ( arg ) );
			}

			else phi = ( d < 0.0 ) ? 1.0 : 0.0;

			const double rho_tot = parameters.rho_ini_B + phi * ( parameters.rho_ini_R - parameters.rho_ini_B );

			const double rho_R = phi * rho_tot;
			const double rho_B = ( 1.0 - phi ) * rho_tot;

			for ( int z = 0; z < nz; z++ )
			{
				const int pos = x + y * nx + z * nx * ny;

				if ( geo[pos] == 0 )
				{
					if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = parameters.wett_R;

					continue;
				}

				const int pto = geo[pos] - 1;

				double *fR = lattice.inif_R + pto * nvel;
				double *fB = lattice.inif_B + pto * nvel;

				//  Mesmos tres modos de equilibrio de initial_conditions_flat().

				if ( eq_modo == EQ_ALFA )
				{
					dist_eq_LYL ( fR, vx, vy, vz, rho_R, parameters.alfa_R, 0.0, lattice );
					dist_eq_LYL ( fB, vx, vy, vz, rho_B, parameters.alfa_B, 0.0, lattice );
				}

				else if ( eq_modo == EQ_SAITO )
				{
					const double p_t = rho_R * parameters.cs2_R + rho_B * parameters.cs2_B;

					double f_eq[nvel];

					dist_eq_saito ( f_eq, vx, vy, vz, rho_tot, p_t, lattice );

					const double cR = ( rho_tot > 0.0 ) ? rho_R / rho_tot : 0.0;

					for ( int i = 0; i < nvel; i++ )
					{
						fR[i] = cR * f_eq[i];
						fB[i] = ( 1.0 - cR ) * f_eq[i];
					}
				}

				else
				{
					dist_eq ( fR, vx, vy, vz, rho_R, lattice );
					dist_eq ( fB, vx, vy, vz, rho_B, lattice );
				}

				const double fase_med = fase_simetrica ? ( 2.0 * phi - 1.0 ) : phi;

				if ( lattice.inif_m != nullptr )
					emite_mediadores ( lattice.inif_m + pto * nvel, fase_med, lattice );

				if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = 2.0 * phi - 1.0;
			}
		}
	}
}

//====================================================================================================================//


//=============================== Read the inicialization file =======================================================//
//
//      Input:
//      Output: geometry file name, pixel dimension, steps, files, relaxation time,
//              initial density
//
//====================================================================================================================//

void read_data ( GEOMETRY& geometry, PARAMETERS& parameters )
{
    //----------------------------------------------------------------------------------------------------------------//

    string name_in = "data_in.txt";

    ifstream f_in( name_in );

	//----------------------------------------------------------------------------------------------------------------//

    string name_out = "dat_out.txt";

    cout << "\nName of the output file: " << name_out << endl;

    //----------------------------------------------------------------------------------------------------------------//

    ofstream fdat( name_out );

    fdat << "Name of the output file: " << name_out << endl;

    //----------------------------------------------------------------------------------------------------------------//

    f_in >> geometry.file;

    cout << "\nName of the geometry file: " << geometry.file << endl;
    fdat << "\nName of the geometry file: " << geometry.file << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_ftesc;
    
    f_in >> st_ftesc;

    f_in >> geometry.ftesc;  // read the pixel's dimension ( m )

    cout << "\nPixel dimension = " << geometry.ftesc << " m" << endl; 
    fdat << "\nPixel dimension = " << geometry.ftesc << " m" << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_steps;
    
    f_in >> st_steps;

    f_in >> parameters.n_steps;

    cout << "\nNumber of steps: " << parameters.n_steps << endl;
    fdat << "\nNumber of steps: " << parameters.n_steps << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_files;
    
    f_in >> st_files;

    f_in >> parameters.n_files;

    cout << "\nNumber of files: " << parameters.n_files << endl;
    fdat << "\nNumber of files: " << parameters.n_files << endl;
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_threads;
    
    f_in >> st_threads;

    f_in >> parameters.n_threads;

    cout << "\nNumber of threads: " << parameters.n_threads << endl;
    fdat << "\nNumber of threads: " << parameters.n_threads << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_tau;
    
    f_in >> st_tau;

    f_in >> parameters.tau;
    
    parameters.visc = (1./3.) * ( parameters.tau - 0.5 );

    cout << "\nRelaxation time: " << parameters.tau << endl;
    fdat << "\nRelaxation time: " << parameters.tau << endl;
    
    cout << "\nKinematic viscosity = " << parameters.visc << endl;
    fdat << "\nKinematic viscosity = " << parameters.visc << endl;
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_nh;
    
    f_in >> st_tau_nh;

    f_in >> parameters.tau_nh;

    cout << "\nRelaxation time (non-hidrodynamic): " << parameters.tau_nh << endl;
    fdat << "\nRelaxation time (non-hidrodynamic): " << parameters.tau_nh << endl;

	//----------------------------------------------------------------------------------------------------------------//
	
	string st_rho;
    
    f_in >> st_rho;

    f_in >> parameters.rho_ini;

    cout << "\nDensity: " << parameters.rho_ini << endl;
    fdat << "\nDensity: " << parameters.rho_ini << endl;
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_R;
    
    f_in >> st_tau_R;

    f_in >> parameters.tau_R;

    parameters.visc_R = ( 1.0 / 3.0 ) * ( parameters.tau_R - 0.5 );

    if ( parameters.tau_R )
    {
		cout << "\nRelaxation time (Red): " << parameters.tau_R << " => viscosity R = " << parameters.visc_R << endl;		
		fdat << "\nRelaxation time (Red): " << parameters.tau_R << " => viscosity R = " << parameters.visc_R << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_B;
    
    f_in >> st_tau_B;

    f_in >> parameters.tau_B;

    parameters.visc_B = ( 1.0 / 3.0 ) * ( parameters.tau_B - 0.5 );

    if ( parameters.tau_B )
    {
		cout << "\nRelaxation time (Blue): " << parameters.tau_B << " => viscosity B = " << parameters.visc_B << endl;		
		fdat << "\nRelaxation time (Blue): " << parameters.tau_B << " => viscosity B = " << parameters.visc_B << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_m;
    
    f_in >> st_tau_m;

    f_in >> parameters.tau_m;
    
    if ( parameters.tau_R )
    {
		cout << "\nRelaxation time (mixture): " << parameters.tau_m << endl;		
		fdat << "\nRelaxation time (mixture): " << parameters.tau_m << endl;
	}
	
    //----------------------------------------------------------------------------------------------------------------//

    string st_rho_ini_R;
    
    f_in >> st_rho_ini_R;

    f_in >> parameters.rho_ini_R;
    
    if ( parameters.tau_R )
    {
		cout << "\nInitial density (Red): " << parameters.rho_ini_R << endl;		
		fdat << "\nInitial density (Red): " << parameters.rho_ini_R << endl;
	}
	
	//----------------------------------------------------------------------------------------------------------------//

    string st_rho_ini_B;
    
    f_in >> st_rho_ini_B;

    f_in >> parameters.rho_ini_B;
    
    if ( parameters.tau_B )
    {
		cout << "\nInitial density (Blue): " << parameters.rho_ini_B << endl;		
		fdat << "\nInitial density (Blue): " << parameters.rho_ini_B << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_fat_int;
    
    f_in >> st_fat_int;

    f_in >> parameters.A_fact;
    
    if ( parameters.tau_R )
    {
		cout << "\nInterfacial factor[0;0.4]: " << parameters.A_fact << endl;		
		fdat << "\nInterfacial factor[0;0.4]: " << parameters.A_fact << endl;
	}
	
	//----------------------------------------------------------------------------------------------------------------//

    string st_fat_int_R;
    
    f_in >> st_fat_int_R;

    f_in >> parameters.A_fact_R;
    
    if ( parameters.tau_R )
    {
		cout << "\nInteraction factor (Red) [0;0.4]: " << parameters.A_fact_R << endl;		
		fdat << "\nInteraction factor (Red) [0;0.4]: " << parameters.A_fact_R << endl;
	}
	
	//----------------------------------------------------------------------------------------------------------------//

    string st_fat_int_B;
    
    f_in >> st_fat_int_B;

    f_in >> parameters.A_fact_B;
    
    if ( parameters.tau_B )
    {
		cout << "\nInteraction factor (Blue) [0;0.4]: " << parameters.A_fact_B << endl;		
		fdat << "\nInteraction factor (Blue) [0;0.4]: " << parameters.A_fact_B << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

	string st_recoll;

    f_in >> st_recoll;

    f_in >> parameters.recoll;
    
    if ( parameters.tau_R )
    {
		cout << "\nRecolloring factor[0;1]: " << parameters.recoll << endl;		
		fdat << "\nRecolloring factor[0;1]: " << parameters.recoll << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_wett_R;
    
    f_in >> st_wett_R;
    
    f_in >> parameters.wett_R;
    
    if ( parameters.tau_R )
    {
		cout << "\nRed_wall interaction [0.0 ; 1.0] : " << parameters.wett_R << endl;		
		fdat << "\nRed_wall interaction [0.0 ; 1.0] : " << parameters.wett_R << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    f_in.close();

    fdat.close();
}

//====================================================================================================================//




//=============================== Read the geometry file =============================================================//
//
//      Input: geometry file name, pointer to adress of begining of the geometry
//      Output: number of fluid points
//
//====================================================================================================================//

int read_geo ( string nome_geo, int *meio, int pts_in, int pts_out, bool wall )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nTamanho:  x = " << nx << ";  y = "   << ny << ";  z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] )
                    {
                        poros++;

                        meio[pos] = poros;
                    }
                }
                else if ( x == 0 && wall )
                {
                    meio[pos] = 0;
                }
                
                else if ( x == nx - 1 && wall )
                {
                    meio[pos] = 0;
                }
                
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//====================================================================================================================//




//=============================== Read the geometry file =============================================================//
//
//      Input: geometry file name, pointer to adress of begining of the geometry
//      Output: number of fluid points
//
//====================================================================================================================//

int read_geo ( string nome_geo, int *meio, int pts_in, int pts_out )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nTamanho:  x = " << nx << ";  y = "   << ny << ";  z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] )
                    {
                        poros++;

                        meio[pos] = poros;
                    }
                }
                else if ( x == 0 )
                {
                    meio[pos] = 0;
                }
                
                else if ( x == nx - 1 )
                {
                    meio[pos] = 0;
                }
                
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//====================================================================================================================//



//=============================== Define a D3Q15 lattice =============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d3q15 ( LATTICE& lattice )
{	
    //------------- |ci| = 0 ----------------------//
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;
    lattice.c_i[ i * dim + 2 ] = 0;

    //------------- |ci| = 1 ----------------------//
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

    //------------ |c| = sqrt(3) -----------------//

	i = 7;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
	
	i = 8;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 9;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 10;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 11;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 12;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 13;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 14;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    lattice.w[0] = 2. / 9.;

    for ( int i = 1; i < 7; i++ ) lattice.w[i] = 1. / 9.;

    for ( int i = 7; i < nvel; i++ ) lattice.w[i] = 1. / 72.;
    
    lattice.c_s2 = 1./3.;	// sound velocity
    
    lattice.one_over_c_s2 = 3.; // inverse of the sound velocity
}

//====================================================================================================================//



//=============================== Defines a D3Q19 lattice ============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d3q19 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;
    lattice.c_i[ i * dim + 2 ] = 0;

    //------------- |ci| = 1 ----------------------//
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

    //------------ |c| = sqrt(2) -----------------//

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 9;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 10;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 11;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 12;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 13;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 14;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    i = 15;
    
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 16;
    	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
	i = 17;    
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 18;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    lattice.w[0] = 1. / 3.;

    for ( int i = 1; i < 7; i++ ) lattice.w[i] =  1. / 18.;

    for ( int i = 7; i < nvel; i++ ) lattice.w[i] = 1. / 36.; 
    
    lattice.c_s2 = 1./3.; // sound velocity
    
    lattice.one_over_c_s2 = 3.; // inverse of the sound velocity
}

//====================================================================================================================//



//=============================== Defines a D3Q19 lattice ============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d3q27 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;
    lattice.c_i[ i * dim + 2 ] = 0;

    //------------- |ci| = 1 ----------------------//
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

    //------------ |c| = sqrt(2) -----------------//

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 9;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 10;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 11;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 12;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 13;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 14;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    i = 15;
    
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 16;
    	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
	i = 17;    
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 18;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    //------------ |c| = sqrt(3) -----------------//

	i = 19;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 20;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 21;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 22;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 23;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 24;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 25;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 26;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    lattice.w[0] = 8. / 27.;

    for ( int i = 1; i < 7; i++ ) lattice.w[i] =  2. / 27.;

    for ( int i = 7; i < 19; i++ ) lattice.w[i] = 1. / 54.; 
    
    for ( int i = 19; i < nvel; i++ ) lattice.w[i] = 1. / 216.; 
    
    lattice.c_s2 = 1./3.; // sound velocity
    
    lattice.one_over_c_s2 = 3.; // inverse of the sound velocity
}

//====================================================================================================================//




//=============================== Defines a D3Q77 lattice ============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d3q77 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;
    lattice.c_i[ i * dim + 2 ] = 0;

    //------------- |ci| = 1 ----------------------//
    
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

    //------------ |c| = sqrt(2) -----------------//

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 9;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 10;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 11;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 12;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 13;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 14;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    i = 15;
    
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 16;
    	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
	i = 17;    
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 18;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

	//------------ |c| = sqrt(3) -----------------//

	i = 19;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 20;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 21;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 22;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 23;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 24;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 25;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 26;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    //------------ |c| = sqrt(6) -----------------//

	i = 27;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 28;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 29;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 30;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 31;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 32;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 33;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 34;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    i = 35;
    
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 36;
    	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
	i = 37;    
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 38;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 39;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 40;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 41;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 42;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 43;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 44;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 45;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 46;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    i = 47;
    
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 48;
    	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
	i = 49;    
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 50;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    //------------ |c| = sqrt(8) -----------------//

	i = 51;

    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 52;

    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 53;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 54;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 55;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 56;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 57;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 58;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  2;
    
    i = 59;
    
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] = -2;
    
    i = 60;
    	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  2;
    
	i = 61;    
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 62;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] = -2;

	//------------- |ci| = 3 ----------------------//
    
	i = 63;
	
    lattice.c_i[ i * dim + 0 ] =  3;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 64;
	
    lattice.c_i[ i * dim + 0 ] = -3;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 65;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  3;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 66;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -3;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 67;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  3;

	i = 68;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -3;
    
    //------------ |c| = sqrt(12) -----------------//

	i = 69;

    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 70;

    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 71;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 72;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 73;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  2;

	i = 74;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] = -2;
    
    i = 75;
	
    lattice.c_i[ i * dim + 0 ] =  2;
    lattice.c_i[ i * dim + 1 ] =  2;
    lattice.c_i[ i * dim + 2 ] = -2;

	i = 76;
	
    lattice.c_i[ i * dim + 0 ] = -2;
    lattice.c_i[ i * dim + 1 ] = -2;
    lattice.c_i[ i * dim + 2 ] =  2;
	 
    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    lattice.w[0] = 5221./20480.;
    
	if ( nvel == 77 )
	{
		for ( int i = 1; i < 7; i++ ) lattice.w[i] = 2907./40960.;

		for ( int i = 7; i < 19; i++ ) lattice.w[i] = 939./40960.;
		
		for ( int i = 19; i < 27; i++ ) lattice.w[i] = 13./2560.;
		
		for ( int i = 27; i < 51; i++ ) lattice.w[i] = 1./16384.; 
		
		for ( int i = 51; i < 63; i++ ) lattice.w[i] = 3./81920.; 
		
		for ( int i = 63; i < 69; i++ ) lattice.w[i] = 31./122880.; 
		
		for ( int i = 69; i < nvel; i++ ) lattice.w[i] = 1./81920.; 
	}
	else cout << "\nErro! Verificar nvel = 77" << endl;

    lattice.c_s2 = 3./8.; // sound velocity
    
    lattice.one_over_c_s2 = 8./3.; // inverse of the sound velocity       
}

//====================================================================================================================//




//=============================== Defines a D2Q9 lattice =============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d2q9 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;

    //------------- |ci| = 1 ----------------------//
    
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;

	//------------ |c| = sqrt(2) -----------------//

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    lattice.w[0] = 4. / 9.;

    for ( int i = 1; i < 5; i++ ) lattice.w[i] =  1. / 9.;

    for ( int i = 5; i < nvel; i++ ) lattice.w[i] = 1. / 36.; 
    
    lattice.c_s2 = 1./3.; // sound velocity
    
    lattice.one_over_c_s2 = 3.; // inverse of the sound velocity
}

//====================================================================================================================//




//=============================== Defines a D2V37 lattice ============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d2v37 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;

    //------------- |ci| = 1 ----------------------//
    
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;

	//------------ |c| = sqrt(2) -----------------//

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    
	//------------ |c| = 2 -----------------------//

	i = 9;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 10;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 11;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] =  2;

	i = 12;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] = -2;
	
	//------------ |c| = sqrt(5) -----------------//

	i = 13;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] =  1;

	i = 14;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] = -1;

	i = 15;
	
	lattice.c_i[ i * dim + 0 ] =  1;	
	lattice.c_i[ i * dim + 1 ] =  2;

	i = 16;
	
	lattice.c_i[ i * dim + 0 ] = -1;	
	lattice.c_i[ i * dim + 1 ] = -2;

	i = 17;
	
	lattice.c_i[ i * dim + 0 ] = -1;	
	lattice.c_i[ i * dim + 1 ] =  2;

	i = 18;
	
	lattice.c_i[ i * dim + 0 ] =  1;	
	lattice.c_i[ i * dim + 1 ] = -2;

	i = 19;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] =  1;

	i = 20;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] = -1;

	//------------ |c| = sqrt(8) -----------------//

	i = 21;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] =  2;

	i = 22;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] = -2;

	i = 23;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] = -2;

	i = 24;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] =  2;
	
	//------------ |c| = 3 -----------------------//

	i = 25;
	
	lattice.c_i[ i * dim + 0 ] =  3;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 26;
	
	lattice.c_i[ i * dim + 0 ] = -3;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 27;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] =  3;

	i = 28;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] = -3;

	//------------ |c| = sqrt(10) ----------------//

	i = 29;
	
	lattice.c_i[ i * dim + 0 ] =  3;	
	lattice.c_i[ i * dim + 1 ] =  1;

	i = 30;
	
	lattice.c_i[ i * dim + 0 ] = -3;	
	lattice.c_i[ i * dim + 1 ] = -1;

	i = 31;
	
	lattice.c_i[ i * dim + 0 ] =  1;	
	lattice.c_i[ i * dim + 1 ] =  3;

	i = 32;
	
	lattice.c_i[ i * dim + 0 ] = -1;	
	lattice.c_i[ i * dim + 1 ] = -3;

	i = 33;
	
	lattice.c_i[ i * dim + 0 ] = -1;	
	lattice.c_i[ i * dim + 1 ] =  3;

	i = 34;
	
	lattice.c_i[ i * dim + 0 ] =  1;	
	lattice.c_i[ i * dim + 1 ] = -3;

	i = 35;
	
	lattice.c_i[ i * dim + 0 ] = -3;	
	lattice.c_i[ i * dim + 1 ] =  1;

	i = 36;
	
	lattice.c_i[ i * dim + 0 ] =  3;	
	lattice.c_i[ i * dim + 1 ] = -1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//
    
     double a_2 = (1./36) * ( 49. - (17.* pow(7.,(2./3.))) / pow(67.+36*sqrt(30.),(1./3.))
                + pow((469 + 252 * sqrt(30.)),(1./3.)) );
    
    lattice.c_s2 = 1. / a_2; // sound velocity squared
    
    lattice.one_over_c_s2 = a_2; // inverse of the sound velocity squared
    
    double a_4 = a_2 * a_2;
    
    double a_8 = a_4 * a_4;
  
    double w_0 = ( -2827 +5330 * a_2 -3023 * a_4 + 432 * a_8 ) / ( 432 * a_8 );
    
	double w_1 = ( 93 - 174 * a_2 + 92 * a_4 ) / ( 72 * a_8 );
	
	double w_2 = ( 81 - 141 * a_2 + 76 * a_4 ) / ( 144 * a_8 );
		
	double w_3 = ( 7. ) / ( 240 * a_4 );	
	
	double w_4 = -( 45 - 60 * a_2 + 16 * a_4 ) / ( 360 * a_8 );	
	
	double w_5 = ( -3. + a_2 )*( -3. + a_2 ) / ( 576 * a_8 );

	double w_6 = -( 5 - 10 * a_2 + 4 * a_4 ) / ( 1080 * a_8 );	
	
	double w_7 = ( 15 - 15 * a_2 + 4 * a_4 ) / ( 1440 * a_8 );	
	
    lattice.w[0] = w_0;
    
    if ( nvel == 37 )
	{
		for ( int i = 1; i < 5; i++ ) lattice.w[i] = w_1;

		for ( int i = 5; i < 9; i++ ) lattice.w[i] = w_2; 
		
		for ( int i = 9; i < 13; i++ ) lattice.w[i] = w_3; 
		
		for ( int i = 13; i < 21; i++ ) lattice.w[i] = w_4; 
		
		for ( int i = 21; i < 25; i++ ) lattice.w[i] = w_5; 
		
		for ( int i = 25; i < 29; i++ ) lattice.w[i] = w_6; 
		
		for ( int i = 29; i < nvel; i++ ) lattice.w[i] = w_7; 
	}
	else cout << "\nErro! Verificar nvel = 37" << endl;
       
}

//====================================================================================================================//




//=============================== Defines a D2V37 lattice ============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d2v25 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;

    //------------- |ci| = 1 ----------------------//
    
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;

	//------------ |c| = sqrt(2) -----------------//

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    
	//------------ |c| = 2 -----------------------//

	i = 9;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 10;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 11;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] =  2;

	i = 12;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] = -2;
	
	//------------ |c| = sqrt(8) -----------------//

	i = 13;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] =  2;

	i = 14;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] = -2;

	i = 15;
	
	lattice.c_i[ i * dim + 0 ] =  2;	
	lattice.c_i[ i * dim + 1 ] = -2;

	i = 16;
	
	lattice.c_i[ i * dim + 0 ] = -2;	
	lattice.c_i[ i * dim + 1 ] =  2;

	//------------ |c| = 3 -----------------------//

	i = 17;
	
	lattice.c_i[ i * dim + 0 ] =  3;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 18;
	
	lattice.c_i[ i * dim + 0 ] = -3;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 19;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] =  3;

	i = 20;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] = -3;

	
	//------------ |c| = 4 -----------------------//

	i = 21;
	
	lattice.c_i[ i * dim + 0 ] =  4;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 22;
	
	lattice.c_i[ i * dim + 0 ] = -4;	
	lattice.c_i[ i * dim + 1 ] =  0;

	i = 23;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] =  4;

	i = 24;
	
	lattice.c_i[ i * dim + 0 ] =  0;	
	lattice.c_i[ i * dim + 1 ] = -4;
            
    //-------------------------------- Weights of the Lattice --------------------------------//
    
    double sqr33 = sqrt( 33 );
    
    double pot4 = pow( ( 15 - sqr33 ), 4 );
    
    lattice.w[0] = 16 * ( 6849 - 1135 * sqr33 ) / ( 3 * pot4 );
    
    if ( nvel  == 25 )
    {    
		for ( int i = 1; i < 5; i++ ) lattice.w[i] = 64 * ( 2619 - 437 * sqr33 ) / ( 15 * pot4 );
		
		for ( int i = 5; i < 9; i++ ) lattice.w[i] = 512 * ( 7 - sqr33 ) / pot4;
		
		for ( int i = 9; i < 13; i++ ) lattice.w[i] = 8 * ( 159 + 47 * sqr33 ) / ( 15 * pot4 );
		
		for ( int i = 13; i < 17; i++ ) lattice.w[i] = 2 * ( 17 + sqr33 ) / pot4;
		
		for ( int i = 17; i < 21; i++ ) lattice.w[i] = 64 * ( 99 - 13 * sqr33 ) / ( 105 * pot4 );
		
		for ( int i = 21; i < nvel; i++ ) lattice.w[i] = 4 * ( -93 +19 * sqr33 ) / ( 105 * pot4 );
	}
	
	else cout << "\nErro!	Verificar nvel = 25" << endl;
        
    double a = 0.5 * sqrt( 0.5 * 15 - sqr33 );
    
    lattice.c_s2 = 1. / ( a * a ); // sound velocity squared
    
    lattice.one_over_c_s2 = a * a; // inverse of the sound velocity squared
}

//====================================================================================================================//




//=============================== Define sites used in the streaming process =========================================//
//
//      Input: geometry, lattice
//      Output: addresses, *ini_dir
//
//====================================================================================================================//

void def_dir_prop ( GEOMETRY geometry, LATTICE& lattice )
{
	int *ini_meio = geometry.ini; 
	int *ini_dir = lattice.ini_stream;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	bool* solid = nullptr;
	
	double *qlost_x = nullptr;
	double *qlost_y = nullptr;
	double *qlost_z = nullptr;
	
	double *ini_qlost_x = nullptr;
	double *ini_qlost_y = nullptr;
	double *ini_qlost_z = nullptr;
	
	if ( lattice.ini_mom_x != nullptr ) ini_qlost_x = lattice.ini_mom_x;
	if ( lattice.ini_mom_y != nullptr ) ini_qlost_y = lattice.ini_mom_y;
	if ( lattice.ini_mom_z != nullptr ) ini_qlost_z = lattice.ini_mom_z;

	if ( dim == 2 )
	{
		//---------------- Define os passos para propagação ------------------------------------------//

		double  step_x[nvel];
		double  step_y[nvel];
		
		int steps[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			double cx = lattice.c_i[ i * dim + 0 ];
			double cy = lattice.c_i[ i * dim + 1 ];

			int cx_int =  round_number ( cx );
			int cy_int =  round_number ( cy );

			steps[i] = abs ( cx_int );

			if ( abs ( cy_int ) > steps[i] )
			{
				steps[i] = abs ( cy_int );
			}

			step_x[i] = cx / steps[i];
			step_y[i] = cy / steps[i];
		}

		//-------------- Encontra as direções contrárias ---------------------------------------------//

		int i_op[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			int cx_i = round_number ( lattice.c_i[ i * dim + 0 ] );
			int cy_i = round_number ( lattice.c_i[ i * dim + 1 ] );

			for ( int j = 1; j < nvel; j++ )
			{
				int cx_j = round_number ( lattice.c_i[ j * dim + 0 ] );
				int cy_j = round_number ( lattice.c_i[ j * dim + 1 ] );

				if ( cx_j == -cx_i  && cy_j == -cy_i )
				{
					i_op[i] = j;
				}
			}
		}

		//--------------------------------------------------------------------------------------------//

		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int* meio_local = ini_meio + x + y * nx;

				if ( *meio_local )
				{
					if ( lattice.ini_mom_x != nullptr )
					{
						qlost_x = ini_qlost_x + ( *meio_local - 1 ) * nvel;
						qlost_x[0] = 0;
					}
					if ( lattice.ini_mom_y != nullptr )
					{
						qlost_y = ini_qlost_y + ( *meio_local - 1 ) * nvel;
						qlost_y[0] = 0;
					}
					
					if ( lattice.ini_solid != nullptr )
					{
						solid = lattice.ini_solid + ( *meio_local - 1 ) * nvel;
						solid[0] = 0;
					}					 
										
					int* dir = ini_dir + ( *meio_local - 1 ) * nvel;

					int* meio_prop = meio_local;
					
					dir[0] = ( *meio_prop - 1 ) * nvel;

					//--------------------------------------------------------------------------------//

					for ( int i = 1; i < nvel; i++ )
					{
						if ( lattice.ini_solid != nullptr ) solid[i] = 0; // fluid
						
						int inv = 0; // número de inversões

						double  stpx = step_x[i];
						double  stpy = step_y[i];

						int x_prop = x;
						int y_prop = y;

						double  x_prop_f = ( double ) x;
						double  y_prop_f = ( double ) y;

						for ( int stp = 1; stp < steps[i] + 1; stp++ )
						{
							x_prop_f = x_prop_f + stpx;
							y_prop_f = y_prop_f + stpy;

							x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
							y_prop = ( round_number ( y_prop_f ) + ny ) % ny;

							meio_prop = ini_meio + x_prop + y_prop * nx;

							if ( *meio_prop == 0 )
							{
								stpx = -stpx;
								stpy = -stpy;

								x_prop_f = x_prop_f + stpx;
								y_prop_f = y_prop_f + stpy;

								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;

								meio_prop = ini_meio + x_prop + y_prop * nx;

								inv++;
								
								if ( lattice.ini_solid != nullptr ) solid[i] = 1; // solid
							}
						}

						if ( inv % 2 == 0 )
						{
							dir[i] = i + ( *meio_prop - 1 ) * nvel;
						}
						else
						{
							dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							
							if ( lattice.ini_mom_x != nullptr )	qlost_x[i] = - 2. * stpx;
							if ( lattice.ini_mom_y != nullptr )	qlost_y[i] = - 2. * stpy;
						}
					}

					//----------------------------------------------------------------------------//
				}
			}
		}
	}
	
	else
	{
		//---------------- Define os passos para propagação ------------------------------------------//

		double  step_x[nvel];
		double  step_y[nvel];
		double  step_z[nvel];

		int steps[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			double cx = lattice.c_i[ i * dim + 0 ];
			double cy = lattice.c_i[ i * dim + 1 ];
			double cz = lattice.c_i[ i * dim + 2 ];

			int cx_int =  round_number ( cx );
			int cy_int =  round_number ( cy );
			int cz_int =  round_number ( cz );

			steps[i] = abs ( cx_int );

			if ( abs ( cy_int ) > steps[i] )
			{
				steps[i] = abs ( cy_int );
			}

			if ( abs ( cz_int ) > steps[i] )
			{
				steps[i] = abs ( cz_int );
			}

			step_x[i] = cx / steps[i];
			step_y[i] = cy / steps[i];
			step_z[i] = cz / steps[i];
		}

		//-------------- Encontra as direções contrárias ---------------------------------------------//

		int i_op[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			int cx_i = round_number ( lattice.c_i[ i * dim + 0 ] );
			int cy_i = round_number ( lattice.c_i[ i * dim + 1 ] );
			int cz_i = round_number ( lattice.c_i[ i * dim + 2 ] );

			for ( int j = 1; j < nvel; j++ )
			{
				int cx_j = round_number ( lattice.c_i[ j * dim + 0 ] );
				int cy_j = round_number ( lattice.c_i[ j * dim + 1 ] );
				int cz_j = round_number ( lattice.c_i[ j * dim + 2 ] );

				if ( cx_j == -cx_i  && cy_j == -cy_i && cz_j == -cz_i )
				{
					i_op[i] = j;
				}
			}
		}

		//--------------------------------------------------------------------------------------------//

		for ( int z = 0; z < nz; z++ )
		{
			for ( int y = 0; y < ny; y++ )
			{
				for ( int x = 0; x < nx; x++ )
				{
					int* meio_local = ini_meio + x + y * nx + z * nx * ny;

					if ( *meio_local )
					{
						if ( lattice.ini_mom_x != nullptr )
						{
							qlost_x = ini_qlost_x + ( *meio_local - 1 ) * nvel;
							qlost_x[0] = 0;
						}
						if ( lattice.ini_mom_y != nullptr )
						{
							qlost_y = ini_qlost_y + ( *meio_local - 1 ) * nvel;
							qlost_y[0] = 0;
						}
						if ( lattice.ini_mom_z != nullptr )
						{
							qlost_z = ini_qlost_z + ( *meio_local - 1 ) * nvel;
							qlost_z[0] = 0;
						}
						
						if ( lattice.ini_solid != nullptr )
						{
							solid = lattice.ini_solid + ( *meio_local - 1 ) * nvel;
							solid[0] = 0;
						}
					
						int* dir = ini_dir + ( *meio_local - 1 ) * nvel;

						int* meio_prop = meio_local;
					
						dir[0] = ( *meio_prop - 1 ) * nvel;

						//----------------------------------------------------------------------------//

						for ( int i = 1; i < nvel; i++ )
						{
							if ( lattice.ini_solid != nullptr ) solid[i] = 0; // fluid
							
							int inv = 0; // número de inversões

							double  stpx = step_x[i];
							double  stpy = step_y[i];
							double  stpz = step_z[i];

							int x_prop = x;
							int y_prop = y;
							int z_prop = z;

							double  x_prop_f = ( double ) x;
							double  y_prop_f = ( double ) y;
							double  z_prop_f = ( double ) z;

							for ( int stp = 1; stp < steps[i] + 1; stp++ )
							{
								x_prop_f = x_prop_f + stpx;
								y_prop_f = y_prop_f + stpy;
								z_prop_f = z_prop_f + stpz;

								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
								z_prop = ( round_number ( z_prop_f ) + nz ) % nz;

								meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

								if ( *meio_prop == 0 )
								{
									stpx = -stpx;
									stpy = -stpy;
									stpz = -stpz;

									x_prop_f = x_prop_f + stpx;
									y_prop_f = y_prop_f + stpy;
									z_prop_f = z_prop_f + stpz;

									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;

									meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

									inv++;
					
									if ( lattice.ini_solid != nullptr ) solid[i] = 1;
								}
							}

							if ( inv % 2 == 0 )
							{
								dir[i] = i + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
																
								if ( lattice.ini_mom_x != nullptr )	qlost_x[i] = - 2. * stpx;
								if ( lattice.ini_mom_y != nullptr )	qlost_y[i] = - 2. * stpy;
								if ( lattice.ini_mom_z != nullptr )	qlost_z[i] = - 2. * stpz;
							}
						}

						//----------------------------------------------------------------------------//
					}
				}
			}
		}
	}
}

//====================================================================================================================//




//=============================== Define sites used in the streaming process (with mirror x ) ========================//
//
//      Input: geometry, lattice
//      Output: addresses, *ini_dir
//
//====================================================================================================================//

void def_dir_prop_mirror_x ( GEOMETRY geometry, LATTICE& lattice )
{
	int *ini_meio = geometry.ini; 
	int *ini_dir = lattice.ini_stream;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	bool* solid = nullptr;
	
	//---------------- Define os passos para propagação ------------------------------------------//

	double  step_x[nvel];
	double  step_y[nvel];
	double  step_z[nvel];

	int steps[nvel];

	for ( int i = 1; i < nvel; i++ )
	{
		double cx = lattice.c_i[ i * dim + 0 ];
		double cy = lattice.c_i[ i * dim + 1 ];
		double cz = lattice.c_i[ i * dim + 2 ];

		int cx_int =  round_number ( cx );
		int cy_int =  round_number ( cy );
		int cz_int =  round_number ( cz );

		steps[i] = abs ( cx_int );

		if ( abs ( cy_int ) > steps[i] )
		{
			steps[i] = abs ( cy_int );
		}

		if ( abs ( cz_int ) > steps[i] )
		{
			steps[i] = abs ( cz_int );
		}

		step_x[i] = cx / steps[i];
		step_y[i] = cy / steps[i];
		step_z[i] = cz / steps[i];
	}

	//-------------- Encontra as direções contrárias ---------------------------------------------//

	int i_op[nvel];
	
	int i_op_mirror_x[nvel];

	for ( int i = 1; i < nvel; i++ )
	{
		int cx_i = round_number ( lattice.c_i[ i * dim + 0 ] );
		int cy_i = round_number ( lattice.c_i[ i * dim + 1 ] );
		int cz_i = round_number ( lattice.c_i[ i * dim + 2 ] );

		for ( int j = 1; j < nvel; j++ )
		{
			int cx_j = round_number ( lattice.c_i[ j * dim + 0 ] );
			int cy_j = round_number ( lattice.c_i[ j * dim + 1 ] );
			int cz_j = round_number ( lattice.c_i[ j * dim + 2 ] );

			if ( cx_j == -cx_i && cy_j == -cy_i && cz_j == -cz_i ) i_op[i] = j;
			
			if ( cx_j == -cx_i && cy_j == cy_i && cz_j == cz_i )
			{
				if ( cx_i > 0 ) i_op_mirror_x[i] = j;
				
				else i_op_mirror_x[i] = i;
			}
		}
	}
	
	//--------------------------------------------------------------------------------------------//

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int* meio_local = ini_meio + x + y * nx + z * nx * ny;

				if ( *meio_local )
				{					
					if ( lattice.ini_solid != nullptr )
					{
						solid = lattice.ini_solid + ( *meio_local - 1 ) * nvel;
						solid[0] = 0;
					}
				
					int* dir = ini_dir + ( *meio_local - 1 ) * nvel;

					int* meio_prop = meio_local;
				
					dir[0] = ( *meio_prop - 1 ) * nvel;

					//----------------------------------------------------------------------------//

					for ( int i = 1; i < nvel; i++ )
					{
						if ( lattice.ini_solid != nullptr ) solid[i] = 0; // fluid
						
						int inv = 0; // número de inversões

						double  stpx = step_x[i];
						double  stpy = step_y[i];
						double  stpz = step_z[i];

						int x_prop = x;
						int y_prop = y;
						int z_prop = z;

						double  x_prop_f = ( double ) x;
						double  y_prop_f = ( double ) y;
						double  z_prop_f = ( double ) z;

						for ( int stp = 1; stp < steps[i] + 1; stp++ )
						{
							x_prop_f = x_prop_f + stpx;
							y_prop_f = y_prop_f + stpy;
							z_prop_f = z_prop_f + stpz;
							
							if ( x != nx - 1 )
							{
								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
								z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
							}
							
							else
							{
								if ( stpx > 0 )
								{
									x_prop = nx - 1;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
								else
								{
									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
							}

							meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

							if ( *meio_prop == 0 && x != nx - 1 )
							{
								stpx = -stpx;
								stpy = -stpy;
								stpz = -stpz;

								x_prop_f = x_prop_f + stpx;
								y_prop_f = y_prop_f + stpy;
								z_prop_f = z_prop_f + stpz;

								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
								z_prop = ( round_number ( z_prop_f ) + nz ) % nz;

								meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

								inv++;
				
								if ( lattice.ini_solid != nullptr ) solid[i] = 1;
							}
							
							if ( *meio_prop == 0 && x == nx - 1 )
							{
								meio_prop = ini_meio + x + y * nx + z * nx * ny;

								inv++;
				
								if ( lattice.ini_solid != nullptr ) solid[i] = 1;
							}
						}
						
						if ( x != nx - 1 )
						{
							if ( inv % 2 == 0 )
							{
								dir[i] = i + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							}
						}
						
						else
						{
							if ( inv % 2 == 0 )
							{
								dir[i] = i_op_mirror_x[i] + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							}
						}						
					}

					//----------------------------------------------------------------------------//
				}
			}
		}
	}
}

//====================================================================================================================//




//=============================== Define sites used in the streaming process (with mirror xy ) ========================//
//
//      Input: geometry, lattice
//      Output: addresses, *ini_dir
//
//====================================================================================================================//

void def_dir_prop_mirror_xy ( GEOMETRY geometry, LATTICE& lattice )
{
	int *ini_meio = geometry.ini; 
	int *ini_dir = lattice.ini_stream;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	bool* solid = nullptr;
	
	//---------------- Define os passos para propagação ------------------------------------------//

	double  step_x[nvel];
	double  step_y[nvel];
	double  step_z[nvel];

	int steps[nvel];

	for ( int i = 1; i < nvel; i++ )
	{
		double cx = lattice.c_i[ i * dim + 0 ];
		double cy = lattice.c_i[ i * dim + 1 ];
		double cz = lattice.c_i[ i * dim + 2 ];

		int cx_int =  round_number ( cx );
		int cy_int =  round_number ( cy );
		int cz_int =  round_number ( cz );

		steps[i] = abs ( cx_int );

		if ( abs ( cy_int ) > steps[i] )
		{
			steps[i] = abs ( cy_int );
		}

		if ( abs ( cz_int ) > steps[i] )
		{
			steps[i] = abs ( cz_int );
		}

		step_x[i] = cx / steps[i];
		step_y[i] = cy / steps[i];
		step_z[i] = cz / steps[i];
	}

	//-------------- Encontra as direções contrárias ---------------------------------------------//

	int i_op[nvel];
	
	int i_op_mirror_x[nvel];
	int i_op_mirror_y[nvel];
	
	int i_op_mirror_xy[nvel];
	
	for ( int i = 1; i < nvel; i++ ) i_op_mirror_xy[i] = -1;

	for ( int i = 1; i < nvel; i++ )
	{
		int cx_i = round_number ( lattice.c_i[ i * dim + 0 ] );
		int cy_i = round_number ( lattice.c_i[ i * dim + 1 ] );
		int cz_i = round_number ( lattice.c_i[ i * dim + 2 ] );

		for ( int j = 1; j < nvel; j++ )
		{
			int cx_j = round_number ( lattice.c_i[ j * dim + 0 ] );
			int cy_j = round_number ( lattice.c_i[ j * dim + 1 ] );
			int cz_j = round_number ( lattice.c_i[ j * dim + 2 ] );

			if ( cx_j == -cx_i && cy_j == -cy_i && cz_j == -cz_i ) i_op[i] = j;
			
			if ( cx_j == -cx_i && cy_j == cy_i && cz_j == cz_i )
			{
				if ( cx_i > 0 ) i_op_mirror_x[i] = j;
				
				else i_op_mirror_x[i] = i;
			}
			
			if ( cx_j == cx_i && cy_j == -cy_i && cz_j == cz_i )
			{
				if ( cy_i > 0 ) i_op_mirror_y[i] = j;
				
				else i_op_mirror_y[i] = i;
			}
					
			if ( cx_j == -cx_i && cy_j == -cy_i && cz_j == cz_i &&  cx_i > 0 &&  cy_i > 0 ) i_op_mirror_xy[i] = j;

		}
	}
	
	for ( int i = 1; i < nvel; i++ )
	{
		if ( i_op_mirror_xy[i] == -1 )
		{
			i_op_mirror_xy[i] = i;
		
			if ( i_op_mirror_x[i] != i && i_op_mirror_y[i] == i ) i_op_mirror_xy[i] = i_op_mirror_x[i];
		
			if ( i_op_mirror_x[i] == i && i_op_mirror_y[i] != i ) i_op_mirror_xy[i] = i_op_mirror_y[i];
		}
		
	}

	//--------------------------------------------------------------------------------------------//

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int* meio_local = ini_meio + x + y * nx + z * nx * ny;

				if ( *meio_local )
				{					
					if ( lattice.ini_solid != nullptr )
					{
						solid = lattice.ini_solid + ( *meio_local - 1 ) * nvel;
						solid[0] = 0;
					}
				
					int* dir = ini_dir + ( *meio_local - 1 ) * nvel;

					int* meio_prop = meio_local;
				
					dir[0] = ( *meio_prop - 1 ) * nvel;

					//----------------------------------------------------------------------------//

					for ( int i = 1; i < nvel; i++ )
					{
						if ( lattice.ini_solid != nullptr ) solid[i] = 0; // fluid
						
						int inv = 0; // número de inversões

						double  stpx = step_x[i];
						double  stpy = step_y[i];
						double  stpz = step_z[i];

						int x_prop = x;
						int y_prop = y;
						int z_prop = z;

						double  x_prop_f = ( double ) x;
						double  y_prop_f = ( double ) y;
						double  z_prop_f = ( double ) z;

						for ( int stp = 1; stp < steps[i] + 1; stp++ )
						{
							x_prop_f = x_prop_f + stpx;
							y_prop_f = y_prop_f + stpy;
							z_prop_f = z_prop_f + stpz;
							
							if ( x != nx - 1 && y != ny - 1 )
							{
								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
								z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
							}
							
							if ( x == nx - 1 && y != ny - 1 )
							{
								if ( stpx > 0 )
								{
									x_prop = nx - 1;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
								else
								{
									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
							}
							
							if ( x != nx - 1 && y == ny - 1 )
							{
								if ( stpy > 0 )
								{
									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ny - 1;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
								else
								{
									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
							}
							
							if ( x == nx - 1 && y == ny - 1 )
							{
								if ( stpx > 0 && stpy > 0 )
								{
									x_prop = nx - 1;
									y_prop = ny - 1;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
								else if ( stpx > 0 )
								{
									x_prop = nx - 1;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
								else if ( stpy > 0 )
								{
									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ny - 1;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;
								}
							}
							
							meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

							if ( *meio_prop == 0 &&  x != nx - 1 &&  y != ny - 1 )
							{
								stpx = -stpx;
								stpy = -stpy;
								stpz = -stpz;

								x_prop_f = x_prop_f + stpx;
								y_prop_f = y_prop_f + stpy;
								z_prop_f = z_prop_f + stpz;

								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
								z_prop = ( round_number ( z_prop_f ) + nz ) % nz;

								meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

								inv++;
				
								if ( lattice.ini_solid != nullptr ) solid[i] = 1;
							}
							
							else if ( *meio_prop == 0  )
							{
								meio_prop = ini_meio + x + y * nx + z * nx * ny;

								inv++;
				
								if ( lattice.ini_solid != nullptr ) solid[i] = 1;
							}
						}
						
						if ( x != nx - 1 && y != ny - 1 )
						{
							if ( inv % 2 == 0 )
							{
								dir[i] = i + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							}
						}
						
						if ( x == nx - 1  && y != ny - 1  )
						{
							if ( inv % 2 == 0 )
							{
								dir[i] = i_op_mirror_x[i] + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							}
						}
						
						if ( x != nx - 1  && y == ny - 1  )
						{
							if ( inv % 2 == 0 )
							{
								dir[i] = i_op_mirror_y[i] + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							}
						}
						
						if ( x == nx - 1  && y == ny - 1  )
						{
							if ( inv % 2 == 0 )
							{
								dir[i] = i_op_mirror_xy[i] + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							}
						}							
					}

					//----------------------------------------------------------------------------//
				}
			}
		}
	}
}

//====================================================================================================================//



//=============================== Return the dot product =============================================================//
//
//      Input: two vetors
//      Output: dot product
//
//====================================================================================================================//

#pragma acc routine seq
double  dot_product ( double  vet_1[dim], double  vet_2[dim] )
{

    double  result = 0.0;

    for ( int i = 0; i < dim; i++ )
    {
        result = result + vet_1[i] * vet_2[i];
    }

    return result;
}

//====================================================================================================================//




//=============================== Return the dot product (of tensors) ================================================//
//
//      Input: two tensors (second order)
//      Output: dot product
//
//====================================================================================================================//

#pragma acc routine seq
double  dot_dot_product ( double* tens_A, double* tens_B )
{

    double  result = 0.0;

    for ( int i = 0; i < dim; i++ )
    {
		for ( int j = 0; j < dim; j++ )
		{
			result = result + tens_A[ i + j * dim ] * tens_B[ i + j * dim ];
		}
	}

    return result;
}

//====================================================================================================================//




//=============================== Return the outer product ===========================================================//
//
//      Input: two vetors
//      Output: outer product
//
//====================================================================================================================//

#pragma acc routine seq
void  outer_product ( double vet_1[dim], double vet_2[dim], double* result )
{

    for ( int i = 0; i < dim; i++ )
    {
		for ( int j = 0; j < dim; j++ )
		{
			result[ i + j * dim ] = vet_1[i] * vet_2[j];
		}
    }

}

//====================================================================================================================//



//=============================== Return the equilibrium distribution function =======================================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq ( double* feq, double vx, double vy, double vz, double rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;

		double vquad = ( vx * vx + vy * vy );
		
		double c_i[dim];
		
		double one_over_c_s2 = lattice.one_over_c_s2;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2
									+ 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
									- 0.5 * vquad * one_over_c_s2 );
		}
	}
	
	else
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;

		double vquad = ( vx * vx + vy * vy + vz * vz );
		
		double c_i[dim];
		
		double one_over_c_s2 = lattice.one_over_c_s2;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];        
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2
									+ 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
									- 0.5 * vquad * one_over_c_s2 );
		}
	}
}

//====================================================================================================================//




//=============================== Return the equilibrium distribution function (heat conduction) =====================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq ( double* feq, double T, LATTICE lattice )
{
	
	for ( int i = 0; i < nvel; i++ )
	{
		
		feq[i] = lattice.w[i] * T ;
	}
	
}

//====================================================================================================================//



//=============================== Return the equilibrium distribution function for solving Poisson's Equation=========//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

void dist_eq_Poisson ( double* feq, double vx, double vy, double vz, double delta_rho, double rho_avg, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;

		double vquad = ( vx * vx + vy * vy );
		
		double c_i[dim];
		
		double one_over_c_s2 = lattice.one_over_c_s2;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * ( delta_rho + cv * one_over_c_s2 * rho_avg
									+ 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2 * rho_avg
									- 0.5 * vquad * one_over_c_s2 * rho_avg );
		}
	}
	
	else
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;

		double vquad = ( vx * vx + vy * vy + vz * vz );
		
		double c_i[dim];
		
		double one_over_c_s2 = lattice.one_over_c_s2;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];        
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * ( delta_rho + cv * one_over_c_s2 * rho_avg
									+ 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2 * rho_avg
									- 0.5 * vquad * one_over_c_s2 * rho_avg );
		}
	}
}

//====================================================================================================================//




//=============================== Return the equilibrium distribution function =======================================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq_fourth ( double* feq, double vx, double vy, double vz, double rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		
		double c_i[dim];
		
		double u_u = ( vx * vx + vy * vy );
		
		double one_over_c_s2 = lattice.one_over_c_s2;    
		
		double one_over_c_s4 = one_over_c_s2 * one_over_c_s2;
		
		double one_over_c_s6 = one_over_c_s4 * one_over_c_s2;
		
		double one_over_c_s8 = one_over_c_s4 * one_over_c_s4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2 + 0.5 * cv*cv * one_over_c_s4
							  - 0.5 * u_u * one_over_c_s2 + (1./6.)*cv*cv*cv*one_over_c_s6
							  - 0.5 * cv * u_u * one_over_c_s4 + (1./24.)*cv*cv*cv*cv * one_over_c_s8 
							  - 0.25 * cv*cv*u_u * one_over_c_s6 + (1./8.)*u_u*u_u * one_over_c_s4 );
		}
	}
	
	else
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double c_i[dim];
		
		double u_u = ( vx * vx + vy * vy + vz * vz );
		
		double one_over_c_s2 = lattice.one_over_c_s2;    
		
		double one_over_c_s4 = one_over_c_s2 * one_over_c_s2;
		
		double one_over_c_s6 = one_over_c_s4 * one_over_c_s2;
		
		double one_over_c_s8 = one_over_c_s4 * one_over_c_s4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];        
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2 + 0.5 * cv*cv * one_over_c_s4
							  - 0.5 * u_u * one_over_c_s2 + (1./6.)*cv*cv*cv*one_over_c_s6
							  - 0.5 * cv * u_u * one_over_c_s4 + (1./24.)*cv*cv*cv*cv * one_over_c_s8 
							  - 0.25 * cv*cv*u_u * one_over_c_s6 + (1./8.)*u_u*u_u * one_over_c_s4 );
		}
	}
	
}

//====================================================================================================================//




//=============================== Return the equilibrium distribution function =======================================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq_sixth ( double* feq, double vx, double vy, double vz, double rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		
		double c_i[dim];
		
		double u_u = ( vx * vx + vy * vy );
		
		double one_over_c_s2 = lattice.one_over_c_s2;    
		
		double one_over_c_s4 = one_over_c_s2 * one_over_c_s2;
		
		double one_over_c_s6 = one_over_c_s4 * one_over_c_s2;
		
		double one_over_c_s8 = one_over_c_s4 * one_over_c_s4;

		double one_over_c_s10 = one_over_c_s8 * one_over_c_s2;
		
		double one_over_c_s12 = one_over_c_s8 * one_over_c_s4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2 + 0.5 * cv*cv * one_over_c_s4
							- 0.5 * u_u * one_over_c_s2 + (1./6.)*cv*cv*cv*one_over_c_s6
							- 0.5 * cv * u_u * one_over_c_s4 + (1./24.)*cv*cv*cv*cv * one_over_c_s8 
							- 0.25 * cv*cv*u_u * one_over_c_s6 + (1./8.)*u_u*u_u * one_over_c_s4 
							+ (1./120.)*cv*cv*cv*cv*cv * one_over_c_s10- (1./12.) * u_u * cv*cv*cv * one_over_c_s8 
							+ (1./8.)*u_u*u_u *cv * one_over_c_s6	+ (1./720.) * cv*cv*cv*cv*cv*cv * one_over_c_s12 
							- (1./48.) * u_u * cv*cv*cv*cv * one_over_c_s10 + (1./16.) * u_u*u_u * cv*cv * one_over_c_s8 
							- (1./48.) * u_u*u_u*u_u * one_over_c_s6 );
		}
	}
	
	else
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double c_i[dim];
		
		double u_u = ( vx * vx + vy * vy + vz * vz );
		
		double one_over_c_s2 = lattice.one_over_c_s2;    
		
		double one_over_c_s4 = one_over_c_s2 * one_over_c_s2;
		
		double one_over_c_s6 = one_over_c_s4 * one_over_c_s2;
		
		double one_over_c_s8 = one_over_c_s4 * one_over_c_s4;
		
		double one_over_c_s10 = one_over_c_s8 * one_over_c_s2;
		
		double one_over_c_s12 = one_over_c_s8 * one_over_c_s4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];        
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2 + 0.5 * cv*cv * one_over_c_s4
							- 0.5 * u_u * one_over_c_s2 + (1./6.)*cv*cv*cv*one_over_c_s6
							- 0.5 * cv * u_u * one_over_c_s4 + (1./24.)*cv*cv*cv*cv * one_over_c_s8 
							- 0.25 * cv*cv*u_u * one_over_c_s6 + (1./8.)*u_u*u_u * one_over_c_s4 
							+ (1./120.)*cv*cv*cv*cv*cv * one_over_c_s10- (1./12.) * u_u * cv*cv*cv * one_over_c_s8 
							+ (1./8.)*u_u*u_u *cv * one_over_c_s6	+ (1./720.) * cv*cv*cv*cv*cv*cv * one_over_c_s12 
							- (1./48.) * u_u * cv*cv*cv*cv * one_over_c_s10 + (1./16.) * u_u*u_u * cv*cv * one_over_c_s8 
							- (1./48.) * u_u*u_u*u_u * one_over_c_s6 );							  
		}
	}
	
}

//====================================================================================================================//



//=============================== Return the equilibrium distribution function =======================================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq_thermal ( double* feq, double vx, double vy, double vz, double rho, double theta, LATTICE lattice )
{
	if ( dim == 3 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double c_i[dim];
			
		double u_2 = ( vx * vx + vy * vy + vz * vz );
		
		double u_4 = u_2 * u_2;
		
		double theta_2 = theta * theta;
		
		const double a_2 = lattice.one_over_c_s2;    
		
		const double a_4 = a_2 * a_2;
		
		const double a_6 = a_4 * a_2;
		
		const double a_8 = a_4 * a_4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];   
			
			double ci_2 = dot_product ( c_i, c_i ); 
			
			double ci_4 = ci_2 * ci_2; 
			
			double ci_u = dot_product ( c_i, v );
			
			double ci_u2 = ci_u * ci_u;
			
			double ci_u3 = ci_u2 * ci_u;
			
			double ci_u4 = ci_u2 * ci_u2;

			feq[i] = lattice.w[i] * rho * ( 1.0 + a_2 * ci_u 
											+ 0.5 * ( a_4 * ci_u2 - a_2 * u_2 + theta *( a_2 * ci_2 - dim ) )
					+ (1./6.) * ( a_6 * ci_u3 - 3. * a_4 * u_2 * ci_u + 3. * theta * a_2 * ci_u *( a_2 * ci_2 - dim - 2. ) )
					+ (1./24.) * ( a_8 * ci_u4 - 6.* a_6 * u_2 * ci_u2 + 3. * a_4 * u_4 
					+ 3. * theta_2 * ( a_4 * ci_4 - ( dim + 2 )*( 2. * a_2 * ci_2 - dim) )
			+ 6. * theta * ( a_6 * ci_2 *ci_u2  - ( dim + 4 ) * a_4 * ci_u2 + (dim + 2) * a_2 * u_2 - a_4 * ci_2 * u_2 ) ) );
		}
	}
	//  CORRECAO: esta funcao e marcada com "acc routine seq", e um fluxo de saida nao existe do
	//  lado do dispositivo -- com -fopenacc o g++ recusa a compilacao inteira por causa desta
	//  linha. A mensagem so faz sentido no hospedeiro, entao fica condicionada a ele.

#ifndef _OPENACC
	else cout << "\nErro! Thermal Equilibrium 2D not implemented!" << endl;
#endif
}

//====================================================================================================================//



//===================== Calcula a distribuição de equilíbrio para a rede D2V37 =======================================//
//
//      Input: velocities, density, temperature, lattice vectors, distribution function
//      Output: equilibrium distribution
//
//====================================================================================================================//

void dist_eq_D2V25 ( double* feq, double vx, double vy, double rho, double theta, LATTICE lattice )
{	
	theta = 0; // Temperature
	
	if ( dim == 2 )
	{
		double v[2];
		v[0] = vx;
		v[1] = vy;
		
		double e_i[2];
		
		const double a2 = lattice.one_over_c_s2;
		const double a4 = a2 * a2;
		const double a6 = a4 * a2;
		const double a8 = a6 * a2;
		
		double theta2 = theta * theta;
		double vv = ( vx * vx + vy * vy );
		double vv2 = vv * vv;

		for ( int i = 0; i < nvel; i++ )
		{
			e_i[0] = lattice.c_i[ i * dim + 0 ];
			e_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cc = dot_product ( e_i, e_i );
			double cv = dot_product ( e_i, v );
			double cc2 = cc * cc;
			double cv2 = cv * cv;
			double cv3 = cv2 * cv;

			feq[i] = rho * lattice.w[i] * ( 1 - theta + a2 * cv + 0.5 * a4 * cv2 - 0.5 * a2 * vv 
					+ 0.5 * a2 * cc * theta + 0.5 * a4 * theta * cc * cv + ( 1./6.) * cv3 * a6
					- 2 * theta * a2 * cv - 0.5 * a4 * vv * cv )
					+ rho * lattice.w[i] * ( 0.125 * theta2 * cc2 * a4 + a2 * theta2 * cc + a2 * theta * vv
					+ 0.25 * a4 * vv * cc * theta - (3./2.) * a4 * theta * cv2 
					+ 0.25 * a6 * theta * cv2 * cc + 0.125 * a4 * vv2 - 0.25 * cv2 * vv * a6
					- (1./192) * vv2 * cc * a6 + (1./24) * a8 * vv * cc * cv2 ); 
					
		}
	}
	
}

//===================================================================================================================//



//===================== Calcula a distribuição de equilíbrio para a rede D2V37 =======================================//
//
//      Input: velocities, density, temperature, lattice vectors, distribution function
//      Output: equilibrium distribution
//
//====================================================================================================================//

void dist_eq_D2V37 ( double* feq, double vx, double vy, double rho, double theta, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];
		double c_i[dim];
		double feq2[nvel];
		double feq3[nvel];
		
		v[0] = vx;
		v[1] = vy;

		double vquad = ( vx * vx + vy * vy );
		double theta2 = theta * theta;
		double a2 = lattice.one_over_c_s2;
		double a4 = a2 * a2;	
		double a6 = a4 * a2;
		double a8 = a4 * a4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );
			double ci_ci = dot_product ( c_i, c_i );

			feq2[i] = lattice.w[i] * rho * ( 1.0 + cv * a2
									+ 0.5 * cv * cv * a4
									- 0.5 * vquad * a2 + 0.5 * theta * ( ci_ci - 2. )  );
									
			feq3[i] = feq2[i] + lattice.w[i] * rho * (1./6.)* ( a6 * cv*cv*cv - 3. * a4 * vquad * cv 
											+ 3. * a4 * theta * cv * ci_ci );
					 
			feq[i] = feq3[i] + lattice.w[i] * rho * (1./24) * ( a8 * cv * cv * cv * cv - 6. * a6 * vquad * cv * cv
											+ 3. * a4 * vquad * vquad 
											+ 3 * theta2 * ( a4 * ci_ci * ci_ci - 4 * ( 2 * a2 * ci_ci - 2 ) )
						+ 6. * theta * ( a6 * ci_ci * cv * cv - 6 * a4 * cv * cv + 4 * a2 * vquad - a4 * ci_ci * vquad ) );                        
		}
	}
}

//===================================================================================================================//



//=============================== Compute momentum ===================================================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//====================================================================================================================//

#pragma acc routine seq
void momentum ( double *f, double& mx, double& my, double& mz, LATTICE lattice )
{
	if ( dim == 2 )
	{
		mz = 0.0;
		
		mx = 0.0;
		my = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{
			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
		}
	}
	 
	else
	{
		mx = 0.0;
		my = 0.0;
		mz = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{
			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
			mz = mz + lattice.c_i[ i * dim + 2 ] * f[i];
		}
	}

}

//====================================================================================================================//



//=============================== Calculate density and velocities ===================================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//====================================================================================================================//

#pragma acc routine seq
void calcula ( double *f, double& vx, double& vy, double& vz, double& rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		vz = 0.0;
		
		double mx = 0.0;
		double my = 0.0;

		double one_over_rho;

		rho = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{

			rho = rho + f[i];

			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
		}

		if ( rho )
		{
			one_over_rho = 1.0 / ( rho );

			vx = mx * one_over_rho;
			vy = my * one_over_rho;
		}
		else
		{
			vx = 0.;
			vy = 0.;
		}
	}
	 
	else
	{
		double mx = 0.0;
		double my = 0.0;
		double mz = 0.0;

		double one_over_rho;

		rho = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{

			rho = rho + f[i];

			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
			mz = mz + lattice.c_i[ i * dim + 2 ] * f[i];
		}

		if ( rho )
		{
			one_over_rho = 1.0 / ( rho );

			vx = mx * one_over_rho;
			vy = my * one_over_rho;
			vz = mz * one_over_rho;
		}
		else
		{
			vx = 0.;
			vy = 0.;
			vz = 0.;
		}
	}
}

//====================================================================================================================//




//=============================== Returns the direction of a vector ==================================================//
//
//      Input: vector
//      Output: unit vector, *ux, *uy, *uz
//
//====================================================================================================================//

#pragma acc routine seq
void unit_vector ( double vx, double vy, double vz, double& ux, double& uy, double& uz )
{
    double modulo = ( sqrt ( vx * vx + vy * vy + vz * vz ) );

    double one_over_modulo;

    if ( modulo )
    {
        one_over_modulo = 1.0 / modulo;

        ux = vx * one_over_modulo;
        uy = vy * one_over_modulo;
        uz = vz * one_over_modulo;
    }
    else
    {
        ux = uy = uz = 0.0;
    }
}

//====================================================================================================================//



//=============================== Calculate density, velocities and temperature ======================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density, temperature
//
//====================================================================================================================//

#pragma acc routine seq
void calcula_th ( double *f, double& vx, double& vy, double& vz, double& rho, double& Tmp, LATTICE lattice )
{
	if ( dim == 2 )
	{
		vz = 0.0;
		
		double mx = 0.0;
		double my = 0.0;

		double one_over_rho;

		rho = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{

			rho = rho + f[i];

			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
		}

		if ( rho )
		{
			one_over_rho = 1.0 / ( rho );

			vx = mx * one_over_rho;
			vy = my * one_over_rho;
		}
		else
		{
			vx = 0.;
			vy = 0.;
		}
		
		double flut[dim];
		double u[dim];
	
		double rhoe = 0.;

		u[0] = vx;
		u[1] = vy;

		for (int i = 0 ;i < nvel; i++ )
		{
			flut[0] = lattice.c_i[ i * dim + 0 ] - u[0];
			flut[1] = lattice.c_i[ i * dim + 1 ] - u[1];

			rhoe = rhoe + f[i] * 0.5 * dot_product( flut, flut );
		}

		if ( rho ) Tmp = rhoe / rho;
		else Tmp = 0.0;
	}
	 
	else
	{
		double mx = 0.0;
		double my = 0.0;
		double mz = 0.0;

		double one_over_rho;

		rho = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{

			rho = rho + f[i];

			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
			mz = mz + lattice.c_i[ i * dim + 2 ] * f[i];
		}

		if ( rho )
		{
			one_over_rho = 1.0 / ( rho );

			vx = mx * one_over_rho;
			vy = my * one_over_rho;
			vz = mz * one_over_rho;
		}
		else
		{
			vx = 0.;
			vy = 0.;
			vz = 0.;
		}
		
		double flut[dim];
		double u[dim];
	
		double rhoe = 0;

		u[0] = vx;
		u[1] = vy;
		u[2] = vz;

		for (int i = 0 ;i < nvel; i++ )
		{
			flut[0] = lattice.c_i[ i * dim + 0 ] - u[0];
			flut[1] = lattice.c_i[ i * dim + 1 ] - u[1];
			flut[2] = lattice.c_i[ i * dim + 2 ] - u[2];

			rhoe = rhoe + f[i] * 0.5 * dot_product( flut, flut );
		}

		if ( rho ) Tmp = rhoe / rho;
		else Tmp = 0.0;
	}
}

//====================================================================================================================//




//=============================== Calculates the vorticity field =====================================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//====================================================================================================================//

void vorticity ( double *f, double& vort_x, double& vort_y, double& vort_z, int x, int y, int z, GEOMETRY geometry, 
					LATTICE lattice )
{
	double rho = 1.0;
	
    int x_p = ( x + 1 ) % geometry.nx;
    int x_m = ( x - 1 + geometry.nx ) % geometry.nx;

    int y_p = ( y + 1 ) % geometry.ny;
    int y_m = ( y - 1 + geometry.ny ) % geometry.ny;

    int z_p = ( z + 1 ) % geometry.nz;
    int z_m = ( z - 1 + geometry.nz ) % geometry.nz;

    int *x_mais, *x_menos, *y_mais, *y_menos, *z_mais, *z_menos;

    x_mais = geometry.ini + x_p + y * geometry.nx + z * geometry.ny * geometry.nx;
    x_menos = geometry.ini + x_m + y * geometry.nx + z * geometry.ny * geometry.nx;

    y_mais = geometry.ini + x + y_p * geometry.nx + z * geometry.ny * geometry.nx;
    y_menos = geometry.ini + x + y_m * geometry.nx + z * geometry.ny * geometry.nx;

    z_mais = geometry.ini + x + y * geometry.nx + z_p * geometry.ny * geometry.nx;
    z_menos = geometry.ini + x + y * geometry.nx + z_m * geometry.ny * geometry.nx;

    double vx_x_mais, vx_x_menos, vy_x_mais, vy_x_menos, vz_x_mais, vz_x_menos;

    if ( *x_mais )
    {
        f = lattice.inif + ( *x_mais - 1 ) * nvel;

        calcula ( f, vx_x_mais, vy_x_mais, vz_x_mais, rho, lattice );
    }
    else
    {
        vx_x_mais = vy_x_mais = vz_x_mais = 0.0;
    }

    if ( *x_menos )
    {
        f = lattice.inif + ( *x_menos - 1 ) * nvel;

        calcula ( f, vx_x_menos, vy_x_menos, vz_x_menos, rho, lattice );
    }
    else
    {
        vx_x_menos = vy_x_menos = vz_x_menos = 0.0;
    }

    double vx_y_mais, vx_y_menos, vy_y_mais, vy_y_menos, vz_y_mais, vz_y_menos;

    if ( *y_mais )
    {
        f = lattice.inif + ( *y_mais - 1 ) * nvel;

        calcula ( f, vx_y_mais, vy_y_mais, vz_y_mais, rho, lattice );
    }
    else
    {
        vx_y_mais = vy_y_mais = vz_y_mais = 0.0;
    }

    if ( *y_menos )
    {
        f = lattice.inif + ( *y_menos - 1 ) * nvel;

        calcula ( f, vx_y_menos, vy_y_menos, vz_y_menos, rho, lattice );
    }
    else
    {
        vx_y_menos = vy_y_menos = vz_y_menos = 0.0;
    }

    double vx_z_mais, vx_z_menos, vy_z_mais, vy_z_menos, vz_z_mais, vz_z_menos;

    if ( *z_mais )
    {
        f = lattice.inif + ( *z_mais - 1 ) * nvel;

        calcula ( f, vx_z_mais, vy_z_mais, vz_z_mais, rho, lattice );
    }
    else
    {
        vx_z_mais = vy_z_mais = vz_z_mais = 0.0;
    }

    if ( *z_menos )
    {
        f = lattice.inif + ( *z_menos - 1 ) * nvel;

        calcula ( f, vx_z_menos, vy_z_menos, vz_z_menos, rho, lattice );
    }
    else
    {
        vx_z_menos = vy_z_menos = vz_z_menos = 0.0;
    }

    double deriv_vz_y = ( vz_y_mais - vz_y_menos ) / 2.0;
    double deriv_vy_z = ( vy_z_mais - vy_z_menos ) / 2.0;

    vort_x = deriv_vz_y - deriv_vy_z;

    double deriv_vx_z = ( vx_z_mais - vx_z_menos ) / 2.0;
    double deriv_vz_x = ( vz_x_mais - vz_x_menos ) / 2.0;

    vort_y = deriv_vx_z - deriv_vz_x;

    double deriv_vy_x = ( vy_x_mais - vy_x_menos ) / 2.0;
    double deriv_vx_y = ( vx_y_mais - vx_y_menos ) / 2.0;

    vort_z = deriv_vy_x - deriv_vx_y;
}

//====================================================================================================================//



//=============================== Propagation step for one site ======================================================//
//
//      Input: lattice, site
//      Output: 
//
//====================================================================================================================//

#pragma acc routine seq
void propag_site ( LATTICE lattice, int pto, double &mx, double &my, double &mz )
{

	double *qlost_x = lattice.ini_mom_x;
	double *qlost_y = lattice.ini_mom_y;
	double *qlost_z = lattice.ini_mom_z;
	
	double sum_qlost_x = 0.0;
	double sum_qlost_y = 0.0;
	double sum_qlost_z = 0.0;

	if ( lattice.inif != nullptr ) // Monophasic
	{
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f = lattice.inif + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;
		
		if ( lattice.ini_mom_x != nullptr ) qlost_x = lattice.ini_mom_x + ( pto ) * nvel;
		if ( lattice.ini_mom_y != nullptr ) qlost_y = lattice.ini_mom_y + ( pto ) * nvel;
		if ( lattice.ini_mom_z != nullptr ) qlost_z = lattice.ini_mom_z + ( pto ) * nvel;

		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_new = lattice.inif_new + dir[i];

			*f_new = f[i];
			
			if ( lattice.ini_mom_x != nullptr ) sum_qlost_x = sum_qlost_x + f[i] * qlost_x[i];
			if ( lattice.ini_mom_y != nullptr ) sum_qlost_y = sum_qlost_y + f[i] * qlost_y[i];
			if ( lattice.ini_mom_z != nullptr ) sum_qlost_z = sum_qlost_z + f[i] * qlost_z[i];
		}
	}
    else 
    {
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f_R = lattice.inif_R + ( pto ) * nvel;
		double *f_B = lattice.inif_B + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;
		
		if ( lattice.ini_mom_x != nullptr ) qlost_x = lattice.ini_mom_x + ( pto ) * nvel;
		if ( lattice.ini_mom_y != nullptr ) qlost_y = lattice.ini_mom_y + ( pto ) * nvel;
		if ( lattice.ini_mom_z != nullptr ) qlost_z = lattice.ini_mom_z + ( pto ) * nvel;

		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_R_new = lattice.inif_R_new + dir[i];
			double* f_B_new = lattice.inif_B_new + dir[i];

			*f_R_new = f_R[i];
			*f_B_new = f_B[i];
			
			if ( lattice.ini_mom_x != nullptr ) sum_qlost_x = sum_qlost_x + ( f_R[i] + f_B[i] ) * qlost_x[i];
			if ( lattice.ini_mom_y != nullptr ) sum_qlost_y = sum_qlost_y + ( f_R[i] + f_B[i] ) * qlost_y[i];
			if ( lattice.ini_mom_z != nullptr ) sum_qlost_z = sum_qlost_z + ( f_R[i] + f_B[i] ) * qlost_z[i];
		}
	}
    
    mx = sum_qlost_x;
    my = sum_qlost_y;
    mz = sum_qlost_z;
	
}

//====================================================================================================================//



//=============================== Propagation step for one site (with mirror xy) =====================================//
//
//      Input: lattice, site
//      Output: 
//
//====================================================================================================================//

void propag_site_mirror ( LATTICE lattice, int pto )
{

	if ( lattice.inif != nullptr ) // Monophasic
	{
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f = lattice.inif + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;

		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_new = lattice.inif_new + dir[i];

			*f_new = f[i];
		}
	}
    else 
    {
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f_R = lattice.inif_R + ( pto ) * nvel;
		double *f_B = lattice.inif_B + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;

		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_R_new = lattice.inif_R_new + dir[i];
			double* f_B_new = lattice.inif_B_new + dir[i];

			*f_R_new = f_R[i];
			*f_B_new = f_B[i];
		}
	}

}

//====================================================================================================================//



//=============================== Emission of the mediators for one site ==============================================//
//
//      Input: mediator populations of the site, the phase carried by them, lattice
//      Output: mediator populations, ready to be propagated
//
//      Emite  M_i = w_i fase / cs^2 .  Depois da propagacao, M_i no sitio x vale  w_i fase( x - c_i ) ,
//      e entao
//
//          soma_i M_i c_i  =  ( 1 / cs^2 ) soma_i w_i fase( x - c_i ) c_i  =  - grad( fase )
//
//      porque  soma_i w_i c_i_alfa c_i_beta = cs^2 delta_alfa_beta .  O 1/cs^2 e o que faz o momento
//      dos mediadores ser o GRADIENTE, e nao cs^2 vezes ele.
//
//      ATENCAO -- correcao de 2026-08-26.  Ate esta data a emissao era  M_i = w_i fase , sem o
//      1/cs^2, e o momento dos mediadores saia cs^2 = 1/3 vezes menor que o gradiente.  Como as
//      sobrecargas dos operadores que recebem ( x, y, z, geometry ) usam gradient(), que ja devolve
//      o gradiente verdadeiro, as duas rotas do MESMO operador davam tensoes interfaciais na razao
//      de 3 para 1 com o mesmo fat_R_B.  Medido no coll_RPL_D3Q19: sigma = 0.014944 pelos
//      mediadores contra 0.044832 por gradient(), razao 3.0002.
//
//      Consequencia: fat_R_B passou a valer 3 vezes mais nos programas com mediadores.  Para
//      reproduzir uma tensao interfacial obtida antes desta data, divida o fat_R_B por 3.
//
//      O fator sai de lattice.one_over_c_s2, entao vale para qualquer rede -- D2Q9, D3Q19, D3Q27 e
//      as multivelocidade, onde cs^2 nao e 1/3.
//
//====================================================================================================================//

#pragma acc routine seq
void emite_mediadores ( double *f_m, double fase, LATTICE lattice )
{
	double fat = fase * lattice.one_over_c_s2;

	for ( int i = 0; i < nvel; i++ )   f_m[i] = lattice.w[i] * fat;
}

//====================================================================================================================//



//=============================== Propagation step for one site (mediators) ==========================================//
//
//      Input: lattice, site
//      Output: 
//
//      Nas direcoes que dao em solido, injeta a molhabilidade no lugar do vizinho.  O valor injetado
//      segue a MESMA normalizacao da emissao ( ver emite_mediadores ): a parede se comporta como um
//      fluido de fase wett_R.  Se um dos dois nao tiver o 1/cs^2, wett_R muda de significado e o
//      angulo de contato vai junto.
//
//====================================================================================================================//

#pragma acc routine seq
void propag_site_med ( LATTICE lattice, PARAMETERS parameters, int pto )
{

	if ( lattice.inif_m != nullptr ) // Mediators
	{
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f_m = lattice.inif_m + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;
		
		bool *solid = lattice.ini_solid + ( pto ) * nvel;
		
		double *W = lattice.w;
		
		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_m_new = lattice.inif_m_new + dir[i];

			if ( solid[i] )	*f_m_new = W[i] * parameters.wett_R * lattice.one_over_c_s2;  

			else *f_m_new = f_m[i];
			
		}
	}
}

//====================================================================================================================//



//=============================== Record the velocity field (monophasic) =============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_velocity ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    string nomevel = "vel_" + to_string( passo ) + ".vtk";

    ofstream fvel ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx * ny * nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
				//int x_desloc = ( x + nx / 2 ) % nx;
				//int y_desloc = ( y + 2*ny / 3 ) % ny;
				
				int x_desloc = x;
				int y_desloc = y;
				
                int *meio = geometry.ini + x_desloc + y_desloc * nx + z * ny * nx;

                if ( *meio )
                {
                    if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						double vx, vy, vz, rho;

						calcula ( f, vx, vy, vz, rho, lattice );

						fvel << vx  << " " << vy << " " << vz << " ";
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						double vx_R, vy_R, vz_R, rho_R;

						calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
						
						double vx_B, vy_B, vz_B, rho_B;

						calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );
						
						double concR = rho_R / ( rho_R + rho_B );
						double concB = 1. - concR;
						
						double vx = concR * vx_R + concB * vx_B;
						double vy = concR * vy_R + concB * vy_B;
						double vz = concR * vz_R + concB * vz_B;

						//  Correcao da meia-forca de Guo, se o programa tiver gravado a forca:
						//  pos-colisao,  rho u = soma_i f_i c_i - F / 2 .  Sem isto o campo
						//  gravado e dominado por |F|/2rho na interface, que nao e velocidade.

						if ( lattice.ini_force != nullptr )
						{
							double rho = rho_R + rho_B;

							const double *F = lattice.ini_force + ( long long ) ( *meio - 1 ) * 3;

							if ( rho > 0.0 )
							{
								vx = vx - 0.5 * F[0] / rho;
								vy = vy - 0.5 * F[1] / rho;
								vz = vz - 0.5 * F[2] / rho;
							}
						}

						fvel << vx  << " " << vy << " " << vz << " ";
					}
                }
                else
                {
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}

//====================================================================================================================//



//=============================== Record the velocity field (monophasic) =============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_velocity ( GEOMETRY geometry, LATTICE lattice, double tempo )
{
    string nomevel = "vel_" + to_string( tempo ) + ".vtk";

    ofstream fvel ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio )
                {
					if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						double vx, vy, vz, rho;

						calcula ( f, vx, vy, vz, rho, lattice );

						fvel << vx  << " " << vy << " " << vz << " ";
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						double vx_R, vy_R, vz_R, rho_R;

						calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
						
						double vx_B, vy_B, vz_B, rho_B;

						calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );
						
						double concR = rho_R / ( rho_R + rho_B );
						double concB = 1. - concR;
						
						double vx = concR * vx_R + concB * vx_B;
						double vy = concR * vy_R + concB * vy_B;
						double vz = concR * vz_R + concB * vz_B;

						fvel << vx  << " " << vy << " " << vz << " ";
					}
					
                }
                else
                {
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}

//====================================================================================================================//




//=============================== Record the vorticity field (monophasic) ============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_vorticity ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
    string name_vor = "vor_" + to_string( passo ) + ".vtk";

    ofstream fvor ( name_vor );

    fvor << "# vtk DataFile Version 2.0" << endl;
    fvor << "Velocidade" << endl;
    fvor << "ASCII" << endl;
    fvor << "DATASET STRUCTURED_POINTS" << endl;
    fvor << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    fvor << "ASPECT_RATIO 1 1 1" << endl;
    fvor << "ORIGIN 0 0 0" << endl;
    fvor << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    fvor << "VECTORS vorticity double" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio )
                {
                    double *f = lattice.inif + ( *meio - 1 ) * nvel;

                    double vort_x, vort_y, vort_z;

                    vorticity ( f, vort_x, vort_y, vort_z, x, y, z, geometry, lattice );

                    fvor << vort_x  << " " << vort_y << " " << vort_z << " ";
                }
                else
                {
                    fvor << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvor << endl;
        }

    }
    fvor.close();
}

//====================================================================================================================//



//=============================== Record the vorticity field (monophasic) ============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_vorticity ( GEOMETRY geometry, LATTICE lattice, double tempo )
{
    string name_vor = "vor_" + to_string( tempo ) + ".vtk";

    ofstream fvor ( name_vor );

    fvor << "# vtk DataFile Version 2.0" << endl;
    fvor << "Velocidade" << endl;
    fvor << "ASCII" << endl;
    fvor << "DATASET STRUCTURED_POINTS" << endl;
    fvor << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    fvor << "ASPECT_RATIO 1 1 1" << endl;
    fvor << "ORIGIN 0 0 0" << endl;
    fvor << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    fvor << "VECTORS vorticity double" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio )
                {
                    double *f = lattice.inif + ( *meio - 1 ) * nvel;

                    double vort_x, vort_y, vort_z;

                    vorticity ( f, vort_x, vort_y, vort_z, x, y, z, geometry, lattice );

                    fvor << vort_x  << " " << vort_y << " " << vort_z << " ";
                }
                else
                {
                    fvor << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvor << endl;
        }

    }
    fvor.close();
}

//====================================================================================================================//



//=============================== Record the density field (monophasic) ==============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_density ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
    string name_rho = "rho_" + to_string( passo ) + ".vtk";

    ofstream frho ( name_rho );

    frho << "# vtk DataFile Version 2.0" << endl;
    frho << "Densidade" << endl;
    frho << "ASCII" << endl;
    frho << "DATASET STRUCTURED_POINTS" << endl;
    frho << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    frho << "ASPECT_RATIO 1 1 1" << endl;
    frho << "ORIGIN 0 0 0" << endl;
    frho << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    frho << "SCALARS densidade double" << endl;
    frho << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio )
                {
					if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						double rho = density( f );

						frho << rho << " ";
					}
					else 
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						double rho_R = density( f_R );
						double rho_B = density( f_B );
						
						double rho = rho_R + rho_B;
						
						frho << rho << " ";
					}
                }
                else
                {
                    frho << 0.0 << " ";
                }
            }

            frho << endl;
        }

    }
    frho.close();
}
//====================================================================================================================//




//=============================== Record the density field (one of two fluids) =======================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_density ( string name, GEOMETRY geometry, double* ini_f, unsigned int passo )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    string name_rho = name + to_string( passo ) + ".vtk";

    ofstream frho ( name_rho );

    frho << "# vtk DataFile Version 2.0" << endl;
    frho << "Densidade" << endl;
    frho << "ASCII" << endl;
    frho << "DATASET STRUCTURED_POINTS" << endl;
    frho << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    frho << "ASPECT_RATIO 1 1 1" << endl;
    frho << "ORIGIN 0 0 0" << endl;
    frho << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    frho << "SCALARS densidade double" << endl;
    frho << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
				//int x_desloc = ( x + nx / 2 ) % nx;
				//int y_desloc = ( y + 2*ny / 3 ) % ny;
				
				int x_desloc = x;
				int y_desloc = y;
				
                int *meio = geometry.ini + x_desloc + y_desloc * nx + z * ny * nx;

                if ( *meio )
                {
					double *f = ini_f + ( *meio - 1 ) * nvel;

					double rho = density( f );

					frho << rho << " ";					
                }
                else
                {
                    frho << 0.0 << " ";
                }
            }

            frho << endl;
        }

    }
    frho.close();
}
//====================================================================================================================//



//=============================== Record the phase (rho_R - rho_B) ===================================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_phase ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
    string name_phase = "phase_" + to_string( passo ) + ".vtk";

    ofstream file_phase ( name_phase );

    file_phase << "# vtk DataFile Version 2.0" << endl;
    file_phase << "Densidade" << endl;
    file_phase << "ASCII" << endl;
    file_phase << "DATASET STRUCTURED_POINTS" << endl;
    file_phase << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    file_phase << "ASPECT_RATIO 1 1 1" << endl;
    file_phase << "ORIGIN 0 0 0" << endl;
    file_phase << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    file_phase << "SCALARS densidade double" << endl;
    file_phase << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio )
                {
                    double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
                    double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

                    double rho_R = density( f_R );
                    double rho_B = density( f_B );
                    
                    double phase = ( rho_R - rho_B ) / ( rho_R + rho_B );

                    file_phase << phase << "	";
                }
                else
                {
                    file_phase << 0.0 << "	";
                }
            }

            file_phase << endl;
        }

    }
    file_phase.close();
}
//====================================================================================================================//



//=============================== Record the geometry of the phases (rho_R - rho_B) ==================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_phase_geometry ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
    string name_phase = "geometry_phase_" + to_string( passo ) + ".vtk";

    ofstream file_phase ( name_phase );

    file_phase << "# vtk DataFile Version 2.0" << endl;
    file_phase << "Densidade" << endl;
    file_phase << "ASCII" << endl;
    file_phase << "DATASET STRUCTURED_POINTS" << endl;
    file_phase << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    file_phase << "ASPECT_RATIO 1 1 1" << endl;
    file_phase << "ORIGIN 0 0 0" << endl;
    file_phase << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    file_phase << "SCALARS densidade int" << endl;
    file_phase << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio )
                {
                    double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
                    double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

                    double rho_R = density( f_R );
                    double rho_B = density( f_B );
                    
                    double phase = ( rho_R - rho_B ) / ( rho_R + rho_B );

                    if ( phase > 0 ) file_phase << 2 << "	";
                    
                    else file_phase << 1 << "	";
                }
                else
                {
                    file_phase << 0 << "	";
                }
            }

            file_phase << endl;
        }

    }
    file_phase.close();
}
//====================================================================================================================//




//=============================== Record the density field (monophasic) ==============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_geo ( GEOMETRY geometry, string name_file )
{
    ofstream file_geo ( name_file );

    file_geo << "# vtk DataFile Version 2.0" << endl;
    file_geo << "Densidade" << endl;
    file_geo << "ASCII" << endl;
    file_geo << "DATASET STRUCTURED_POINTS" << endl;
    file_geo << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    file_geo << "ASPECT_RATIO 1 1 1" << endl;
    file_geo << "ORIGIN 0 0 0" << endl;
    file_geo << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    file_geo << "SCALARS densidade double" << endl;
    file_geo << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < geometry.nz; z++ )
    {
        for ( int y = 0; y < geometry.ny; y++ )
        {
            for ( int x = 0; x < geometry.nx; x++ )
            {
                int *meio = geometry.ini + x + y * geometry.nx + z * geometry.ny * geometry.nx;

                if ( *meio ) file_geo << 1 << " ";
                
                else file_geo  << 0.0 << " ";
            }

            file_geo << endl;
        }

    }
    file_geo.close();
}
//====================================================================================================================//




//=============================== Record the density a scalar field ==================================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_scalar_field ( GEOMETRY geometry, string name_file, double* scalar )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    ofstream file_scalar ( name_file );

    file_scalar << "# vtk DataFile Version 2.0" << endl;
    file_scalar << "Densidade" << endl;
    file_scalar << "ASCII" << endl;
    file_scalar << "DATASET STRUCTURED_POINTS" << endl;
    file_scalar << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    file_scalar << "ASPECT_RATIO 1 1 1" << endl;
    file_scalar << "ORIGIN 0 0 0" << endl;
    file_scalar << "POINT_DATA " << nx * ny * nz << endl;
    file_scalar << "SCALARS densidade double" << endl;
    file_scalar << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pto = x + y * nx + z * nx * ny;
                
				file_scalar << scalar[pto] << "	";
            }

            file_scalar << endl;
        }

    }
    file_scalar.close();
}
//====================================================================================================================//



//====================================================================================================================//
//                                                                                                                    //
//   G R A V A C A O   E M   V T K   B I N A R I O                                                                    //
//                                                                                                                    //
//   As funcoes rec_*_bin abaixo escrevem exatamente os mesmos campos das rec_* em ASCII, no formato                   //
//   VTK legacy binario.  Regras do formato (VTK File Formats, legacy):                                                //
//                                                                                                                    //
//      * o cabecalho continua sendo texto, uma linha por item, terminada em '\n';                                     //
//      * a terceira linha e BINARY em vez de ASCII;                                                                   //
//      * logo depois de "LOOKUP_TABLE default" (escalares) ou de "VECTORS ..." (vetores) vem o bloco                  //
//        de dados cru, sem separadores, seguido de um '\n';                                                           //
//      * os dados sao gravados em BIG-ENDIAN, independentemente da maquina;                                           //
//      * o tipo declarado (float / double / int) tem de bater com o numero de bytes gravados.                         //
//                                                                                                                    //
//   Precisao: por omissao grava-se float (4 bytes).  Passando double_prec = true grava-se double                      //
//   (8 bytes).  Para visualizar, float e mais do que suficiente -- o campo tem 7 digitos significativos,              //
//   enquanto o ASCII gravava 6.  Ou seja, o arquivo binario em float e ao mesmo tempo menor e um pouco                //
//   mais preciso do que o ASCII que substitui.                                                                        //
//                                                                                                                    //
//   Tamanho, por voxel:                                                                                               //
//                                                                                                                    //
//      escalar ASCII   ~ 9 a 13 bytes        escalar binario float   4 bytes                                          //
//      vetor  ASCII    ~ 27 a 39 bytes       vetor  binario float   12 bytes                                          //
//                                                                                                                    //
//   Os arquivos sao lidos pelo ParaView, VisIt, meshio e pelo vtk do Python sem nenhuma mudanca --                     //
//   o leitor reconhece o formato pela linha BINARY.                                                                    //
//                                                                                                                    //
//====================================================================================================================//


//------------------ A maquina e little-endian? ---------------------------------------------------------------------//

inline bool vtk_host_little_endian ()
{
	const unsigned int um = 1u;

	return *( const unsigned char * ) &um == 1u;
}

//------------------ Empilha um valor no buffer, em big-endian ------------------------------------------------------//

template < class T >
inline void vtk_push ( vector<char> &buf, T valor )
{
	const char *p = ( const char * ) &valor;

	if ( vtk_host_little_endian() )
	{
		for ( int i = ( int ) sizeof ( T ) - 1; i >= 0; i-- ) buf.push_back ( p[ i ] );
	}
	else
	{
		for ( size_t i = 0; i < sizeof ( T ); i++ ) buf.push_back ( p[ i ] );
	}
}

//------------------ Empilha um valor real, na precisao escolhida ---------------------------------------------------//

inline void vtk_push_real ( vector<char> &buf, double valor, bool double_prec )
{
	if ( double_prec ) vtk_push< double > ( buf, valor );

	else               vtk_push< float  > ( buf, ( float ) valor );
}

//------------------ Cabecalho comum --------------------------------------------------------------------------------//

inline void vtk_bin_header ( ofstream &f, const string &titulo, int nx, int ny, int nz )
{
	//  Nada de 'endl' aqui: o formato pede '\n', e endl ainda forcaria um flush por linha.

	f << "# vtk DataFile Version 2.0\n";
	f << titulo << "\n";
	f << "BINARY\n";
	f << "DATASET STRUCTURED_POINTS\n";
	f << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
	f << "ASPECT_RATIO 1 1 1\n";
	f << "ORIGIN 0 0 0\n";
	f << "POINT_DATA " << ( long long ) nx * ny * nz << "\n";
}

inline void vtk_bin_scalars ( ofstream &f, const string &nome, bool double_prec )
{
	f << "SCALARS " << nome << ( double_prec ? " double" : " float" ) << "\n";
	f << "LOOKUP_TABLE default\n";
}

inline void vtk_bin_scalars_int ( ofstream &f, const string &nome )
{
	f << "SCALARS " << nome << " int\n";
	f << "LOOKUP_TABLE default\n";
}

//  Mascaras ( geometria, rotulo de fase ) cabem num byte por voxel.  E o unico caso em que o ASCII
//  chegava a ser competitivo -- gravava "0 " ou "1 ", dois caracteres -- e em que float ou int
//  dariam um arquivo maior do que ele.

inline void vtk_bin_scalars_byte ( ofstream &f, const string &nome )
{
	f << "SCALARS " << nome << " unsigned_char\n";
	f << "LOOKUP_TABLE default\n";
}

inline void vtk_bin_vectors ( ofstream &f, const string &nome, bool double_prec )
{
	f << "VECTORS " << nome << ( double_prec ? " double" : " float" ) << "\n";
}

//------------------ Descarrega o buffer de um plano z --------------------------------------------------------------//
//
//  O buffer guarda um plano z de cada vez: para uma geometria de 300^3 sao 300*300*3*4 = 1.08 MB,
//  em vez dos 324 MB que custaria montar o campo inteiro na memoria.

inline void vtk_flush ( ofstream &f, vector<char> &buf )
{
	if ( ! buf.empty() )
	{
		f.write ( buf.data(), ( streamsize ) buf.size() );

		buf.clear();
	}
}

//====================================================================================================================//



//=============================== Record the velocity field (binary) =================================================//
//
//      Input: geometry, lattice, step, precisao
//      Output: vel_<passo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_velocity_bin ( GEOMETRY geometry, LATTICE lattice, unsigned int passo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string nomevel = "vel_" + to_string( passo ) + ".vtk";

	ofstream fvel ( nomevel, ios::binary );

	vtk_bin_header  ( fvel, "Velocidade", nx, ny, nz );
	vtk_bin_vectors ( fvel, "velocidade", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * 3 * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double vx = 0.0, vy = 0.0, vz = 0.0;

				if ( *meio )
				{
					if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						double rho;

						calcula ( f, vx, vy, vz, rho, lattice );
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						double vx_R, vy_R, vz_R, rho_R;
						double vx_B, vy_B, vz_B, rho_B;

						calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
						calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

						double concR = rho_R / ( rho_R + rho_B );
						double concB = 1. - concR;

						vx = concR * vx_R + concB * vx_B;
						vy = concR * vy_R + concB * vy_B;
						vz = concR * vz_R + concB * vz_B;
					}
				}

				vtk_push_real ( buf, vx, double_prec );
				vtk_push_real ( buf, vy, double_prec );
				vtk_push_real ( buf, vz, double_prec );
			}
		}

		vtk_flush ( fvel, buf );
	}

	fvel << "\n";

	fvel.close();
}

//====================================================================================================================//



//=============================== Record the velocity field (binary, tempo) ==========================================//
//
//      Input: geometry, lattice, tempo, precisao
//      Output: vel_<tempo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_velocity_bin ( GEOMETRY geometry, LATTICE lattice, double tempo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string nomevel = "vel_" + to_string( tempo ) + ".vtk";

	ofstream fvel ( nomevel, ios::binary );

	vtk_bin_header  ( fvel, "Velocidade", nx, ny, nz );
	vtk_bin_vectors ( fvel, "velocidade", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * 3 * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double vx = 0.0, vy = 0.0, vz = 0.0;

				if ( *meio )
				{
					if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						double rho;

						calcula ( f, vx, vy, vz, rho, lattice );
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						double vx_R, vy_R, vz_R, rho_R;
						double vx_B, vy_B, vz_B, rho_B;

						calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
						calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );

						double concR = rho_R / ( rho_R + rho_B );
						double concB = 1. - concR;

						vx = concR * vx_R + concB * vx_B;
						vy = concR * vy_R + concB * vy_B;
						vz = concR * vz_R + concB * vz_B;
					}
				}

				vtk_push_real ( buf, vx, double_prec );
				vtk_push_real ( buf, vy, double_prec );
				vtk_push_real ( buf, vz, double_prec );
			}
		}

		vtk_flush ( fvel, buf );
	}

	fvel << "\n";

	fvel.close();
}

//====================================================================================================================//



//=============================== Record the vorticity field (binary) ================================================//
//
//      Input: geometry, lattice, step, precisao
//      Output: vor_<passo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_vorticity_bin ( GEOMETRY geometry, LATTICE lattice, unsigned int passo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string name_vor = "vor_" + to_string( passo ) + ".vtk";

	ofstream fvor ( name_vor, ios::binary );

	vtk_bin_header  ( fvor, "Vorticidade", nx, ny, nz );
	vtk_bin_vectors ( fvor, "vorticity", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * 3 * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double vort_x = 0.0, vort_y = 0.0, vort_z = 0.0;

				if ( *meio )
				{
					double *f = lattice.inif + ( *meio - 1 ) * nvel;

					vorticity ( f, vort_x, vort_y, vort_z, x, y, z, geometry, lattice );
				}

				vtk_push_real ( buf, vort_x, double_prec );
				vtk_push_real ( buf, vort_y, double_prec );
				vtk_push_real ( buf, vort_z, double_prec );
			}
		}

		vtk_flush ( fvor, buf );
	}

	fvor << "\n";

	fvor.close();
}

//====================================================================================================================//



//=============================== Record the vorticity field (binary, tempo) =========================================//
//
//      Input: geometry, lattice, tempo, precisao
//      Output: vor_<tempo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_vorticity_bin ( GEOMETRY geometry, LATTICE lattice, double tempo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string name_vor = "vor_" + to_string( tempo ) + ".vtk";

	ofstream fvor ( name_vor, ios::binary );

	vtk_bin_header  ( fvor, "Vorticidade", nx, ny, nz );
	vtk_bin_vectors ( fvor, "vorticity", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * 3 * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double vort_x = 0.0, vort_y = 0.0, vort_z = 0.0;

				if ( *meio )
				{
					double *f = lattice.inif + ( *meio - 1 ) * nvel;

					vorticity ( f, vort_x, vort_y, vort_z, x, y, z, geometry, lattice );
				}

				vtk_push_real ( buf, vort_x, double_prec );
				vtk_push_real ( buf, vort_y, double_prec );
				vtk_push_real ( buf, vort_z, double_prec );
			}
		}

		vtk_flush ( fvor, buf );
	}

	fvor << "\n";

	fvor.close();
}

//====================================================================================================================//



//=============================== Record the density field (binary) ==================================================//
//
//      Input: geometry, lattice, step, precisao
//      Output: rho_<passo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_density_bin ( GEOMETRY geometry, LATTICE lattice, unsigned int passo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string name_rho = "rho_" + to_string( passo ) + ".vtk";

	ofstream frho ( name_rho, ios::binary );

	vtk_bin_header  ( frho, "Densidade", nx, ny, nz );
	vtk_bin_scalars ( frho, "densidade", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double rho = 0.0;

				if ( *meio )
				{
					if ( lattice.inif != nullptr )
					{
						rho = density ( lattice.inif + ( *meio - 1 ) * nvel );
					}
					else
					{
						rho = density ( lattice.inif_R + ( *meio - 1 ) * nvel )
						    + density ( lattice.inif_B + ( *meio - 1 ) * nvel );
					}
				}

				vtk_push_real ( buf, rho, double_prec );
			}
		}

		vtk_flush ( frho, buf );
	}

	frho << "\n";

	frho.close();
}

//====================================================================================================================//



//=============================== Record the density field of one fluid (binary) =====================================//
//
//      Input: nome, geometry, distribuicoes, step, precisao
//      Output: <nome><passo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_density_bin ( string name, GEOMETRY geometry, double* ini_f, unsigned int passo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string name_rho = name + to_string( passo ) + ".vtk";

	ofstream frho ( name_rho, ios::binary );

	vtk_bin_header  ( frho, "Densidade", nx, ny, nz );
	vtk_bin_scalars ( frho, "densidade", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double rho = 0.0;

				if ( *meio ) rho = density ( ini_f + ( *meio - 1 ) * nvel );

				vtk_push_real ( buf, rho, double_prec );
			}
		}

		vtk_flush ( frho, buf );
	}

	frho << "\n";

	frho.close();
}

//====================================================================================================================//



//=============================== Record the phase field (binary) ====================================================//
//
//      Input: geometry, lattice, step, precisao
//      Output: phase_<passo>.vtk em VTK binario
//
//====================================================================================================================//

void rec_phase_bin ( GEOMETRY geometry, LATTICE lattice, unsigned int passo, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string name_phase = "phase_" + to_string( passo ) + ".vtk";

	ofstream file_phase ( name_phase, ios::binary );

	vtk_bin_header  ( file_phase, "Fase", nx, ny, nz );
	vtk_bin_scalars ( file_phase, "fase", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				double phase = 0.0;

				if ( *meio )
				{
					double rho_R = density ( lattice.inif_R + ( *meio - 1 ) * nvel );
					double rho_B = density ( lattice.inif_B + ( *meio - 1 ) * nvel );

					phase = ( rho_R - rho_B ) / ( rho_R + rho_B );
				}

				vtk_push_real ( buf, phase, double_prec );
			}
		}

		vtk_flush ( file_phase, buf );
	}

	file_phase << "\n";

	file_phase.close();
}

//====================================================================================================================//



//=============================== Record the geometry of the phases (binary) =========================================//
//
//      Input: geometry, lattice, step
//      Output: geometry_phase_<passo>.vtk em VTK binario, unsigned_char ( 0 solido, 1 azul, 2 vermelho )
//
//====================================================================================================================//

void rec_phase_geometry_bin ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	string name_phase = "geometry_phase_" + to_string( passo ) + ".vtk";

	ofstream file_phase ( name_phase, ios::binary );

	vtk_bin_header       ( file_phase, "Geometria das fases", nx, ny, nz );
	vtk_bin_scalars_byte ( file_phase, "fase" );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				unsigned char valor = 0;

				if ( *meio )
				{
					double rho_R = density ( lattice.inif_R + ( *meio - 1 ) * nvel );
					double rho_B = density ( lattice.inif_B + ( *meio - 1 ) * nvel );

					valor = ( ( rho_R - rho_B ) / ( rho_R + rho_B ) > 0 ) ? 2 : 1;
				}

				buf.push_back ( ( char ) valor );
			}
		}

		vtk_flush ( file_phase, buf );
	}

	file_phase << "\n";

	file_phase.close();
}

//====================================================================================================================//



//=============================== Record the geometry (binary) =======================================================//
//
//      Input: geometry, nome do arquivo
//      Output: <nome> em VTK binario, unsigned_char ( 1 fluido, 0 solido )
//
//====================================================================================================================//

void rec_geo_bin ( GEOMETRY geometry, string name_file )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	ofstream file_geo ( name_file, ios::binary );

	vtk_bin_header       ( file_geo, "Geometria", nx, ny, nz );
	vtk_bin_scalars_byte ( file_geo, "geometria" );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int *meio = geometry.ini + x + y * nx + z * ny * nx;

				buf.push_back ( ( *meio ) ? ( char ) 1 : ( char ) 0 );
			}
		}

		vtk_flush ( file_geo, buf );
	}

	file_geo << "\n";

	file_geo.close();
}

//====================================================================================================================//



//=============================== Record a scalar field (binary) =====================================================//
//
//      Input: geometry, nome do arquivo, campo escalar ( um valor por voxel ), precisao
//      Output: <nome> em VTK binario
//
//====================================================================================================================//

void rec_scalar_field_bin ( GEOMETRY geometry, string name_file, double* scalar, bool double_prec )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	ofstream file_scalar ( name_file, ios::binary );

	vtk_bin_header  ( file_scalar, "Campo escalar", nx, ny, nz );
	vtk_bin_scalars ( file_scalar, "escalar", double_prec );

	vector<char> buf;

	buf.reserve ( ( size_t ) nx * ny * ( double_prec ? 8 : 4 ) );

	for ( int z = 0; z < nz; z++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				vtk_push_real ( buf, scalar[ x + y * nx + z * nx * ny ], double_prec );
			}
		}

		vtk_flush ( file_scalar, buf );
	}

	file_scalar << "\n";

	file_scalar.close();
}

//====================================================================================================================//



//=============================== Records data recovery file =========================================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//====================================================================================================================//

void rec_recovery ( GEOMETRY geometry, LATTICE lattice, unsigned int passo, int precision )
{
    cout << "\n\nGravando arquivo de recuperação..." << endl;

    ofstream f_rec ( "Arq_rec.dat", ios::binary );

    f_rec << passo << endl;

    for ( int pto = 0; pto < geometry.fluid; pto++ )
    {
        double *f = lattice.inif + ( pto ) * nvel;

        double vx, vy, vz, rho;

        calcula ( f, vx, vy, vz, rho, lattice  );

        f_rec << setprecision( precision ) << vx << " ";

        f_rec << setprecision( precision ) << vy << " ";

        f_rec << setprecision( precision ) << vz << " ";

        f_rec << setprecision( precision ) << rho << " ";

    }

    f_rec.close();

    cout << "... ... ... !" << endl << endl;
}

//====================================================================================================================//




//=================================== Returns the density of a site ==================================================//
//
//      Input: distribution function
//      Output: density
//
//====================================================================================================================//

#pragma acc routine seq
double density ( double *f )
{

    double rho = 0.0;

    for ( int i = 0 ; i < nvel; i++ ) rho = rho + f[i];

    return rho;
}

//====================================================================================================================//







//=============================== Calculates the momentum in the x direction =========================================//
//
//      Input: distribution function f[nvel], lattice vectors
//      Output: momentum in the x direction
//
//====================================================================================================================//

#pragma acc routine seq
double quant_mov_x ( double f[nvel], LATTICE lattice )
{
    double mx = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
    }

    return mx;

}

//====================================================================================================================//




//=============================== Calculates the momentum in the y direction =========================================//
//
//      Input: distribution function f[nvel], lattice vectors
//      Output: momentum in the y direction
//
//====================================================================================================================//

#pragma acc routine seq
double quant_mov_y ( double f[nvel], LATTICE lattice )
{
    double my = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        my = my + lattice.c_i[ i * dim + 1 ] * f[i];
    }

    return my;

}

//====================================================================================================================//



//=============================== Calculates the momentum in the z direction =========================================//
//
//      Input: distribution function f[nvel], lattice vectors
//      Output: momentum in the z direction
//
//====================================================================================================================//

#pragma acc routine seq
double quant_mov_z ( double f[nvel], LATTICE lattice )
{
    double mz = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        mz = mz + lattice.c_i[ i * dim + 2 ] * f[i];
    }

    return mz;

}

//====================================================================================================================//



//================================ Calculates the force term =========================================================//
//
//      Input: velocities, Forces, velocities, density, relaxation time, lattice
//      Output: S (source)
//
//====================================================================================================================//

#pragma acc routine seq
void source ( double Fx, double Fy, double Fz, double vx, double vy, double vz, double rho, double tau, double S[nvel], 
				LATTICE lattice )
{
	if ( dim == 2 )
	{
		double one_ov_cs_sqd = lattice.one_over_c_s2;
		
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		
		double F[dim];
		
		F[0] = Fx;
		F[1] = Fy;
		
		double factor = ( 1. - 1. / ( 2. * tau ) );
			
		double vF = dot_product ( F, v );
		
		double c_i[dim];

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			  
			double cF = dot_product ( c_i, F );
			
			double cv = dot_product ( c_i, v );

			S[i] = factor * lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
		}
	}
	else
	{
		double one_ov_cs_sqd = lattice.one_over_c_s2;
		
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double F[dim];
		
		F[0] = Fx;
		F[1] = Fy;
		F[2] = Fz;
		
		double factor = ( 1. - 1. / ( 2. * tau ) );
			
		double vF = dot_product ( F, v );
		
		double c_i[dim];

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];           
			  
			double cF = dot_product ( c_i, F );
			
			double cv = dot_product ( c_i, v );

			S[i] = factor * lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
		}
	}
}

//====================================================================================================================//




//================================ Calculates the force term =========================================================//
//
//      Input: velocities, Forces, velocities, density, relaxation time, lattice
//      Output: S (source)
//
//====================================================================================================================//

#pragma acc routine seq
void force_Guo ( double Fx, double Fy, double Fz, double vx, double vy, double vz, double tau, double G[nvel], 
				LATTICE lattice )
{
	if ( dim == 2 )
	{
		double one_ov_cs_sqd = lattice.one_over_c_s2;
		
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		
		double F[dim];
		
		F[0] = Fx;
		F[1] = Fy;
		
		double factor = ( 1. - 1. / ( 2. * tau ) );
			
		double vF = dot_product ( F, v );
		
		double c_i[dim];

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			  
			double cF = dot_product ( c_i, F );
			
			double cv = dot_product ( c_i, v );

			G[i] = factor * lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
		}
	}
	else
	{
		double one_ov_cs_sqd = lattice.one_over_c_s2;
		
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double F[dim];
		
		F[0] = Fx;
		F[1] = Fy;
		F[2] = Fz;
		
		double factor = ( 1. - 1. / ( 2. * tau ) );
			
		double vF = dot_product ( F, v );
		
		double c_i[dim];

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];           
			  
			double cF = dot_product ( c_i, F );
			
			double cv = dot_product ( c_i, v );

			G[i] = factor * lattice.w[i] *  one_ov_cs_sqd * ( ( cF - vF) + one_ov_cs_sqd * cF * cv );
		}
	}
}

//====================================================================================================================//




//============================= Calculates the tensor Q_i = c_alpha c_beta - delta_alpha_beta ========================//
//
//      Input: lattice vectors
//      Output: tensor Q_i
//
//====================================================================================================================//

void calc_Q ( LATTICE lattice )
{
    for ( int i = 0 ; i < nvel; i++ )
    {
        for ( int alpha = 0; alpha < dim; alpha++ )
        {
            for ( int beta = 0; beta < dim; beta++ )
            {

                lattice.Q_i[ alpha + beta * dim + i * dim * dim ] 
                = lattice.c_i[ i*dim + alpha ] * lattice.c_i[ i*dim + beta ] - ( lattice.c_s2 ) * ( alpha == beta );
            }
        }
    }
}

//====================================================================================================================//




//============================= Calculates the tensor K2_i = c_alpha c_beta - delta_alpha_beta ========================//
//
//      Input: lattice vectors
//      Output: tensor K2_i
//
//====================================================================================================================//

void calc_K2 ( LATTICE lattice )
{
	double c_s2 = lattice.c_s2;
	
    for ( int i = 0 ; i < nvel; i++ )
    {
        for ( int alpha = 0; alpha < dim; alpha++ )
        {
            for ( int beta = 0; beta < dim; beta++ )
            {
                lattice.K2_i[ alpha + beta * dim + i * dim * dim ] = lattice.w[i] / ( 2. * c_s2 * c_s2 )
                 * ( lattice.c_i[ i*dim + alpha ] * lattice.c_i[ i*dim + beta ] - ( c_s2 ) * ( alpha == beta ) );
            }
        }
    }
}

//====================================================================================================================//




//===================== Imposição de tensão interfacial (Spencer, Halliday & Care) ===================================//
//
//      Input: ddistribution function, lattice vectors, factor dependent on the model,
//              unit vectors (normal to the interface)
//      Output: distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void interf_tension_SHC ( double *f, double fat, double* un, LATTICE lattice )
{
	//  A constante K do colchete e fixada pela conservacao de massa:
	//
	//      soma_i w_i [ ( n.c_i )^2 - c_i^2 + K ]  =  cs^2 - D cs^2 + K  =  0   =>   K = ( D - 1 ) cs^2
	//
	//  porque  soma_i w_i c_i_alfa c_i_beta = cs^2 delta_alfa_beta  ( traco = D cs^2 ) e |n| = 1.
	//  Em D3Q19, ( D - 1 ) cs^2 = 2/3 -- que era o valor cravado aqui.  Escrita assim a funcao
	//  tambem vale em 2D, onde a constante e 1/3, e o antigo  if ( dim == 3 )  deixa de ser preciso:
	//  antes, em D2Q9, esta funcao nao fazia NADA e o modelo rodava sem tensao interfacial, calado.
	//
	//  O momento se conserva para qualquer K: os momentos impares de terceira ordem se anulam.

	const double K = ( ( double ) dim - 1.0 ) * lattice.c_s2;

	for ( int i = 0; i < nvel; i++ )
	{
		//  c_i e lido direto de lattice.c_i.  Antes havia um  new double[dim]  aqui, com a copia
		//  das tres componentes: uma alocacao dinamica por sitio e por passo de tempo, dentro de
		//  uma  acc routine seq .

		const double *c_i = lattice.c_i + i * dim;

		const double prod_n_ci = dot_product ( ( double* ) c_i, un );

		const double c_2 = dot_product ( ( double* ) c_i, ( double* ) c_i );

		f[i] = f[i] + fat * lattice.w[i] * ( prod_n_ci * prod_n_ci - c_2 + K );
	}
}

//====================================================================================================================//




//===================== Imposição de tensão interfacial (Reis-Phillips/Leclaire) =====================================//
//
//      Input: distribution function, perturbation amplitude A, gradient magnitude, gradient vector and lattice
//      Output: distribution function with the RPL interfacial perturbation applied
//
//      D3Q19 form:
//          Omega_i^(2) = 0.5 * A * |grad(phi)| * [ w_i (c_i.grad(phi))^2 / |grad(phi)|^2 - B_i ]
//
//      With the usual D3Q19 normalization, sigma_lu ~= (2/9) A.
//
//====================================================================================================================//

#pragma acc routine seq
void interf_tension_RPL ( double *f, double A, double mod_grad, double* grad, LATTICE lattice )
{
    constexpr double eps = 1.0e-30;

    if ( nvel != 19 || dim != 3 ) return;
    if ( mod_grad <= eps || A == 0.0 ) return;

    constexpr double B[19] = {
        -2.0 / 9.0,
         1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
         1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
         1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0,
         1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0,
         1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0,
         1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0
    };

    double inv_mod_grad2 = 1.0 / ( mod_grad * mod_grad );
    double prefactor = 0.5 * A * mod_grad;

    for ( int i = 0; i < nvel; i++ )
    {
        double ci_dot_grad =
            lattice.c_i[ i * dim + 0 ] * grad[0] +
            lattice.c_i[ i * dim + 1 ] * grad[1] +
            lattice.c_i[ i * dim + 2 ] * grad[2];

        double directional_term = lattice.w[i] * ci_dot_grad * ci_dot_grad * inv_mod_grad2;

        f[i] += prefactor * ( directional_term - B[i] );
    }
}

//====================================================================================================================//



//===================== Imposição de tensão interfacial (Keijo, Philippi, Cirilo, Emerich) ===========================//
//
//      Input: ddistribution function, lattice vectors, factor dependent on the model,
//              unit vectors (normal to the interface)
//      Output: distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void interf_tension ( double *f, double fat, double tau, double* un, double mod_grad, LATTICE lattice )
{
	//  Mesma constante de conservacao de massa que em interf_tension_SHC:  K = ( D - 1 ) cs^2 .
	//  Em D3Q19 isso vale 2/3, que e o  2 * lattice.c_s2  que estava escrito -- so que aquele valor
	//  era o de tres dimensoes, cravado, e a funcao inteira estava dentro de um  if ( dim == 3 ) .
	//
	//  ATENCAO -- havia aqui um  new double[dim]  SEM o delete correspondente: um vazamento de
	//  dim * 8 bytes por sitio e por passo de tempo.  Numa caixa de 10^4 sitios rodando 10^5 passos,
	//  sao ~24 GB.  A alocacao foi removida ( c_i e lido direto de lattice.c_i ), o que resolve o
	//  vazamento e tira uma alocacao dinamica de dentro de uma  acc routine seq .

	const double rho = density ( f );

	const double fator = ( rho / ( 2. * tau ) ) * fat * mod_grad;

	const double K = ( ( double ) dim - 1.0 ) * lattice.c_s2;

	for ( int i = 0; i < nvel; i++ )
	{
		const double *c_i = lattice.c_i + i * dim;

		const double prod_n_ci = dot_product ( ( double* ) c_i, un );

		const double c_2 = dot_product ( ( double* ) c_i, ( double* ) c_i );

		f[i] = f[i] + fator * lattice.w[i] * ( prod_n_ci * prod_n_ci - c_2 + K );
	}
}

//====================================================================================================================//




//===================== Imposição de tensão interfacial (Gunstensen & Rothman) =======================================//
//
//      Input: ddistribution function, lattice vectors, factor dependent on the model,
//              unit vectors (normal to the interface)
//      Output: distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void interf_tension_GR ( double *f, double fat, double tau, double* un, double mod_grad, LATTICE lattice )
{		
	double fator = fat * mod_grad;
	
	if ( dim == 3 )
	{
		double c_i[dim];
		
		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];           
		
			double prod_n_ci = dot_product ( c_i, un );
			
			double gamma_i = fator * lattice.w[i] * ( prod_n_ci * prod_n_ci - lattice.c_s2 );
			
			f[i] = f[i] + gamma_i;
		}
	}
}

//====================================================================================================================//



//===================== Etapa de recoloração (Latva-Koko) ============================================================//
//
//      Input: distribution functions, lattice vectors, gradient, magnitude of gradient,
//              mass fractions, density, recolloring factor
//      Output: distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void recolloring ( double* f, double* f_R, double* f_B, double* grad, double mod_mM, double conc_R, double conc_B, 
					double rho, double beta, LATTICE lattice )
{
			
    double fator = beta * conc_R * conc_B;

	double c_i[3];
	
	double cos_phi = 0.;
	
    //----------------------------------------------------------------------------------------------------------------//
	
    f_R[0] = conc_R * f[0];
    f_B[0] = conc_B * f[0];

    for ( int i = 1; i < nvel; i++ )
    {

        c_i[0] = lattice.c_i[ i * dim + 0 ];
		c_i[1] = lattice.c_i[ i * dim + 1 ];
		c_i[2] = ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0;           
			  
		double c_2 = dot_product( c_i, c_i );

        double prod_vm_ci = dot_product( c_i, grad );

        if ( mod_mM ) cos_phi = prod_vm_ci / ( mod_mM * sqrt ( c_2 ) );

        else cos_phi = 0.0;

        f_R[i] = conc_R * f[i] - fator * lattice.w[i] * rho * cos_phi;
        f_B[i] = conc_B * f[i] + fator * lattice.w[i] * rho * cos_phi;
    }
}
//====================================================================================================================//



//======================== Computes the gradient of a scalar field ===================================================//
//
//      Input: distribution functions, lattice vectors, gradient, magnitude of gradient,
//              mass fractions, density, recolloring factor
//      Output: distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void gradient ( double* ini_psi, int x, int y, int z, double& grad_x, double& grad_y, double& grad_z, LATTICE lattice,
				GEOMETRY geometry )
                 
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    grad_x = 0.0;
    grad_y = 0.0;
    grad_z = 0.0;
    
    if ( dim == 2 )
    {
		for ( int i = 1; i < nvel; i++ )
		{
			int ci_x = ( int ) lattice.c_i[ i * dim + 0 ];
			int ci_y = ( int ) lattice.c_i[ i * dim + 1 ];

			int  pos_x = ( x + ci_x + nx ) % nx;
			int  pos_y = ( y + ci_y + ny ) % ny;

			double *psi = ini_psi + pos_x  +  pos_y * nx;

			grad_x = grad_x + *psi * ci_x * lattice.w[i];
			grad_y = grad_y + *psi * ci_y * lattice.w[i];
		}

		grad_x = grad_x * lattice.one_over_c_s2;
		grad_y = grad_y * lattice.one_over_c_s2;
	}
	else
	{
		for ( int i = 1; i < nvel; i++ )
		{
			int ci_x = ( int ) lattice.c_i[ i * dim + 0 ];
			int ci_y = ( int ) lattice.c_i[ i * dim + 1 ];
			int ci_z = ( int ) lattice.c_i[ i * dim + 2 ];

			int  pos_x = ( x + ci_x + nx ) % nx;
			int  pos_y = ( y + ci_y + ny ) % ny;
			int  pos_z = ( z + ci_z + nz ) % nz;

			double *psi = ini_psi + pos_x  +  pos_y * nx + pos_z * nx * ny;

			grad_x = grad_x + *psi * ci_x * lattice.w[i];
			grad_y = grad_y + *psi * ci_y * lattice.w[i];
			grad_z = grad_z + *psi * ci_z * lattice.w[i];
		}

		grad_x = grad_x * lattice.one_over_c_s2;
		grad_y = grad_y * lattice.one_over_c_s2;
		grad_z = grad_z * lattice.one_over_c_s2;
	}
}

//====================================================================================================================//



//======================== Computes the force term of the immiscible Shan-Chen model =================================//
//
//      Input: distribution functions, lattice vectors, gradient, magnitude of gradient,
//              mass fractions, density, recolloring factor
//      Output: distribution function
//
//====================================================================================================================//

//      Esta versao continua recebendo GEOMETRY e existe para nao quebrar quem ja a usava. Ela nao
//      pode ser chamada de dentro de um kernel: GEOMETRY contem um std::string e nao pode ser
//      copiada para o dispositivo. Para uso em GPU ha a sobrecarga logo abaixo, que recebe nx, ny,
//      nz e G soltos e e marcada com "acc routine seq".

void force_SC ( double* ini_psi, double* ini_psi_adj, int x, int y, int z, double& F_x, double& F_y, double& F_z, 
				LATTICE lattice, GEOMETRY geometry, PARAMETERS parameters )                 
{
	force_SC ( ini_psi, ini_psi_adj, x, y, z, geometry.nx, geometry.ny, geometry.nz,
	           parameters.A_fact, F_x, F_y, F_z, lattice );
}

//====================================================================================================================//




//=============================== Shan-Chen interaction force (device safe) ==========================================//
//
//      Input: pseudopotential of the component and of the other one, site, dimensions of the grid,
//             interaction strength, lattice
//      Output: force density acting on the component
//
//              F(x) = - G psi(x) soma_i w_i psi_adj(x + c_i) c_i
//
//      Devolve DENSIDADE de forca (forca por volume), que e o que a equacao de estado do modelo
//      pressupoe. coll_SC() recebe densidade de forca; ja coll_BGK() e op_bgk() recebem
//      ACELERACAO e fazem F = acc * rho internamente.
//
//      Os vizinhos sao lidos com envolvimento periodico e sem consultar a geometria: os voxels
//      solidos precisam ter um pseudopotencial definido pelo programa que chama -- e por ali que
//      entra a molhabilidade.
//
//====================================================================================================================//

#pragma acc routine seq
void force_SC ( double* ini_psi, double* ini_psi_adj, int x, int y, int z, int nx, int ny, int nz,
                double G, double& F_x, double& F_y, double& F_z, LATTICE lattice )
{
	
    double sum_psi_adj_cix = 0.;
    double sum_psi_adj_ciy = 0.;
    double sum_psi_adj_ciz = 0.;
    
	//  CORRECAO: faltava o termo z * nx * ny. Em 3D o pseudopotencial local era lido sempre no
	//  plano z = 0, enquanto os vizinhos, logo abaixo, ja eram endereçados corretamente. Compare
	//  com coll_phasetrans_BGK_SC, que faz a mesma conta certa.

    double *psi = ini_psi + x + y * nx + z * nx * ny;
    
    if ( dim == 2 )
    {
		for ( int i = 1; i < nvel; i++ )
		{
			int ci_x = ( int ) lattice.c_i[ i * dim + 0 ];
			int ci_y = ( int ) lattice.c_i[ i * dim + 1 ];

			int  pos_x = ( x + ci_x + nx ) % nx;
			int  pos_y = ( y + ci_y + ny ) % ny;

			double *psi_adj = ini_psi_adj + pos_x  +  pos_y * nx;

			sum_psi_adj_cix = sum_psi_adj_cix + *psi_adj * ci_x * lattice.w[i];
			sum_psi_adj_ciy = sum_psi_adj_ciy + *psi_adj * ci_y * lattice.w[i];
		}
		
		F_x = -psi[0] * G * sum_psi_adj_cix;
		F_y = -psi[0] * G * sum_psi_adj_ciy;
	}
	else
	{
		for ( int i = 1; i < nvel; i++ )
		{
			int ci_x = ( int ) lattice.c_i[ i * dim + 0 ];
			int ci_y = ( int ) lattice.c_i[ i * dim + 1 ];
			int ci_z = ( int ) lattice.c_i[ i * dim + 2 ];

			int  pos_x = ( x + ci_x + nx ) % nx;
			int  pos_y = ( y + ci_y + ny ) % ny;
			int  pos_z = ( z + ci_z + nz ) % nz;

			double *psi_adj = ini_psi_adj + pos_x  +  pos_y * nx + pos_z * nx * ny;

			sum_psi_adj_cix = sum_psi_adj_cix + *psi_adj * ci_x * lattice.w[i];
			sum_psi_adj_ciy = sum_psi_adj_ciy + *psi_adj * ci_y * lattice.w[i];
			sum_psi_adj_ciz = sum_psi_adj_ciz + *psi_adj * ci_z * lattice.w[i];
		}
		
		F_x = -psi[0] * G * sum_psi_adj_cix;
		F_y = -psi[0] * G * sum_psi_adj_ciy;
		F_z = -psi[0] * G * sum_psi_adj_ciz;
	}
}

//====================================================================================================================//




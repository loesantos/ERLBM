
//=============================== round a number =====================================================================//


int round_number ( double  num )
{
    int num_int;

    if ( num < 0 ) num_int = ( int ) ( num - 0.5 );
    else num_int = ( int ) ( num + 0.5 );

    return num_int;
}

//===================== Returns the angular coeficient of the best line (minimizing errors) ==========================//

double coefi ( double* y, int number , int npts, double k_std, double& avg )
{
	// k_std e apenas uma escala para tornar dy/dx adimensional; k_std = 0 (primeira chamada,
	// ou k em m^2 truncado por um antigo parametro inteiro) daria divisao por zero.

	if ( k_std == 0.0 ) k_std = 1.0;

	double sum_y = 0.0;
	double sum_x = 0.0;

	for ( int x = number; x < number + npts; x++ )
	{	
		sum_y = sum_y + y[x] / k_std;
		
		sum_x = sum_x + ( double ) x;
	}
    
	double x_avg = sum_x / ( double ) npts;    
    double y_avg = sum_y / ( double ) npts;
    
    double sum_xy_yavg = 0.0;    
    double sum_xx_xavg = 0.0;
    
    for ( int x = number; x < number + npts; x++ )
	{
		sum_xy_yavg = sum_xy_yavg + x * ( y[x] / k_std - y_avg );
		
		sum_xx_xavg = sum_xx_xavg + x * ( x - x_avg );
    }
            
	double dy_dx = sum_xy_yavg / sum_xx_xavg;
	
	avg = y_avg * k_std;

	return dy_dx;
        
}

//=============================== Calculate Absolut Permeability =====================================================//

double intrinsic_permeability ( double soma_mx, double Q_lost, double visc, int passo, double D_caract, double& k_mts,
								GEOMETRY geo )
{

    double mx_med = soma_mx / (double)( geo.nx );

    //------------------------------------------------------------------------------------//

    double Q_lost_med = Q_lost / (double)( geo.nx );
    //double Q_lost_med = sum_force / ( double ) ( geo.nx );

    //------------------------------------------------------------------------------------//

    k_mts = (geo.ftesc)*(geo.ftesc) * geo.phi * visc * mx_med / Q_lost_med;

    //------------------------------------------------------------------------------------//

    // converte m^2 para mDarcy
    double k_darcy = 10000000.0 * k_mts / ( 0.0000000098697 );

    cout << "\rStep : " << passo << "    k = " << k_darcy << "mDa;    k = " << k_mts
         << "m^2" << endl << endl;

    //------------------------------------------------------------------------------------//

    double Re = ( D_caract ) * ( mx_med / geo.fluid ) / visc;

    cout << "Reynolds = " << Re << endl << endl;

    ofstream fperme ("k_darcy.dat",ios::app);
    fperme << passo << " " << k_darcy << endl;
    fperme.close();

    ofstream fkm ("k_mts.dat",ios::app);
    fkm << passo << " " << k_mts << endl;
    fkm.close();
    
    return k_darcy;
}

//=============================== Computes the kinectic energy =======================================================//


double kinectic_energy ( GEOMETRY geometry, LATTICE lattice )
{
	double sum_energy = 0;
	
	for ( int pto = 0; pto < geometry.fluid; pto++ )
	{
		double *f = lattice.inif + ( pto ) * nvel;
		
		 double vx, vy, vz, rho;

         calcula ( f, vx, vy, vz, rho, lattice );
         
         sum_energy = sum_energy + 0.5 * ( vx*vx + vy*vy + vz*vz );
	}
	return sum_energy;
}

//=============================== Computes the enstrophy =============================================================//

double calc_enstrophy ( GEOMETRY geometry, LATTICE lattice )
{
	double sum_enstr = 0;
	
	for ( int x = 0; x < geometry.nx; x++ )
	{
		for ( int y = 0; y < geometry.ny; y++ )
		{
			for ( int z = 0; z < geometry.nz; z++ )
			{
				unsigned int pto = x + y * geometry.nx + z * geometry.nx * geometry.ny;
				
				double *f = lattice.inif + ( pto ) * nvel;
				
				double vort_x, vort_y, vort_z;
				
				vorticity ( f, vort_x, vort_y,vort_z, x, y, z, geometry, lattice );
				 
				sum_enstr = sum_enstr + 0.5 * ( vort_x*vort_x + vort_y*vort_y + vort_z*vort_z );
			}
		}
	}
	return sum_enstr;
}

//=============================== Returns the z position of the interface ============================================//

double pos_interf_y ( GEOMETRY geometry, LATTICE lattice, int x_h, int z_h )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	
	//------------ Determina a posição h(x) da interface ----------------//
	
	double sum_h = 0.;
	
	double sum_pond = 0.;
	
	int pos_y = 0;
	
	for ( int y = 0; y < ny; y++ ) 
	{		
		int pos = x_h + y * nx + z_h * nx * ny;
					
		int *geo = geometry.ini + pos;

		if ( geo[0] )
		{
			pos_y ++;
			
			//---------- Aponta os ponteiros ----------------------------//
			
			int pto = geo[0] - 1;
			
			double *f_R = lattice.inif_R + ( pto ) * nvel;			
			double *f_B = lattice.inif_B + ( pto ) * nvel;
			
			double rhoR = density ( f_R );
			double rhoB = density ( f_B );
			
			double concR = rhoR / ( rhoR + rhoB );
			double concB = rhoB / ( rhoR + rhoB );
			
			sum_h = sum_h + (double) pos_y * concR * concB; 
			
			sum_pond = sum_pond + concR * concB;
		}
	}
	
	double pos_h = sum_h / sum_pond;
	
	return pos_h;
}

//=============================== Returns the z position of the interface ============================================//

double pos_interf_z ( GEOMETRY geometry, LATTICE lattice, int x_h, int y_h )
{	 
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	//------------ Determina a posição h(z) da interface ----------------//
	
	double sum_h = 0.;
	
	double sum_pond = 0.;
	
	for ( int z = 5; z < ( nz - 5 ); z++ ) 
	{
		int pos = x_h + y_h * nx + z * nx * ny;
					
		int *geo = geometry.ini + pos;

		if ( geo[0] )
		{
			//---------- Aponta os ponteiros ----------------------------//
			
			int pto = geo[0] - 1;
			
			double *f_R = lattice.inif_R + ( pto ) * nvel;			
			double *f_B = lattice.inif_B + ( pto ) * nvel;
			
			double rhoR = density ( f_R );
			double rhoB = density ( f_B );
			
			double concR = rhoR / ( rhoR + rhoB );
			double concB = rhoB / ( rhoR + rhoB );
			
			if ( concR > 0.01 && concB > 0.01 )
			{			
				sum_h = sum_h + (double) z * concR * concB; 
			
				sum_pond = sum_pond + concR * concB;
			}
		}
	}
	
	double pos_h = sum_h / sum_pond;
	
	return pos_h;
}

//========================= Posição inicial da interface =============================================================//

double pos_ini_z ( GEOMETRY geometry, LATTICE lattice, int x_h, int y_h )
{	 
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	double concR_z = 0.;
	double concR_zPlus = 0.;
	
	double concB_z = 0.;
	double concB_zPlus = 0.;
	
	double pos_h = 0.;
	
	//------------ Determina a posição h(z) da interface ----------------//
	
	for ( int z = 1; z < (nz-2); z++ ) 
	{
		int pos = x_h + y_h * nx + z * nx * ny;
					
		int *geo = geometry.ini + pos;

		if ( geo[0] )
		{
			//---------- Aponta os ponteiros ----------------------------//
			
			int pto = geo[0] - 1;
			
			double *f_R = lattice.inif_R + ( pto ) * nvel;			
			double *f_B = lattice.inif_B + ( pto ) * nvel;
			
			double rhoR = density ( f_R );
			double rhoB = density ( f_B );
			
			concR_z = rhoR / ( rhoR + rhoB );
			concB_z = rhoB / ( rhoR + rhoB );
		}
		
		int posPlus = x_h + y_h * nx + ( z + 1 ) * nx * ny;
					
		int *geoPlus = geometry.ini + posPlus;

		if ( geoPlus[0] )
		{
			//---------- Aponta os ponteiros ----------------------------//
			
			int pto = geoPlus[0] - 1;
			
			double *f_R = lattice.inif_R + ( pto ) * nvel;			
			double *f_B = lattice.inif_B + ( pto ) * nvel;
			
			double rhoR = density ( f_R );
			double rhoB = density ( f_B );
			
			concR_zPlus = rhoR / ( rhoR + rhoB );
			concB_zPlus = rhoB / ( rhoR + rhoB );
		}
		
		if ( ( concR_z > concB_z && concR_zPlus < concB_zPlus ) || ( concR_z < concB_z && concR_zPlus > concB_zPlus ) )
		{
			pos_h = ( double ) z + 0.5;
			
			break;
		}
	}
	
	return pos_h;
}

//=============================== Returns the radius of a capillary tube =============================================//

double capillary_R ( GEOMETRY geometry, int z )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	
	int x_0 = nx / 2;
	int y_0 = ny / 2;
	
	double pi = 3.141592;
	
	int n_pix = ( int )( 2. * pi * nx/2 );
	
	double inc_theta = 2. * pi  / n_pix;
	
	double Radius = 0.;
	
	int flag = 0;
		
	for ( int R = 0; R < nx; R ++ )
	{
		for ( double theta = 0; theta <= 2*pi; theta = theta + 2 * inc_theta )
		{
			int x = x_0 + R * cos ( theta );
			int y = y_0 + R * sin ( theta );
				
			int pos = x + y * nx + z * nx * ny;
						
			int *geo = geometry.ini + pos;

			if ( geo[0] == 0 && flag == 0 )
			{
				Radius = R;
				
				flag = 1;
			}
		}
	}

	return Radius;
	
}

//=========================	Returns the contact angle (degrees) ======================================================//

double interf_angle ( GEOMETRY geometry, LATTICE lattice, double Radius, int step, int rest_time, double t_adm, 
						double delta_t, double x0, double y0 )
{
	double t_seg = ( step - rest_time ) * delta_t;

	double pi = acos(-1.);
	
	int ny = geometry.ny;
	
	int yd = ny / 2;
		
	//------------ Determina a posição h(z) da interface (no centro) -----------------------------------//
		
	double h1 = pos_interf_z ( geometry, lattice, x0, y0 );	
	
	//------------ Determina a posição h_d(z) da interface e theta_1 -----------------------------------//
	
	int xd = ( int ) ( Radius - 4. );
			
	double h2 = pos_interf_z ( geometry, lattice, x0 - xd, yd );
	
	double theta1 = contact_angle ( h1, h2, xd, Radius );
	
	xd = ( int ) (  Radius - 5. );
		
	h2 = pos_interf_z ( geometry, lattice, x0 - xd, yd );
	
	double theta2 = contact_angle ( h1, h2, xd, Radius );

	//------------ Determina e grava o ângulo de contato -----------------------------------------------//
	
	double theta = ( theta1 + theta2 ) / 2.; // média de duas medidas
			
	double theta_degrees = theta * ( 180 / pi );
		
	ofstream file_theta ( "theta.csv", ios::app );
	
	if ( step == 0 ) file_theta << "step, t_adm, t_seg, theta(deg), theta(rad)" << endl;
	
	else file_theta << step << "," << t_adm << "," << t_seg << "," << theta_degrees << "," << theta << endl;
	
	file_theta.close();
	
	return theta_degrees;	
}

//=========================	Returns the contact angle and position z of an interface  ================================//

double interf_height ( GEOMETRY geometry, LATTICE lattice, double Radius, int step, int rest_time, 
						double t_adm, double h_theo, double h_0, double delta_x, double delta_t, double x0, double y0 )
{
	double t_seg = ( step - rest_time ) * delta_t;
				
	//------------ Determina a posição h(z) da interface (x_0, y_0) -------------------------------------//
		
	double h1 = pos_interf_z ( geometry, lattice, x0, y0 );	
	
	//cout << "\nh1 = " << h1 << endl;
			
	//-------------- Grava a posição da interface -------------------------------------------------------//
    
    double h_adm = ( h1 - h_0 ) / h_theo;
		
	double h_mm =  ( h1 - h_0 ) * delta_x * 1000.;
		
    if ( step == 0 )  
	{		
		ofstream file_h ( "height.csv" );
		
		file_h << "step, t/T, t(s), h/H, h(mm)" << endl;
				
		file_h.close();
	}
	else
	{
		ofstream file_h ( "height.csv", ios::app );
		
		file_h << step << "," << t_adm << "," << t_seg << "," << h_adm << "," << h_mm  << endl;
		
		file_h.close();
	}
		
	return h_adm;	
}

//====================	Returns the contact angle (auxiliary function of interface)  =================================//
	
double contact_angle ( double h1, double h2, double xd, double Radius )
{
	double pi = 3.141592;
	
	double diff_h = h2 - h1;

	double R = ( xd * xd + diff_h * diff_h ) / ( 2. * diff_h );
		
	double sin_phi = Radius / R;
	
	if ( sin_phi > 1.0 ) sin_phi = 1.0;
	
	double theta = 0.5 * pi - asin( sin_phi );	
	
	if ( theta < 0. ) theta = 0.;
	
	return theta;
}

//============================= Determina a tensão interfacial =======================================================//
		
double calc_inter_tension_y ( GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    double sum_diff = 0.0;

    int x = nx / 2;
    int z = nz / 2;

    for ( int y = 0; y < ny; y++ )
    {
        int pos = x + y * nx + z * nx * ny;
						
		int *geo = geometry.ini + pos;

        if ( geo[0] )
        {
            //----------- Aponta os ponteiros --------------------------------------------//
			
			int pto = geo[0] - 1;
			
			double *R = lattice.inif_R + ( pto ) * nvel;			
			double *B = lattice.inif_B + ( pto ) * nvel;

            double Pi_xx = ( R[1] + R[2] + R[11] + R[12] + R[13] + R[14] )
							+ ( B[1] + B[2] + B[11] + B[12] + B[13]  + B[14] );

            double Pi_yy = ( R[3] + R[4] + R[15] + R[16] + R[17] + R[18] )
							+ ( B[3] + B[4] + B[15] + B[16] + B[17] + B[18] );

            double diff = Pi_yy - Pi_xx;

            sum_diff = sum_diff + diff / 2.;
        }
    }
    
    ofstream file_tens ( "Inter_tens.csv", ios::app );
	
    file_tens << sum_diff << endl;
    
    file_tens.close();
    
    return sum_diff;
}

//============================= Determina velocidade da interface em z ===============================================//
		
double interf_vel_z ( GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    double sum_vz = 0.0;
    
    int sum_ptos = 0;
    
    for ( int x = 0; x < nx; x ++ )
    {
		for ( int y = 0; y < ny; y ++ )
		{
			for ( int z = 5; z < nz - 5; z ++ )
			{
				int pos = x + y * nx + z * nx * ny;
						
				int *geo = geometry.ini + pos;

				if ( geo[0] )
				{
					//----------- Aponta os ponteiros --------------------------------------------//
					
					int pto = geo[0] - 1;
					
					double *f_R = lattice.inif_R + ( pto ) * nvel;
					double *f_B = lattice.inif_B + ( pto ) * nvel;
			   
					double rho_R = density ( f_R );
					double rho_B = density ( f_B );

					double conc_R = rho_R / ( rho_R + rho_B );
					double conc_B = 1.0 - conc_R;
					
					double vz = 0.;
					
					if ( conc_B < 1.e-5 )
					{
						vz = quant_mov_z ( f_R, lattice ) / rho_R;
						
						sum_vz = sum_vz + vz;
						
						sum_ptos ++;
					}
		
					if ( conc_R < 1.e-5 )
					{
						vz = quant_mov_z ( f_B, lattice ) / rho_B;
						
						sum_vz = sum_vz + vz;
						
						sum_ptos ++;
					}
				}
			}
		}
    }
    
    double vel_int = sum_vz / ( double ) sum_ptos; 
    
    return vel_int;
}

//============================= Grava o número capilar ===============================================================//

void rec_Ca_number ( int step, int rest_time, double delta_x, double delta_t, double vel_itf, double sigma_sim, 
						PARAMETERS parameters )
{
	const double mu_sim = parameters.visc_R * parameters.rho_ini_R;
	
	double t_seg = ( step - rest_time ) * delta_t;
	
	double Ca_number = vel_itf * mu_sim / sigma_sim;
	
	double vel_mm = 1000. * vel_itf * delta_x / delta_t; // velocidade em mm/s
	
	if ( step == 0 ) 
	{
		ofstream file_Ca ( "Ca.csv" );
		
		file_Ca << "step, t_seg, Ca, v(mm/s)" << endl;
	    
		file_Ca.close();
	}
	else	
	{	
		ofstream file_Ca ( "Ca.csv", ios::app );
	
		file_Ca << step << "," << t_seg << "," << Ca_number << "," << vel_mm << endl;
    
		file_Ca.close();
	}
}

//============================= Determina a tensão interfacial =======================================================//
		
double calc_inter_tension_x ( GEOMETRY geometry, LATTICE lattice, PARAMETERS parameters )
{
	int *geo = geometry.ini;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	double sum_diff = 0.0;

	int y = ny / 2;
	int z = nz / 2;

	for ( int x = nx / 10; x < 9 * nx / 10; x++ )
	{
		int pos = x + y * nx + z * nx * ny;	
				
		if ( geo[pos] )
		{
			//----------- Aponta os ponteiros --------------------------------------------//

			double *f_R = lattice.inif_R + ( geo[pos] - 1 ) * nvel;
			double *f_B = lattice.inif_B + ( geo[pos] - 1 ) * nvel;	
			
			double Pi_xx = 0.0;
			double Pi_yy = 0.0;				
			
			for ( int i = 0; i < nvel; i++ )
			{
				double* c =  lattice.c_i + i * dim;
			
				Pi_xx += ( f_R[i] + f_B[i] ) * c[0] * c[0];
				Pi_yy += ( f_R[i] + f_B[i] ) * c[1] * c[1];
			}
			
			double diff = Pi_xx - Pi_yy;
			
			sum_diff = sum_diff + diff;
		}
	}
	
	ofstream file_tens ( "Inter_tens.csv", ios::app );
	
    file_tens << sum_diff << endl;
    
    file_tens.close();
	
	return 	sum_diff;	
}


//=============================== Returns the total momentum x - direction ===========================================//

double total_mx ( GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	double sum_mx = 0.;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {				
                int *meio = geometry.ini + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						sum_mx = sum_mx + quant_mov_x ( f, lattice );		
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						sum_mx = sum_mx + quant_mov_x ( f_R, lattice ) + quant_mov_x ( f_B, lattice );
					}
                }
            }
        }
    }
    
    return sum_mx;
}

//=============================== Returns the total momentum y - direction ===========================================//

double total_my ( GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	double sum_my = 0.;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {				
                int *meio = geometry.ini + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						sum_my = sum_my + quant_mov_y ( f, lattice );		
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						sum_my = sum_my + quant_mov_y ( f_R, lattice ) + quant_mov_y ( f_B, lattice );
					}
                }
            }
        }
    }
    
    return sum_my;
}

//=============================== Returns the total momentum z - direction ===========================================//

double total_mz ( GEOMETRY geometry, LATTICE lattice )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	double sum_mz = 0.;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {				
                int *meio = geometry.ini + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						sum_mz = sum_mz + quant_mov_z ( f, lattice );		
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						sum_mz = sum_mz + quant_mov_z ( f_R, lattice ) + quant_mov_z ( f_B, lattice );
					}
                }
            }
        }
    }
    
    return sum_mz;
}

//=============================== Medidas de uma gota ================================================================//
//
//      calc_drop_radius     : centro, raio e correntes espurias de uma gota isolada
//      calc_laplace_tension : salto de pressao e tensao interfacial pela lei de Laplace
//
//      As duas trabalham sobre a struct DROP ( Definitions.cpp ), que leva as entradas
//      ( cilindro, frac_in, frac_out, espessura ) e recebe todas as saidas.  Rodam no HOSPEDEIRO:
//      num programa OpenACC, faca 'acc update self' de inif_R e inif_B antes de chamar.  O custo e
//      desprezivel perto dos milhares de passos entre dois diagnosticos.
//
//====================================================================================================================//

//=============================== Centro, raio e correntes espurias de uma gota ======================================//
//
//      Preenche, em DROP:
//
//          x0, y0, z0    centro da gota, pelo CENTROIDE CIRCULAR de phi = rho_R / ( rho_R + rho_B ).
//                        Circular porque a caixa e periodica: a media aritmetica de x erraria assim
//                        que a gota encostasse na borda.
//
//          area          soma_sitios phi.  Como soma_sitios rho_R e conservada exatamente, esta e a
//                        medida estavel do tamanho da gota.
//
//          raio          raio equimolar:  raiz( area / ( pi nz ) )  ( cilindro )
//                                         ( 3 area / 4 pi )^(1/3)   ( esfera )
//
//          raio_grad     primeiro momento radial de rho_R rho_B sobre o dominio inteiro: o raio
//                        medio da superficie de tensao.  Estimador independente do de cima; os dois
//                        concordam dentro de ( espessura / raio )^2.  Se divergirem, a gota nao esta
//                        redonda ou nao convergiu.
//
//          raio_linha    o estimador de uma linha so ( primeiro momento de rho_R rho_B ao longo da
//                        linha que passa pelo centro ), guardado por compatibilidade com programas
//                        antigos.  E FRAGIL: usa poucos sitios e soma a linha inteira, inclusive o
//                        campo distante, onde rho_R rho_B e minusculo mas o braco de alavanca e
//                        grande.  Nao use para calibrar nada.
//
//          massa_R, massa_B, qx, qy, qz, u_max, u_rms, n_ruim
//
//      Entradas usadas: drop.cilindro.
//
//      Retorna drop.raio.
//
//====================================================================================================================//

double calc_drop_radius ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	int n_voxels = nx * ny * nz;

	const double pi = 3.14159265358979323846;

	const double tx = 2.0 * pi / ( double ) nx;
	const double ty = 2.0 * pi / ( double ) ny;
	const double tz = 2.0 * pi / ( double ) nz;

	//--------------- Passo 1: massas, area, centroide e correntes espurias ---------------------------//

	double massa_R = 0.0, massa_B = 0.0;
	double qx = 0.0, qy = 0.0, qz = 0.0;
	double area = 0.0, u2_sum = 0.0, u_max = 0.0;
	double cx_c = 0.0, cx_s = 0.0, cy_c = 0.0, cy_s = 0.0, cz_c = 0.0, cz_s = 0.0;

	int n_ruim = 0, n_fluid = 0;

	#pragma omp parallel for \
		reduction(+:massa_R,massa_B,qx,qy,qz,area,u2_sum,cx_c,cx_s,cy_c,cy_s,cz_c,cz_s,n_ruim,n_fluid) \
		reduction(max:u_max)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		int pto = geo[pos] - 1;

		double vxR, vyR, vzR, rho_R;
		double vxB, vyB, vzB, rho_B;

		calcula ( lattice.inif_R + pto * nvel, vxR, vyR, vzR, rho_R, lattice );
		calcula ( lattice.inif_B + pto * nvel, vxB, vyB, vzB, rho_B, lattice );

		n_fluid = n_fluid + 1;

		//  Teste de finitude escrito sem isnan() para poder ser reaproveitado em codigo de
		//  dispositivo compilado com matematica rapida.

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) || ! ( vxR == vxR ) || ! ( vxB == vxB ) )
		{
			n_ruim = n_ruim + 1;

			continue;
		}

		massa_R = massa_R + rho_R;
		massa_B = massa_B + rho_B;

		qx = qx + rho_R * vxR + rho_B * vxB;
		qy = qy + rho_R * vyR + rho_B * vyB;
		qz = qz + rho_R * vzR + rho_B * vzB;

		double rho = rho_R + rho_B;

		if ( rho <= 0.0 ) continue;

		double ux = ( rho_R * vxR + rho_B * vxB ) / rho;
		double uy = ( rho_R * vyR + rho_B * vyB ) / rho;
		double uz = ( rho_R * vzR + rho_B * vzB ) / rho;

		//  CORRECAO DA MEIA-FORCA.  Num esquema forcado a Guo a velocidade fisica NAO e o momento
		//  cru das populacoes:
		//
		//      antes da colisao    rho u = soma_i f_i c_i  +  F / 2
		//      depois da colisao   rho u = soma_i f_i c_i  -  F / 2
		//
		//  Esta medida roda DEPOIS da colisao, entao o sinal e o de baixo.  Sem isto, uma gota em
		//  equilibrio perfeito ( u = 0 em toda parte ) aparece com uma "corrente espuria" igual a
		//  | F | / 2 rho -- radial, axissimetrica, concentrada na interface e proporcional a
		//  tensao interfacial.  Nao e corrente nenhuma: e a metade da forca.
		//
		//  Medido na lei de Laplace, SHC com forcamento de Guo, R = 40, beta = 0.8, sigma = 5.36e-03:
		//  |u|max caiu de 3.54e-04 ( = |F|max / 2 rho a 0.7% ) para 1.97e-05 quando esta linha
		//  passou a existir.  Um fator 18.
		//
		//  lattice.ini_force = nullptr ( o padrao ) desliga a correcao: os modelos que nao usam
		//  forca de corpo -- Gunstensen, Santos, a versao H -- nao precisam dela.

		if ( lattice.ini_force != nullptr )
		{
			const double *F = lattice.ini_force + ( long long ) pto * 3;

			ux = ux - 0.5 * F[0] / rho;
			uy = uy - 0.5 * F[1] / rho;
			uz = uz - 0.5 * F[2] / rho;

			qx = qx - 0.5 * F[0];
			qy = qy - 0.5 * F[1];
			qz = qz - 0.5 * F[2];
		}

		double u2 = ux * ux + uy * uy + uz * uz;

		u2_sum = u2_sum + u2;

		double u = sqrt ( u2 );

		if ( u > u_max ) u_max = u;

		double phi = rho_R / rho;

		area = area + phi;

		int z = pos / ( nx * ny );
		int y = ( pos - z * nx * ny ) / nx;
		int x = pos - z * nx * ny - y * nx;

		cx_c = cx_c + phi * cos ( tx * x );   cx_s = cx_s + phi * sin ( tx * x );
		cy_c = cy_c + phi * cos ( ty * y );   cy_s = cy_s + phi * sin ( ty * y );
		cz_c = cz_c + phi * cos ( tz * z );   cz_s = cz_s + phi * sin ( tz * z );
	}

	double x0 = atan2 ( cx_s, cx_c ) / tx;   if ( x0 < 0.0 ) x0 += nx;
	double y0 = atan2 ( cy_s, cy_c ) / ty;   if ( y0 < 0.0 ) y0 += ny;
	double z0 = atan2 ( cz_s, cz_c ) / tz;   if ( z0 < 0.0 ) z0 += nz;

	if ( nz == 1 ) z0 = 0.0;

	double raio;

	if ( drop.cilindro ) raio = ( area > 0.0 ) ? sqrt ( area / ( pi * ( double ) nz ) ) : 0.0;

	else                 raio = ( area > 0.0 ) ? cbrt ( 3.0 * area / ( 4.0 * pi ) )     : 0.0;

	//--------------- Passo 2: raio da interface, pelo momento de rho_R rho_B -------------------------//

	double som_w = 0.0, som_wr = 0.0;
	double sx_e = 0.0, sw_e = 0.0, sx_d = 0.0, sw_d = 0.0;

	//  A LINHA DO ESTIMADOR raio_linha.
	//
	//  Num CILINDRO o campo e uniforme em z e o centroide circular em z e indefinido: cz_s e cz_c
	//  sao exatamente zero no papel e ruido de arredondamento ( ~1e-13 ) no computador, e o atan2
	//  de dois ruidos devolve um angulo qualquer em [ 0, nz ).  Com z0 > nz - 0.5 o arredondamento
	//  dava  z_lin = nz , que nao e plano nenhum: nenhum sitio casava, sw_e = sw_d = 0 e tanto
	//  raio_linha quanto sigma_2pontos saiam ZERO -- de forma intermitente, so em alguns passos.
	//  Num cilindro qualquer plano serve, entao fixamos z = 0.  Na esfera, o clamp evita o mesmo
	//  arredondamento na borda.

	int y_lin = ( int ) ( y0 + 0.5 );   if ( y_lin >= ny ) y_lin -= ny;

	int z_lin = ( nz == 1 || drop.cilindro ) ? 0 : ( int ) ( z0 + 0.5 );   if ( z_lin >= nz ) z_lin -= nz;

	#pragma omp parallel for reduction(+:som_w,som_wr,sx_e,sw_e,sx_d,sw_d)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		int pto = geo[pos] - 1;

		double rho_R = density ( lattice.inif_R + pto * nvel );
		double rho_B = density ( lattice.inif_B + pto * nvel );

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) ) continue;

		int z = pos / ( nx * ny );
		int y = ( pos - z * nx * ny ) / nx;
		int x = pos - z * nx * ny - y * nx;

		//  Imagem minima: a caixa e periodica.

		double dx = ( double ) x - x0;   if ( dx >  0.5 * nx ) dx -= nx;   if ( dx < -0.5 * nx ) dx += nx;
		double dy = ( double ) y - y0;   if ( dy >  0.5 * ny ) dy -= ny;   if ( dy < -0.5 * ny ) dy += ny;
		double dz = ( double ) z - z0;   if ( dz >  0.5 * nz ) dz -= nz;   if ( dz < -0.5 * nz ) dz += nz;

		if ( drop.cilindro ) dz = 0.0;

		double r = sqrt ( dx * dx + dy * dy + dz * dz );

		double w = rho_R * rho_B;

		som_w  = som_w  + w;
		som_wr = som_wr + w * r;

		if ( y == y_lin && z == z_lin )
		{
			if ( dx < 0.0 ) { sx_e = sx_e + w * ( double ) x;   sw_e = sw_e + w; }

			else            { sx_d = sx_d + w * ( double ) x;   sw_d = sw_d + w; }
		}
	}

	//--------------- Devolve --------------------------------------------------------------------------//

	drop.x0 = x0;   drop.y0 = y0;   drop.z0 = z0;

	drop.area = area;

	drop.raio = raio;

	drop.raio_grad = ( som_w > 0.0 ) ? som_wr / som_w : 0.0;

	drop.raio_linha = ( sw_e > 0.0 && sw_d > 0.0 ) ? 0.5 * ( sx_d / sw_d - sx_e / sw_e ) : 0.0;

	drop.massa_R = massa_R;   drop.massa_B = massa_B;

	drop.qx = qx;   drop.qy = qy;   drop.qz = qz;

	drop.u_max = u_max;

	drop.u_rms = ( n_fluid > 0 ) ? sqrt ( u2_sum / ( double ) n_fluid ) : 0.0;

	drop.n_fluid = n_fluid;

	drop.n_ruim = n_ruim;

	return raio;
}

//=============================== Tensao interfacial pela lei de Laplace =============================================//
//
//      Chama calc_drop_radius e acrescenta o salto de pressao e a tensao:
//
//          cilindro ( 2D ) :  sigma = R delta_p          esfera ( 3D ) :  sigma = R delta_p / 2
//
//      A pressao e p = cs^2 ( rho_R + rho_B ).  No seio de cada fase a outra especie decai
//      exponencialmente, entao o termo cruzado se anula e essa e mesmo a pressao termodinamica.
//
//          p_in  : media sobre  r <  frac_in  * raio
//          p_out : media sobre  r >= frac_out * raio
//
//      MEDIAS SOBRE REGIOES, e nao dois pontos: as correntes espurias deixam uma pressao dinamica da
//      ordem de rho |u|^2 -- dentro de uma gota pequena isso chega a alguns por cento do proprio
//      salto -- e o centro da gota costuma ser justamente um ponto de estagnacao, onde ela e maxima.
//      Os desvios padrao das duas regioes ( sd_in, sd_out ) sao gravados: sao a barra de erro.
//
//      A regiao externa e limitada duas vezes: por baixo, para nao encostar na interface
//      ( raio + 4 * drop.espessura ); por cima, para caber na caixa ( 0.45 * min(nx,ny) ).
//
//      p_centro e p_canto guardam a medida de DOIS PONTOS ( centro da gota e voxel 0 ), e
//      sigma_linha o sigma que sairia dela com raio_linha -- so para comparar com resultados antigos
//      dentro da mesma rodada.
//
//      Entradas usadas: drop.cilindro, drop.frac_in, drop.frac_out, drop.espessura.
//
//      Retorna drop.sigma.
//
//====================================================================================================================//

//=============================== Tensao interfacial de uma interface PLANA =========================================//
//
//      Definicao mecanica ( Rowlinson & Widom; Latva-Kokko & Rothman, eq. 16 ):
//
//          sigma = integral ( P_N - P_T ) dx
//
//      com a interface normal a x,  P_N = P_xx  e  P_T = ( P_yy + P_zz ) / 2 .  O tensor de pressao
//      e o segundo momento das populacoes,  P_ab = soma_i f_i c_ia c_ib , que em repouso ( u = 0 )
//      e exatamente o tensor de pressao.
//
//      E a medida COMPLEMENTAR a lei de Laplace, e melhor em dois pontos:
//
//        * nao ha curvatura, entao nao ha correcao de Tolman nem extrapolacao em 1/R -- o valor sai
//          direto, de uma rodada so;
//        * a integral de | grad( fase ) | atraves da interface e uma REGRA DE SOMA exata: vale a
//          variacao total da fase entre os dois lados ( 2, para rho^N indo de -1 a +1 ), qualquer
//          que seja o perfil.  Latva-Kokko & Rothman usam isso na deducao da eq. (21).  Aqui ela
//          serve de conferencia: integral_grad tem de dar 2 na precisao da maquina.
//
//      ATENCAO -- E ELA QUE DECIDE SE ESTA MEDIDA SE APLICA:
//
//        * modelos em que a tensao entra NAS POPULACOES, por operador de perturbacao
//          ( Latva-Kokko, Liu-Valocchi-Kang, Reis-Phillips, Li-Yu-Luo, SHC_FONTE, SHC_H ):
//          a anisotropia esta em f, e  sigma_pop  e a resposta;
//
//        * modelos de FORCA DE CORPO ( Saito, SHC_ARTIGO, SHC_EQ56 ): numa interface plana a
//          curvatura e ZERO.  Na rota do continuum surface force isso significa forca zero, e as
//          populacoes nao carregam anisotropia nenhuma -- sigma_pop da zero, e esta certo.  A tensao
//          desses modelos vive no TENSOR, e quem a integra e o proprio programa, que sabe o que
//          seus campos significam.  A biblioteca nao adivinha: ela reporta sigma_pop e a regra de
//          soma, e o programa preenche drop.sigma quando tem outra rota.
//
//      A janela de integracao e  x em [ nx/4 , 3nx/4 ) : a caixa periodica tem DUAS interfaces, e
//      esta janela contem so a central.  Tudo e media sobre a secao transversal ( ny x nz ), o que
//      com ny = nz = 1 recai no caso unidimensional exato.
//
//      Preenche em DROP:  sigma_pop, integral_grad, sigma ( = sigma_pop ), massa_R/B, qx/qy/qz,
//      u_max, u_rms, n_fluid, n_ruim.  O perfil sai de calc_perfil_logistico() com drop.plana.
//
//      Devolve sigma_pop.
//
//====================================================================================================================//

double calc_flat_tension ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	constexpr double eps = 1.0e-30;

	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	const int n_voxels = nx * ny * nz;

	const int *geo = geometry.ini;

	//------------------ Massas, momento e correntes, como em calc_drop_radius ----------------------//

	double massa_R = 0.0, massa_B = 0.0;
	double qx = 0.0, qy = 0.0, qz = 0.0;
	double u2_sum = 0.0, u_max = 0.0;

	int n_fluid = 0, n_ruim = 0;

	#pragma omp parallel for \
		reduction(+:massa_R,massa_B,qx,qy,qz,u2_sum,n_fluid,n_ruim) reduction(max:u_max)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		const int pto = geo[pos] - 1;

		double vxR, vyR, vzR, rho_R;
		double vxB, vyB, vzB, rho_B;

		calcula ( lattice.inif_R + pto * nvel, vxR, vyR, vzR, rho_R, lattice );
		calcula ( lattice.inif_B + pto * nvel, vxB, vyB, vzB, rho_B, lattice );

		n_fluid = n_fluid + 1;

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) || ! ( vxR == vxR ) || ! ( vxB == vxB ) )
		{
			n_ruim = n_ruim + 1;

			continue;
		}

		massa_R = massa_R + rho_R;
		massa_B = massa_B + rho_B;

		const double rho = rho_R + rho_B;

		if ( rho <= eps ) continue;

		double ux = ( rho_R * vxR + rho_B * vxB ) / rho;
		double uy = ( rho_R * vyR + rho_B * vyB ) / rho;
		double uz = ( rho_R * vzR + rho_B * vzB ) / rho;

		//  Mesma correcao da meia-forca de calc_drop_radius(): pos-colisao,
		//  rho u = soma_i f_i c_i - F/2 .  Sem ela o que aparece como corrente espuria e F/2rho.

		if ( lattice.ini_force != nullptr )
		{
			const double *F = lattice.ini_force + ( long long ) pto * 3;

			ux = ux - 0.5 * F[0] / rho;
			uy = uy - 0.5 * F[1] / rho;
			uz = uz - 0.5 * F[2] / rho;
		}

		qx = qx + rho * ux;
		qy = qy + rho * uy;
		qz = qz + rho * uz;

		const double u2 = ux * ux + uy * uy + uz * uz;

		u2_sum = u2_sum + u2;

		const double u = sqrt ( u2 );

		if ( u > u_max ) u_max = u;
	}

	//------------------ A integral mecanica, so na interface central -------------------------------//

	const int x_ini = nx / 4;
	const int x_fim = 3 * nx / 4;

	double soma_aniso = 0.0, soma_grad = 0.0;

	int n_col = 0;

	#pragma omp parallel for reduction(+:soma_aniso,soma_grad,n_col)

	for ( int x = x_ini; x < x_fim; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				const int pos = x + y * nx + z * nx * ny;

				if ( geo[pos] == 0 ) continue;

				const int pto = geo[pos] - 1;

				const double *f_R = lattice.inif_R + pto * nvel;
				const double *f_B = lattice.inif_B + pto * nvel;

				double P_xx = 0.0, P_yy = 0.0, P_zz = 0.0;

				for ( int i = 0; i < nvel; i++ )
				{
					const double *c = lattice.c_i + i * dim;

					const double f = f_R[i] + f_B[i];

					P_xx = P_xx + f * c[0] * c[0];
					P_yy = P_yy + f * c[1] * c[1];
					P_zz = P_zz + f * c[2] * c[2];
				}

				soma_aniso = soma_aniso + P_xx - 0.5 * ( P_yy + P_zz );

				//  Regra de soma:  integral | grad( fase ) | dx = variacao total da fase.

				if ( lattice.inif_m != nullptr )
				{
					double gx, gy, gz;

					momentum ( lattice.inif_m + pto * nvel, gx, gy, gz, lattice );

					soma_grad = soma_grad + sqrt ( gx * gx + gy * gy + gz * gz );
				}

				n_col = n_col + 1;
			}
		}
	}

	//  Media sobre a secao transversal: o que se quer e a integral por COLUNA.

	const double n_sec = ( n_col > 0 ) ? ( double ) n_col / ( double ) ( x_fim - x_ini ) : 1.0;

	drop.sigma_pop     = soma_aniso / n_sec;
	drop.integral_grad = soma_grad  / n_sec;

	drop.sigma = drop.sigma_pop;

	//  Positividade das populacoes. A interface plana perfeitamente simetrica pode continuar com
	//  u = 0 mesmo para amplitudes interfaciais absurdas; nesse caso a perda de positividade e um
	//  criterio conservador, mais informativo que a simples ausencia de NaN.

	double min_f_R = 1.0e300, min_f_B = 1.0e300, min_f_total = 1.0e300;

	for ( int pto = 0; pto < n_fluid; pto++ )
	{
		const double *f_R = lattice.inif_R + pto * nvel;
		const double *f_B = lattice.inif_B + pto * nvel;

		for ( int i = 0; i < nvel; i++ )
		{
			if ( f_R[i] < min_f_R ) min_f_R = f_R[i];
			if ( f_B[i] < min_f_B ) min_f_B = f_B[i];

			const double f = f_R[i] + f_B[i];

			if ( f < min_f_total ) min_f_total = f;
		}
	}

	drop.min_f_R = min_f_R;
	drop.min_f_B = min_f_B;
	drop.min_f_total = min_f_total;

	drop.massa_R = massa_R;   drop.massa_B = massa_B;
	drop.qx = qx;   drop.qy = qy;   drop.qz = qz;
	drop.u_max = u_max;
	drop.u_rms = ( n_fluid > 0 ) ? sqrt ( u2_sum / ( double ) n_fluid ) : 0.0;
	drop.n_fluid = n_fluid;
	drop.n_ruim  = n_ruim;

	//  A medida do perfil usa o mesmo centro que a laje: o centro do vermelho e nx/4.

	drop.plana = true;
	drop.x0 = 0.25 * ( double ) nx;
	drop.y0 = 0.0;
	drop.z0 = 0.0;

	return drop.sigma_pop;
}

//====================================================================================================================//


//=============================== Curva interfacial de Latva-Kokko & Rothman ( 2005 ) ================================//
//
//      Latva-Kokko & Rothman, Phys. Rev. E 71, 056702 (2005), eqs. (15)-(17).  A recoloracao
//      daqueles autores nao so separa as cores: ela leva o perfil a uma forma FECHADA.  O balanco
//      entre a difusao de cor e o termo de segregacao da
//
//          d phi / ds = K phi ( 1 - phi )                                        eq. (15)
//
//      cuja solucao e uma logistica
//
//          phi( s ) = 1 / ( 1 + exp[ K ( s - s_0 ) ] )                           eq. (17)
//
//      com  K  uma constante de comprimento INVERSO, proporcional a beta.  O artigo mostra isso no
//      modelo D1Q2 ( onde K = 4 beta / N ) e afirma que a FORMA vale para qualquer rede que use a
//      recoloracao das eqs. (9)-(10) -- so o prefator K( beta ) muda.  E o que esta funcao mede.
//
//      Aqui  s  e a distancia radial ao centro da gota, e o perfil sai da media azimutal de
//      phi = rho_R / rho  sobre cascas de meio sitio de largura.  O ajuste e linear:
//
//          ln[ phi / ( 1 - phi ) ] = - K ( r - r_0 )
//
//      entao a inclinacao da K e o zero da o raio da interface.  So entram as camadas com
//      0.02 < phi < 0.98 : fora dessa faixa phi e dominado pela cauda exponencial e pelo ruido de
//      arredondamento, e o logito estoura.
//
//      O que sai em DROP:
//
//          K_perfil    K medido
//          r_perfil    raio onde phi = 1/2 ( terceiro estimador independente de raio )
//          res_perfil  maior residuo de ln[phi/(1-phi)] -- e ele que diz se a logistica cabe
//          esp_perfil  largura 10%-90% implicada, 2 ln(9) / K
//          n_perfil    quantas camadas entraram no ajuste
//
//      Devolve K_perfil.  Roda no hospedeiro, com OpenMP, como as outras medidas.
//
//====================================================================================================================//

double calc_perfil_logistico ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	constexpr double eps = 1.0e-30;

	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	const int n_voxels = nx * ny * nz;

	const int *geo = geometry.ini;

	//  Cascas de meio sitio.  O alcance vai ate a diagonal da caixa, o que basta e sobra.

	const double d_casca = 0.5;

	const int n_casca = ( int ) ( sqrt ( ( double ) ( nx * nx + ny * ny + nz * nz ) ) / d_casca ) + 2;

	double *soma_phi = new double[ n_casca ]();
	double *soma_n   = new double[ n_casca ]();

	const double x0 = drop.x0, y0 = drop.y0, z0 = drop.z0;

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		const int pto = geo[pos] - 1;

		const double rho_R = density ( lattice.inif_R + pto * nvel );
		const double rho_B = density ( lattice.inif_B + pto * nvel );

		const double rho = rho_R + rho_B;

		if ( rho <= eps ) continue;

		const int z = pos / ( nx * ny );
		const int y = ( pos - z * nx * ny ) / nx;
		const int x = pos - z * nx * ny - y * nx;

		//  Mesma convencao periodica das outras medidas.

		double dx = ( double ) x - x0;   if ( dx >  0.5 * nx ) dx -= nx;   if ( dx < -0.5 * nx ) dx += nx;
		double dy = ( double ) y - y0;   if ( dy >  0.5 * ny ) dy -= ny;   if ( dy < -0.5 * ny ) dy += ny;
		double dz = ( double ) z - z0;   if ( dz >  0.5 * nz ) dz -= nz;   if ( dz < -0.5 * nz ) dz += nz;

		if ( drop.cilindro ) dz = 0.0;

		//  drop.plana = true : a interface e um plano normal a x, e a coordenada do perfil e a
		//  distancia |x - x0| ao centro da laje -- o analogo unidimensional do raio.

		const double r = drop.plana ? fabs ( dx ) : sqrt ( dx * dx + dy * dy + dz * dz );

		const int k = ( int ) ( r / d_casca );

		if ( k < 0 || k >= n_casca ) continue;

		soma_phi[k] += rho_R / rho;
		soma_n  [k] += 1.0;
	}

	//  Ajuste de minimos quadrados de  ln[ phi / ( 1 - phi ) ]  contra r.

	double S = 0.0, Sr = 0.0, Sl = 0.0, Srr = 0.0, Srl = 0.0;

	int n_usadas = 0;

	for ( int k = 0; k < n_casca; k++ )
	{
		if ( soma_n[k] < 1.0 ) continue;

		const double phi = soma_phi[k] / soma_n[k];

		if ( phi <= 0.02 || phi >= 0.98 ) continue;

		const double r = ( ( double ) k + 0.5 ) * d_casca;

		const double l = log ( phi / ( 1.0 - phi ) );

		S += 1.0;  Sr += r;  Sl += l;  Srr += r * r;  Srl += r * l;

		n_usadas++;
	}

	double K = 0.0, r_meio = 0.0, res_max = 0.0;

	if ( n_usadas >= 3 )
	{
		const double den = S * Srr - Sr * Sr;

		if ( fabs ( den ) > eps )
		{
			const double a = ( S * Srl - Sr * Sl ) / den;		// inclinacao
			const double b = ( Sl - a * Sr ) / S;				// intercepto

			K = -a;

			if ( fabs ( a ) > eps ) r_meio = -b / a;			// onde ln[...] = 0, isto e phi = 1/2

			//  Residuo: quanto a logistica deixa de explicar.

			for ( int k = 0; k < n_casca; k++ )
			{
				if ( soma_n[k] < 1.0 ) continue;

				const double phi = soma_phi[k] / soma_n[k];

				if ( phi <= 0.02 || phi >= 0.98 ) continue;

				const double r = ( ( double ) k + 0.5 ) * d_casca;

				const double d = log ( phi / ( 1.0 - phi ) ) - ( a * r + b );

				if ( fabs ( d ) > res_max ) res_max = fabs ( d );
			}
		}
	}

	delete[] soma_phi;
	delete[] soma_n;

	drop.K_perfil   = K;
	drop.r_perfil   = r_meio;

	if ( drop.plana ) drop.x_interface = drop.x0 + r_meio;
	drop.res_perfil = res_max;
	drop.esp_perfil = ( K > eps ) ? 2.0 * log ( 9.0 ) / K : 0.0;
	drop.n_perfil   = n_usadas;

	return K;
}

//====================================================================================================================//


double calc_laplace_tension ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	calc_drop_radius ( geometry, lattice, drop );

	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	int n_voxels = nx * ny * nz;

	double x0 = drop.x0, y0 = drop.y0, z0 = drop.z0;

	double raio = drop.raio;

	//--------------- Limites das duas regioes -------------------------------------------------------//

	double r_in_max = drop.frac_in * raio;

	double r_out_min = drop.frac_out * raio;

	if ( r_out_min < raio + 4.0 * drop.espessura ) r_out_min = raio + 4.0 * drop.espessura;

	double r_caixa = 0.45 * ( double ) ( ( nx < ny ) ? nx : ny );

	if ( r_out_min > r_caixa ) r_out_min = r_caixa;

	//--------------- Varredura ----------------------------------------------------------------------//

	double p_in_s = 0.0, p_out_s = 0.0, p_in_s2 = 0.0, p_out_s2 = 0.0;
	double p_ctr = 0.0, p_cnt = 0.0;

	int n_in = 0, n_out = 0, n_ctr = 0;

	#pragma omp parallel for reduction(+:p_in_s,p_out_s,p_in_s2,p_out_s2,n_in,n_out,p_ctr,p_cnt,n_ctr)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		int pto = geo[pos] - 1;

		double rho_R = density ( lattice.inif_R + pto * nvel );
		double rho_B = density ( lattice.inif_B + pto * nvel );

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) ) continue;

		int z = pos / ( nx * ny );
		int y = ( pos - z * nx * ny ) / nx;
		int x = pos - z * nx * ny - y * nx;

		double dx = ( double ) x - x0;   if ( dx >  0.5 * nx ) dx -= nx;   if ( dx < -0.5 * nx ) dx += nx;
		double dy = ( double ) y - y0;   if ( dy >  0.5 * ny ) dy -= ny;   if ( dy < -0.5 * ny ) dy += ny;
		double dz = ( double ) z - z0;   if ( dz >  0.5 * nz ) dz -= nz;   if ( dz < -0.5 * nz ) dz += nz;

		if ( drop.cilindro ) dz = 0.0;

		double r = sqrt ( dx * dx + dy * dy + dz * dz );

		//  Pressao do sitio.  Ver DROP::cs2_R / cs2_B: modelos de alpha variavel dao a cada fluido
		//  a sua velocidade do som, e ai a pressao nao e c_s^2 rho.

		double p = ( drop.cs2_R >= 0.0 ) ? ( rho_R * drop.cs2_R + rho_B * drop.cs2_B )
		                                 : ( lattice.c_s2 * ( rho_R + rho_B ) );

		if      ( r <  r_in_max  ) { p_in_s  = p_in_s  + p;   p_in_s2  = p_in_s2  + p * p;   n_in  = n_in  + 1; }

		else if ( r >= r_out_min ) { p_out_s = p_out_s + p;   p_out_s2 = p_out_s2 + p * p;   n_out = n_out + 1; }

		if ( r < 1.5 )   { p_ctr = p_ctr + p;   n_ctr = n_ctr + 1; }

		if ( pos == 0 )    p_cnt = p_cnt + p;
	}

	//--------------- Devolve --------------------------------------------------------------------------//

	double p_in  = ( n_in  > 0 ) ? p_in_s  / ( double ) n_in  : 0.0;
	double p_out = ( n_out > 0 ) ? p_out_s / ( double ) n_out : 0.0;

	double var_in  = ( n_in  > 0 ) ? p_in_s2  / ( double ) n_in  - p_in  * p_in  : 0.0;
	double var_out = ( n_out > 0 ) ? p_out_s2 / ( double ) n_out - p_out * p_out : 0.0;

	double fator_forma = drop.cilindro ? 1.0 : 0.5;

	drop.p_in  = p_in;    drop.sd_in  = ( var_in  > 0.0 ) ? sqrt ( var_in  ) : 0.0;
	drop.p_out = p_out;   drop.sd_out = ( var_out > 0.0 ) ? sqrt ( var_out ) : 0.0;

	drop.n_in = n_in;   drop.n_out = n_out;

	drop.delta_p = p_in - p_out;

	drop.sigma = fator_forma * raio * drop.delta_p;

	drop.p_centro = ( n_ctr > 0 ) ? p_ctr / ( double ) n_ctr : 0.0;

	drop.p_canto = p_cnt;

	drop.sigma_linha = fator_forma * drop.raio_linha * ( drop.p_centro - drop.p_canto );

	return drop.sigma;
}

//====================================================================================================================//


//=============================== Modal amplitude of a capillary wave ================================================//
//
//      Localiza a interface em cada coluna y -- o x onde  rho^N = ( rho_R - rho_B ) / rho  cruza
//      zero, por interpolacao linear entre os dois sitios vizinhos -- e projeta o deslocamento
//      resultante no modo excitado:
//
//          delta( y ) = x_int( y ) - x_med
//
//          a_cos = ( 2 / ny ) soma_y delta( y ) cos( k y )
//          a_sin = - ( 2 / ny ) soma_y delta( y ) sin( k y )
//
//          amp = sqrt( a_cos^2 + a_sin^2 )        fase = atan2( a_sin, a_cos )
//
//      de modo que  delta( y ) = amp cos( k y + fase )  reproduz o modo.  A projecao usa as ny
//      colunas: o erro de leitura de cada uma entra dividido por sqrt( ny ).
//
//      'amp_res' e a potencia que sobra fora do modo excitado ( harmonicos gerados pela nao
//      linearidade, mais ruido ).  Enquanto amp_res << amp a onda esta no regime linear.
//
//      A varredura em x vai de fora a fora e pega o PRIMEIRO cruzamento de + para - a partir do
//      lado vermelho, que e a unica interface do ensaio.  Colunas sem cruzamento entram em
//      'n_faltou' em vez de contaminar a media: se n_faltou > 0 a interface se perdeu em algum
//      ponto e a medida nao vale.
//
//      Input : geometry, lattice, wave ( modo, cs2_R, cs2_B )
//      Output: wave.k, x_med, a_cos, a_sin, amp, fase, amp_res, n_linhas, n_faltou, massas,
//              qx/qy/qz, u_max, u_rms, n_fluid, n_ruim.  Devolve wave.amp.
//
//====================================================================================================================//

double calc_onda_capilar ( GEOMETRY geometry, LATTICE lattice, WAVE& onda )
{
	int *geo = geometry.ini;

	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	const double PI = 3.14159265358979323846;

	const double k = 2.0 * PI * ( double ) onda.modo / ( double ) ny;

	onda.k = k;

	//--------------- Passo 1: campo, massas e correntes ---------------------------------------------//

	double massa_R = 0.0, massa_B = 0.0;
	double qx = 0.0, qy = 0.0, qz = 0.0;
	double u2_sum = 0.0, u_max = 0.0;

	int n_fluid = 0, n_ruim = 0;

	//  rho^N de cada sitio, guardado para a busca da interface.  nullptr onde nao ha fluido.

	const long long n_voxels = ( long long ) nx * ny * nz;

	vector<double> fase_no ( ( size_t ) n_voxels, 2.0 );		// 2.0 marca "nao ha fluido aqui"

	#pragma omp parallel for \
		reduction(+:massa_R,massa_B,qx,qy,qz,u2_sum,n_fluid,n_ruim) reduction(max:u_max)

	for ( long long pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		const int pto = geo[pos] - 1;

		double vxR, vyR, vzR, rho_R;
		double vxB, vyB, vzB, rho_B;

		calcula ( lattice.inif_R + pto * nvel, vxR, vyR, vzR, rho_R, lattice );
		calcula ( lattice.inif_B + pto * nvel, vxB, vyB, vzB, rho_B, lattice );

		n_fluid = n_fluid + 1;

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) || ! ( vxR == vxR ) || ! ( vxB == vxB ) )
		{
			n_ruim = n_ruim + 1;

			continue;
		}

		massa_R = massa_R + rho_R;
		massa_B = massa_B + rho_B;

		qx = qx + rho_R * vxR + rho_B * vxB;
		qy = qy + rho_R * vyR + rho_B * vyB;
		qz = qz + rho_R * vzR + rho_B * vzB;

		const double rho = rho_R + rho_B;

		if ( rho <= 0.0 ) continue;

		fase_no[ ( size_t ) pos ] = ( rho_R - rho_B ) / rho;

		double ux = ( rho_R * vxR + rho_B * vxB ) / rho;
		double uy = ( rho_R * vyR + rho_B * vyB ) / rho;
		double uz = ( rho_R * vzR + rho_B * vzB ) / rho;

		//  Correcao da meia-forca, a mesma de calc_drop_radius: depois da colisao o momento cru
		//  traz - F / 2.  Sem isto a "velocidade" da interface sai contaminada pela forca capilar.

		if ( lattice.ini_force != nullptr )
		{
			const double *F = lattice.ini_force + ( long long ) pto * 3;

			ux = ux - 0.5 * F[0] / rho;
			uy = uy - 0.5 * F[1] / rho;
			uz = uz - 0.5 * F[2] / rho;

			qx = qx - 0.5 * F[0];
			qy = qy - 0.5 * F[1];
			qz = qz - 0.5 * F[2];
		}

		const double u2 = ux * ux + uy * uy + uz * uz;

		u2_sum = u2_sum + u2;

		const double u = sqrt ( u2 );

		if ( u > u_max ) u_max = u;
	}

	onda.massa_R = massa_R;   onda.massa_B = massa_B;
	onda.qx = qx;   onda.qy = qy;   onda.qz = qz;
	onda.n_fluid = n_fluid;   onda.n_ruim = n_ruim;
	onda.u_max = u_max;
	onda.u_rms = ( n_fluid > 0 ) ? sqrt ( u2_sum / ( double ) n_fluid ) : 0.0;

	//--------------- Passo 2: posicao da interface em cada coluna -----------------------------------//
	//
	//  Media sobre z quando a caixa tem mais de uma camada; o ensaio padrao usa nz = 1.

	vector<double> h ( ( size_t ) ny, 0.0 );
	vector<int>    n_h ( ( size_t ) ny, 0 );

	int n_faltou = 0;

	for ( int y = 0; y < ny; y++ )
	{
		for ( int z = 0; z < nz; z++ )
		{
			bool achou = false;

			for ( int x = 0; x + 1 < nx; x++ )
			{
				const double a = fase_no[ ( size_t ) ( x     + y * nx + ( long long ) z * nx * ny ) ];
				const double b = fase_no[ ( size_t ) ( x + 1 + y * nx + ( long long ) z * nx * ny ) ];

				if ( a > 1.5 || b > 1.5 ) continue;		// sitio solido

				//  Cruzamento de vermelho ( rho^N > 0 ) para azul ( rho^N < 0 ).

				if ( a >= 0.0 && b < 0.0 )
				{
					const double t = a / ( a - b );		// a - b > 0 aqui

					h[ ( size_t ) y ] += ( double ) x + t;

					n_h[ ( size_t ) y ] += 1;

					achou = true;

					break;
				}
			}

			if ( ! achou ) n_faltou++;
		}
	}

	//--------------- Passo 3: projecao de Fourier ---------------------------------------------------//

	double soma = 0.0;

	int n_linhas = 0;

	for ( int y = 0; y < ny; y++ )
	{
		if ( n_h[ ( size_t ) y ] == 0 ) continue;

		h[ ( size_t ) y ] /= ( double ) n_h[ ( size_t ) y ];

		soma += h[ ( size_t ) y ];

		n_linhas++;
	}

	onda.n_linhas = n_linhas;
	onda.n_faltou = n_faltou;

	if ( n_linhas == 0 )
	{
		onda.x_med = 0.0;   onda.a_cos = 0.0;   onda.a_sin = 0.0;
		onda.amp = 0.0;     onda.fase = 0.0;    onda.amp_res = 0.0;

		return 0.0;
	}

	const double x_med = soma / ( double ) n_linhas;

	double a_cos = 0.0, a_sin = 0.0, pot = 0.0;

	for ( int y = 0; y < ny; y++ )
	{
		if ( n_h[ ( size_t ) y ] == 0 ) continue;

		const double d = h[ ( size_t ) y ] - x_med;

		a_cos += d * cos ( k * ( double ) y );
		a_sin -= d * sin ( k * ( double ) y );

		pot += d * d;
	}

	a_cos *= 2.0 / ( double ) n_linhas;
	a_sin *= 2.0 / ( double ) n_linhas;

	const double amp = sqrt ( a_cos * a_cos + a_sin * a_sin );

	//  Potencia total do deslocamento: para um cosseno puro de amplitude A, media( d^2 ) = A^2 / 2.

	const double resto = 2.0 * pot / ( double ) n_linhas - amp * amp;

	onda.x_med = x_med;
	onda.a_cos = a_cos;
	onda.a_sin = a_sin;
	onda.amp = amp;
	onda.fase = atan2 ( a_sin, a_cos );
	onda.amp_res = ( resto > 0.0 ) ? sqrt ( resto ) : 0.0;

	return onda.amp;
}

//====================================================================================================================//

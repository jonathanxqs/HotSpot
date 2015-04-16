#ifndef _DEF_H
#include "def.h"
#endif

void raytracing(double spin, double epsilon, 
	            double orbit_omega, double r_orbit, 
	            double r_spot, double period, double hstart,
				double tt, double robs, double rstep2, 
				double pobs, double pstep, double iobs, 
				double dobs, double aatol, double rtol, 
				double E_obs[N_E], double fphi[N_T][N_E],
				double fphi_1[N_T][N_E],
				vector<string>  ofnames, int &hit_sum)
{	
	double spin2 = spin*spin;
	int stop_integration;
	int i, j;
	int hitcheck;
	int n_period;
	double xobs, yobs;
	double r0, r02, s0, s02;
	double fact1, fact2, fact3, fact4;
	double t0, x0, y0, phi0;
	double kt0, kx0, ky0, kphi0;
	double tau, xau, yau, zau, kyau;
	double hnext;
	double const0, const1;
	double t, x, y, phi;
	double horizon;
	double pp, qq;
	double xem[4] = {0, 0, 0, 0};
	double gfactor;
	double limbdark;
	double v[4], p[4];
	double Dspot;
	double tstep;
	double n_t_period;
	double *pointer; 
	double *pointer_1; 

	double Estep = double(G_MAX)*1.0/N_E;

	Dspot = 0;  

	xobs = robs*cos(pobs);
	yobs = robs*sin(pobs);
							
	/* ----- compute photon initial conditions ----- */		
	r02 = xobs*xobs + yobs*yobs + dobs*dobs;
	r0  = sqrt(r02);
			
	fact1 = dobs*sin(iobs) - yobs*cos(iobs);
	fact2 = fact1*fact1;
	fact3 = xobs*xobs + fact2;
	fact4 = sqrt(fact3);
				
	t0   = tt;
	x0   = r0;
	y0   = acos((yobs*sin(iobs) + dobs*cos(iobs))/r0);
	phi0 = atan(xobs/fact1);
				
	s0  = sin(y0);
	s02 = s0*s0;
				
	kx0   = - dobs/r0;
	ky0   = (cos(iobs) - dobs*(yobs*sin(iobs) + dobs*cos(iobs))/r02)/fact4;
	kphi0 = xobs*sin(iobs)/fact3;
	kt0   = sqrt(kx0*kx0 + r02*ky0*ky0 + r02*s02*kphi0*kphi0);
										
	/* ----- solve geodesic equations 
	fourth-order runge-kutta-nystrom method 
	see E. Lund et al., 2009 JINST 4 P04001 ----- */							
	v[0] = t0;
	v[1] = x0;
	v[2] = y0;
	v[3] = phi0;
					
	p[0] = kt0;
	p[1] = kx0;
	p[2] = ky0;
	p[3] = kphi0;

	tau = 0;
	xau = 0;
	yau = 0;
	zau = 0;
	kyau = 0;
							
	const0 = kt0;
	const1 = r02*s02*kphi0/kt0;			
	stop_integration = 0;
			
	hnext = hstart;			
	do {							
		RungeKuttaNystrom(spin, epsilon, hnext, 
			v, p, tau, xau, yau, zau, kyau, aatol, rtol);															
	
		t = v[0];
		x = v[1];
		y = v[2];
		phi = v[3];	

		if (y > Pi/2) {
			intersection_new(tau, xau, yau, zau, t, 
				x, y, phi, xem);
			stop_integration = 1;				 				
		}
						
		horizon = horizon_JP(spin, epsilon, x, y);
		if (horizon < 0.05) {
			stop_integration = 4;   /* the photon crosses the horizon */
		} else if (x < 1) {
			stop_integration = 5;   /* the photon hits the singularity */
		} else if (x != x) {
			stop_integration = 6;   /* numerical problems! */
		} else if (t < 0) {
			stop_integration = 7;   /* numerical problems! */
		} else if (x > 1.05*dobs) {
			stop_integration = 8;   /* the photon escapes to infinity */
		}							
	
	} while (stop_integration == 0);

	int temp = N_T/N_gif;
	vector<string>::iterator it = ofnames.begin();
					
	if (stop_integration == 1) {
		tstep = (double) T_MAX/N_T;
		n_period = int (T_MAX/period);
		n_t_period = n_period*period/tstep;

		for (i = 0; i < n_t_period; xem[0] = xem[0] - tstep) {
			spot(orbit_omega, r_orbit, r_spot, period,
	             xem, Dspot, hitcheck);

			if (hitcheck == 1) {
				++hit_sum;
				redshift(spin, orbit_omega, epsilon, xem[1], const0, 
					const1, kyau, gfactor, limbdark);
		        /* Upsilon = 1 for isotropic radiation; 
			    Upsilon = limbdark[0] for limb-darkened radiation */			
		        /*Upsilon = 1;			 
			    Upsilon = limbdark;*/

			    /* --- integration - part 1 --- */
		        pp = gfactor;
		        qq = gfactor*gfactor*gfactor*gfactor;
		        qq = qq*exp(- Dspot*Dspot/(2*r_spot*r_spot));
		        qq = qq*robs*robs*rstep2*pstep*cos(iobs);

		        pointer = fphi[i];
				pointer_1 = fphi_1[i];
			    /* if (i % temp == 0){
					int shift = i / temp;
				    ofstream outfile((*(it+shift)).c_str(), ofstream::app);
					outfile<<setiosflags(ios::fixed)<<i<<" \t "
						<<setprecision(10)<<xobs<<" \t "<<setprecision(10)
						<<yobs<<" \t "<<setprecision(10)<<gfactor<<" \t "
						<<setprecision(10)<<qq<<"\n";
				} 
				*/

				//for (j = 0; j < N_E; j++) {			
				//	if (E_obs[j] < pp && E_obs[j + 1] > pp) {
				//		*(pointer + j) += qq;
				//		*(pointer_1 + j) += qq;
				//		j = N_E;
				//	}							
				//}

				j = (int) floor(pp/Estep);
				j = j < N_E ? j : (N_E - 1);
				*(pointer + j) += qq;
				*(pointer_1 + j) += qq;

			}
			
			++i;
		}
	} 

	//return;

	int cross_check = 0;
	stop_integration = 0;
	do {							
		RungeKuttaNystrom(spin, epsilon, hnext, 
			v, p, tau, xau, yau, zau, kyau, aatol, rtol);															
	
		t = v[0];
		x = v[1];
		y = v[2];
		phi = v[3];	

		if (y > Pi/2 && cross_check == 0) {
			cross_check = 1;
		}
		if (cross_check == 1 && y < Pi/2) {
			 intersection_new(tau, xau, yau, zau, t, 
				 x, y, phi, xem);
			 stop_integration = 1;
			 cross_check = 0;
		}
						
		horizon = horizon_JP(spin, epsilon, x, y);
		if (horizon < 0.05) {
			stop_integration = 4;   /* the photon crosses the horizon */
		} else if (x < 1) {
			stop_integration = 5;   /* the photon hits the singularity */
		} else if (x != x) {
			stop_integration = 6;   /* numerical problems! */
		} else if (t < 0) {
			stop_integration = 7;   /* numerical problems! */
		} else if (x > 1.05*dobs) {
			stop_integration = 8;   /* the photon escapes to infinity */
		}							
	
	} while (stop_integration == 0);

	temp = N_T/N_gif;
	it = ofnames.begin();

	if (stop_integration == 1) {
		tstep = (double) T_MAX/N_T;
		n_period = int (T_MAX/period);
		n_t_period = n_period*period/tstep;

		for (i = 0; i < n_t_period; xem[0] = xem[0] - tstep) {
			spot(orbit_omega, r_orbit, r_spot, period,
	             xem, Dspot, hitcheck);

			if (hitcheck == 1) {
				++hit_sum;
				redshift(spin, orbit_omega, epsilon, xem[1], const0, 
					const1, kyau, gfactor, limbdark);
		        /* Upsilon = 1 for isotropic radiation; 
			    Upsilon = limbdark[0] for limb-darkened radiation */			
		        /*Upsilon = 1;			 
			    Upsilon = limbdark;*/

			    /* --- integration - part 1 --- */
		        pp = gfactor;
		        qq = gfactor*gfactor*gfactor*gfactor;
		        qq = qq*exp(- Dspot*Dspot/(2*r_spot*r_spot));
		        qq = qq*robs*robs*rstep2*pstep*cos(iobs);

		        pointer = fphi[i];
			    /* if (i % temp == 0){
					int shift = i / temp;
				    ofstream outfile((*(it+shift)).c_str(), ofstream::app);
					outfile<<setiosflags(ios::fixed)<<i<<" \t "
						<<setprecision(10)<<xobs<<" \t "<<setprecision(10)
						<<yobs<<" \t "<<setprecision(10)<<gfactor<<" \t "
						<<setprecision(10)<<qq<<"\n";
				} 
				*/

				//for (j = 0; j < N_E; j++) {			
				//	if (E_obs[j] < pp && E_obs[j + 1] > pp) {
				//		*(pointer + j) += qq;
				//		j = N_E;
				//	}							
				//}	

				j = (int) floor(pp/Estep);
				j = j < N_E ? j : (N_E - 1);
				*(pointer + j) += qq;

			}
			
			++i;
		}
	} 

	return;
}

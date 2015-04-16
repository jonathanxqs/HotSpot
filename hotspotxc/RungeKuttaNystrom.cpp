/* ----- solve geodesic equations using fourth-order 
runge-kutta-nystrom method, see E. Lund et al., 2009 
JINST 4 P04001 for details----- */

#ifndef _DEF_H
#include "def.h"
#endif

void RungeKuttaNystrom(double spin, double epsilon, 
	                   double &hnext, double v[4], double p[4],
					   double &tau, double &xau, double &yau, 
					   double &zau, double &kyau, 
					   double aatol, double rtol)
{
	double spin2 = spin*spin;
	int i, n1, n2, n3;
	double h;
	double v1, v2;
	double check;
	double u[4];
	double RK1[4], RK2[4], RK3[4], RK4[4];
	double Gamma[4][4][4];
	double verr[4], vtol[4], err[4];
					
	do {
		h = hnext;

		/* ----- compute RK1 ----- */				
		v1 = v[1];						
		v2 = v[2];		
		Christoffel(spin, epsilon, v1, v2, Gamma);
		for (i = 0; i <= 3; i++) {
			u[i] = p[i];
		}

		for (n1 = 0; n1 <= 3; n1++) {		
			for (n2 = 0; n2 <= 3; n2++) {				
				for (n3 = 0; n3 <= 3; n3++) {				
					if (n2 == 0 && n3 == 0) {				
						RK1[n1] = - Gamma[n1][0][0]*u[0]*u[0];
					} else {
						RK1[n1] -= Gamma[n1][n2][n3]*u[n2]*u[n3];
					}					
				}							
			}						
		}											
						
		/* ----- compute RK2 ----- */										
		v1 = v[1] + h*p[1]/2 + h*h*RK1[1]/16;						
		v2 = v[2] + h*p[2]/2 + h*h*RK1[2]/16;											
		Christoffel(spin, epsilon, v1, v2, Gamma);										
		for (i = 0; i <= 3; i++) {
			u[i] = p[i] + h*RK1[i]/4;
		}
											
		for (n1 = 0; n1 <= 3; n1++) {				
			for (n2 = 0; n2 <= 3; n2++) {				
				for (n3 = 0; n3 <= 3; n3++) {				
					if (n2 == 0 && n3 == 0) {
						RK2[n1] = - Gamma[n1][0][0]*u[0]*u[0];
					} else {
						RK2[n1] -= Gamma[n1][n2][n3]*u[n2]*u[n3];
					}		
				}							
			}					
		}
															
		/* ----- compute RK3 ----- */			
		for (i = 0; i <= 3; i++) {
			u[i] = p[i] + h*RK2[i]/4;
		}
						
		for (n1 = 0; n1 <= 3; n1++) {				
			for (n2 = 0; n2 <= 3; n2++) {				
				for (n3 = 0; n3 <= 3; n3++) {				
					if (n2 == 0 && n3 == 0) {
						RK3[n1] = - Gamma[n1][0][0]*u[0]*u[0];
					} else {
						RK3[n1] -= Gamma[n1][n2][n3]*u[n2]*u[n3];
					}					
				}							
			}						
		}					
					
		/* ----- compute RK4 ----- */				
		v1 = v[1] + h*p[1] + h*h*RK3[1]/4;						
		v2 = v[2] + h*p[2] + h*h*RK3[2]/4;										
		Christoffel(spin, epsilon, v1, v2, Gamma);								
		for (i = 0; i <= 3; i++) {
			u[i] = p[i] + h*RK3[i]/2;
		}
					
		for (n1 = 0; n1 <= 3; n1++) {				
			for (n2 = 0; n2 <= 3; n2++) {		
				for (n3 = 0; n3 <= 3; n3++) {			
					if (n2 == 0 && n3 == 0) {			
						RK4[n1] = - Gamma[n1][0][0]*u[0]*u[0];
					} else {
						RK4[n1] -= Gamma[n1][n2][n3]*u[n2]*u[n3];			
					}							
				}						
			}					
		}
																
		/* ----- local error ----- */									
		for (i = 0; i <= 3; i++) {				
			verr[i] = 0.5*h*h*(RK1[i] - RK2[i] - RK3[i] + RK4[i]);							
			verr[i] *= verr[i];							
			vtol[i] = aatol + fabs(v[i])*rtol;							
			vtol[i] *= vtol[i];				
			err[i] = 0.25*verr[i]/vtol[i]; 			
		}
										
		check = sqrt(err[0] + err[1] + err[2] + err[3]);															
		/* ----- next step ----- */										
		hnext = h*pow(1/check,0.25);
															
		/* ----- limitation criterion ----- */				
		if (hnext < h/4) { 
			hnext = h/4;
		} else if (hnext > 4*h) {
			hnext = 4*h;
		} else if (hnext == h) {
			hnext = 0.9*h;
		}

	} while (check > 1);
											
	/* ----- solutions to the fourth-order RKN method ----- */	
	tau = v[0];
	xau = v[1];					
	yau = v[2];					
	zau = v[3];									
	kyau = p[2];
					
	for (i = 0; i <= 3; i++)  { 
		v[i] += h*p[i] + (RK1[i] + RK2[i] + RK3[i])*h*h/12;								
		p[i] += (RK1[i] + 2*RK2[i] + 2*RK3[i] + RK4[i])*h/12;
	}

	return;
}

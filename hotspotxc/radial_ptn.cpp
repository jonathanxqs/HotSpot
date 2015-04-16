/*The function radialpotential is used to calculate specific
energy E, specific angular momentum L, \partial_{rr}{V_{eff}} 
and \partial_{zz}{V_{eff}} in equatorial plane with given 
spin parameter, deformation parameter b and radius xx.*/
#ifndef _DEF_H
#include "def.h"
#endif

int radial_ptn (double spin, double gg, double r,
	                double &Vrr, double &Vzz)
{
	double spin2 = spin*spin;
	int i;
	double dr, z, dz, temp, theta;
	double x, x2;
	double g001, g031, g331;
	double aux;
	double orbit_omega;
	double gmn[3][4][4];
	double Veff[3];
	double E;
	double L;

	/* check stability along the radial direction */
	dr   = 0.001*r;
	temp = r + dr;
	dr   = temp - r;
	
	for (i = 0; i <= 2; i++) {
		x  = r + 0.5*dr*(i - 1);
		metric_JP(spin, gg, x, Pi/2, gmn[i]);	
	}
	
	g001 = 0.5*(gmn[2][0][0] - gmn[0][0][0])/dr;
	g031 = 0.5*(gmn[2][0][3] - gmn[0][0][3])/dr;
	g331 = 0.5*(gmn[2][3][3] - gmn[0][3][3])/dr;
	orbit_omega  = (-g031 + sqrt(g031*g031 - g001*g331))/g331;
	aux = sqrt(-gmn[1][0][0] -2*gmn[1][0][3]*orbit_omega - gmn[1][3][3]*orbit_omega*orbit_omega);
	
	E = - (gmn[1][0][0] + gmn[1][0][3]*orbit_omega)/aux;
	L = (gmn[1][0][3] + gmn[1][3][3]*orbit_omega)/aux;
	
	for (i = 0; i <= 2; i++) {
		aux = E*E*gmn[i][3][3] + 2*E*L*gmn[i][0][3] + L*L*gmn[i][0][0];
		Veff[i] = aux/(gmn[i][0][3]*gmn[i][0][3] - gmn[i][0][0]*gmn[i][3][3]);	
	}
	
	Vrr = Veff[2] - 2*Veff[1] + Veff[0];
	
	/* check stability along the vertical direction */
	x = r;
	
	for (i = 0; i <= 2; i++) {
		dz = 0.001;
		z = 0.5*dz*(i - 1);
		x2 = x*x + z*z;
		x = sqrt(x2);
		theta = acos(z/x);
		
		metric_JP(spin, gg, x, theta, gmn[i]);
        aux = E*E*gmn[i][3][3] + 2*E*L*gmn[i][0][3] + L*L*gmn[i][0][0];
		Veff[i] = aux/(gmn[i][0][3]*gmn[i][0][3] - gmn[i][0][0]*gmn[i][3][3]);	
	}
	
	Vzz = Veff[2] - 2*Veff[1] + Veff[0];

	return 0;	
}

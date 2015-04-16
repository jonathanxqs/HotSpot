#ifndef _DEF_H
#include "def.h"
#endif

void redshift(double spin, double isco_omega, double epsilon, 
			  double radius, double ktt, double ktkp, 
			  double kyy, double& gg, double& ldr)
{
	int i;
	double x;
	double dr, temp;
	double cc;
	double gmn[3][4][4];
	double g001, g031, g331;
	double omega;
	double uet;
	double uephi;
	double mem;
	
	dr = 0.00001*radius;
	temp = radius + dr;
	dr = temp - radius;
	cc = ktkp;

	for (i = 0; i <= 2; i++) {
		x  = radius + 0.5*dr*(i - 1);
		metric_JP(spin, epsilon, x, Pi/2, gmn[i]);	
	}
	
	g001 = 0.5*(gmn[2][0][0] - gmn[0][0][0])/dr;
	g031 = 0.5*(gmn[2][0][3] - gmn[0][0][3])/dr;
	g331 = 0.5*(gmn[2][3][3] - gmn[0][3][3])/dr;
	omega  = (-g031 + sqrt(g031*g031 - g001*g331))/g331;
	uet = sqrt(-gmn[1][0][0] -2*gmn[1][0][3]*isco_omega - gmn[1][3][3]*isco_omega*isco_omega);

	gg = uet/(1 - cc*isco_omega);
	uet = 1/uet;
	uephi = omega*uet;
	mem = kyy*radius;
	mem = mem/(ktt*uet + ktt*cc*uephi);
	
	if (mem < 0) {
		mem = 0;
	} else if (mem > 1) {
		mem = 1;
	}

	ldr = 0.5 + 0.75*mem;

	return ;
}

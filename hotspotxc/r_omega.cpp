#ifndef _DEF_H
#include "def.h"
#endif

double r_omega(double spin, double epsilon, double radius)
{
	double spin2 = spin*spin;
	int i;
	double dr, temp;
	double x;
	double omega;
	double g001, g031, g331;
	double gmn[3][4][4];

	dr   = 0.00001*radius;
	temp = radius + dr;
	dr   = temp - radius;
	
	for (i = 0; i <= 2; i++) {
		x  = radius + 0.5*dr*(i - 1);
		metric_JP(spin, epsilon, x, Pi/2, gmn[i]);	
	}
	
	g001 = 0.5*(gmn[2][0][0] - gmn[0][0][0])/dr;
	g031 = 0.5*(gmn[2][0][3] - gmn[0][0][3])/dr;
	g331 = 0.5*(gmn[2][3][3] - gmn[0][3][3])/dr;
	omega  = (-g031 + sqrt(g031*g031 - g001*g331))/g331;
	
	return omega;	
}

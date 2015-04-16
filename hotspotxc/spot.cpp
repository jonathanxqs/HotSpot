#ifndef _DEF_H
#include "def.h"
#endif

void spot(double isco_omega, double r_orbit, 
	      double r_spot, double period,
	      double xem[4], double &D, int &check)
{
	double tt;
	double x_1, y_1, r_1, phi_1, x_2, y_2, r_2, phi_2;

	tt = fmod(xem[0], period);
	
	r_1 = r_orbit;
	phi_1 = 2*Pi*tt/period;
	r_2 = xem[1];
	phi_2 = xem[3];

	x_1 = r_1*cos(phi_1);
	y_1 = r_1*sin(phi_1);
	x_2 = r_2*cos(phi_2);
	y_2 = r_2*sin(phi_2);

	D = sqrt((x_1-x_2)*(x_1-x_2) + (y_1-y_2)*(y_1-y_2));

	check = 2;
	if (D < 4*r_spot)  {	
		check = 1;
	}

	return;
}

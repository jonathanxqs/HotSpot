#ifndef _DEF_H
#include "def.h"
#endif

double find_isco(double spin, double gg)
{
	int check;
	double r;
	double Vrr = 0, Vzz = 0;
	double r_isco;
	double r_strt, r_step, r_accy;
	double r_min, r_mid, r_max, r_dlt;

	r_strt = 20;
	r_step = 0.5;
	r_accy = 0.0000001;

	r = r_strt;
	check = 0;
	while (check == 0)  {
		radial_ptn(spin, gg, r, Vrr, Vzz);

		if (Vrr < 0 && Vzz < 0) {
			r -= r_step;
			check = 0;
		} else if (Vrr > 0) {
			check = 1;
		} else if (Vzz > 0) {
			check = 2;
		} else if (r <= r_accy) {
			check = 3;
		} else {
			check = -1;
		}
	}

	r_min = r;
	r_max = r + r_step;
	r_dlt = (r_max - r_min)/2;

	while (r_dlt > r_accy) {
		r_mid = (r_min + r_max)/2;
		radial_ptn(spin, gg, r_mid, Vrr, Vzz);
		r_isco = r_max;

		if (Vrr < 0 && Vzz < 0) {
			r_max = r_mid;
		} else {
			r_min = r_mid;
		}
		r_dlt = (r_max - r_min)/2;
	}

	r_isco = r_max;		
	if (r_isco > r_strt - 0.001) {
		ofstream  oferr("error-isco.dat");
		oferr << "r_isco > r_start for spin = " << spin
			  << " and gg = " << gg;  
		oferr.close();
	}

	return r_isco;	
}

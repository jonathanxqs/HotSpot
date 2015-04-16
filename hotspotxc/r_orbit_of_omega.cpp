#ifndef _DEF_H
#include "def.h"
#endif

double r_orbit_of_omega(double spin, double epsilon, 
						double omega, double isco)
{
	double r_accy;
	double r_min, r_mid, r_max, r_dlt;
	double omega_mid;

	r_accy = 0.0000001;

	r_min = isco;
	r_max = isco + 100;
	r_mid = (r_min + r_max)/2;
	r_dlt = (r_max - r_min)/2;

	while (r_dlt > r_accy)  {
		omega_mid = r_omega(spin, epsilon, r_mid);
		if (omega < omega_mid) {
			r_min = r_mid;
		} else {
			r_max = r_mid;
		}
		r_mid = (r_min + r_max)/2;
		r_dlt = (r_max - r_min)/2;
	}

	return r_mid/isco;
}

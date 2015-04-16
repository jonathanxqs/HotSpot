#ifndef _DEF_H
#include "def.h"
#endif

int outp_spec(double spectra[N_T][N_E], 
			  double outp_spectra[N_T_PRT][N_E_PRT]) 
{
	int E_stp_prt, T_stp_prt;
	int i, j, ii, jj;
	double tmp_spectra[N_T][N_E_PRT];

	E_stp_prt = N_E/N_E_PRT;
	T_stp_prt = N_T/N_T_PRT;

	double *tmp = tmp_spectra[0];
	double *spc = spectra[0]; 
	for (i = 0; i < N_T; i++) {
		for (j = 0; j < N_E_PRT; j++) {
			*tmp = 0;
			for (jj = 0; jj < E_stp_prt; jj++) {
				*tmp += *spc++;
			}
			tmp++;
		}
	}

	double *out = outp_spectra[0];
	tmp = tmp_spectra[0];
	for (i = 0; i < N_T_PRT; i++) {
		for (j = 0; j < N_E_PRT; j++) {
			tmp = tmp_spectra[i*T_stp_prt] + j;
			for (ii = 0; ii < T_stp_prt; ii++) {
				*out += *tmp;
				tmp += N_E_PRT;
			}
			out++;
		}
	}

	return 0;
}

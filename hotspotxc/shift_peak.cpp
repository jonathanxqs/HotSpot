#ifndef _DEF_H
#include "def.h"
#endif

int shift_peak(double prd, double T_obs[], double int_spec[], 
	           double spectra[N_T][N_E])
{
	int i, j, k;
	int i_zero, n_prd;
	int peak_post, peak_shft;
	int extd_n_t;
	double t_stp, n_tt;
	double int_intst;
	double spec_peak;
	vector <double> extd_spec;
	vector < vector<double> > extd_spectra;
	
	t_stp = (double) T_MAX/N_T;
	n_prd = (int) (T_MAX/prd);
	n_tt = n_prd*prd/t_stp;
	i_zero = (int) floor(n_tt);

	int_intst = integral(N_T, T_obs, int_spec);
	for (i = 0; i < N_T; i++) {
    	int_spec[i] = n_prd*prd*int_spec[i]/int_intst;
	}

	j = 0;
	extd_n_t = (N_T/(i_zero + 1) + 1)*(i_zero + 1);
	for (i = 0; i < extd_n_t; ++i)  {
		j = i % (i_zero + 1);
		extd_spec.push_back(int_spec[j]);
		vector <double> tmp_spectra(spectra[j], spectra[j] + N_E);
		extd_spectra.push_back(tmp_spectra);
	}

	spec_peak = 0;
	peak_post = 0;
	for (i = 0; i < extd_n_t; i++) {
    	if (extd_spec[i] >= spec_peak) {
    		spec_peak = extd_spec[i];
    		peak_post = i;
    	}
	}
	peak_shft = N_T/2 - peak_post;
	for (i = 0; i < N_T; i++) {
		j = ( - peak_shft + i + extd_n_t) % extd_n_t;
    	int_spec[i] = extd_spec[j];

    	k = 0;
    	vector<double>::iterator iter = extd_spectra[j].begin();
    	for (; iter != extd_spectra[j].end(); iter++)  {
    		*(spectra[i] + k) = *iter;
    		++k;
    	}
	}

	return 0;
}

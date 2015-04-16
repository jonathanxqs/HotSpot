#ifndef _DEF_H

#define _DEF_H

#define T_MAX 100
#define G_MAX 4

#define R_N 1000
#define PHI_N 1000

#define N_E 400
#define N_T 400

#define N_E_PRT 400
#define N_T_PRT 400

#define N_gif 50

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <time.h>

using namespace std;

const double Pi  = 3.141592653589793;

void Christoffel(double spin, double epsilon, double w1, 
	             double w2, double CS[][4][4]);

void redshift(double spin, double orbit_omega, double epsilon,
			  double R_ORBIT, double ktt, double ktkp, 
			  double kyy, double& gg, double& ldr);

void find_isco_old(double spin, double spin2, 
				   double epsilon, double& r_isco);

double find_isco_old(double spin, double epsilon, 
	             double &orbit_omega);

double find_isco(double spin, double epsilon);

void metric_JP(double spin, double epsilon, double r, 
	           double theta, double mn[4][4]);

void RungeKuttaNystrom(double spin, double epsilon, 
	                   double &hnext, double v[4], double p[4],
					   double &tau, double &xau, double &yau, 
					   double &zau, double &kyau, 
					   double aatol, double rtol);

void raytracing(double spin, double epsilon, 
	            double orbit_omega, double r_orbit, 
	            double r_spot, double period, double hstart,
				double tt, double robs, double rstep2, 
				double pobs, double pstep, double iobs, 
				double dobs, double aatol, double rtol, 
				double E_obs[N_E], double fphi[N_T][N_E],
				double fphi_1[N_T][N_E],
				vector<string>  ofnames, int &hit_sum);

void spot(double orbit_omega, double r_orbit, 
	      double r_spot, double period,
	      double xem[4], double &D, int &check);

double integral(int N, double x[], double y[]);

int timer(double, double, double);

void intersection_new(double t_1, double r_1, double theta_1, 
					  double phi_1, double t_2, double r_2, 
					  double theta_2, double phi_2, double x_eq[]);

double horizon_JP(double spin, double epsilon, double r, 
	              double theta);

int spec(double epsilon, double spin, double iobs_deg, 
		 double r_spot, double r_multi_isco);

double r_omega(double spin, double epsilon, double radius);

int radial_ptn (double spin, double gg, double r,
	            double &Vrr, double &Vzz);

double r_orbit_of_omega(double spin, double epsilon, 
						double omega, double isco);

int shift_peak(double prd, double T_obs[], double int_spec[], 
	           double spectra[N_T][N_E]);

int outp_spec(double spectra[N_T][N_E], 
			  double outp_spectra[N_T_PRT][N_E_PRT]);

#endif

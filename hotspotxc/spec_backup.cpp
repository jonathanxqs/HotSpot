#ifndef _DEF_H
#include "def.h"
#endif

int spec(double epsilon, double spin, double iobs_deg, 
		 double r_spot, double r_multi_isco)
{
	int i, j, N_robs;
	int Ntot = N_T*N_E;
	int hit_sum;
	int hit_tmp;

	double robs_min, robs_max;
	double rstep_check, pstep_check;
	double spin2;
	double D, D2;
	double N_tot[N_T];
	double iobs, dobs;
	double robs, pobs;
	double robs_i, robs_f, rstep, rstep2, pstep, p_i, p_f;
	double isco;
	double tt;
	double ttmin;
	double timeduration = 0; 
	double nettimeduration = 0;
	double remaintime = 0;
	double Estep;
	double intspectra[N_T];
	double intintensity;
	double tempfphi;
	double E_obs[N_E + 1], E_shift[N_E + 1];
	double N_obs[N_T];
	double T_obs[N_T];
	double fphi[N_T][N_E];
	double fphi_1[N_T][N_E];
	double aatol, rtol;
	double orbit_omega, r_orbit, period;
	double hstart;

	// clock_t timestart, timefinish;
	
	/* ------------- The Region for Codes Test --------------- */
	//cout<<sizeof(double)<<endl;
	/*double a[1024][100];
	int ii, jj;

	for (ii = 0; ii < 1000; ii++)
		for (jj = 0; jj < 1000; jj++)  a[ii][jj] = 1.0;*/

	/*int a[3][4]={1,3,5,7,9,11,13,15,17,19,21,23};
    int *p;
    for(p=a[0];p<a[0]+12;p++)
    cout <<*p<<" ";
    cout <<endl;
	system("pause");
    return 0;*/

	/*for(int test = 0; test < 1000; test++);
	cout<<test<<endl;*/

	/*epsilon = 5;
    for (spin = 0; spin <= 2; spin = spin + 0.01) {
    	spin2 = spin*spin;
    	cout << spin << "\t" << epsilon << "\t";
		find_isco(spin, spin2, epsilon, isco);
		cout << "\t" << isco << endl;
	}
	*/
	/* ------------------------------------------------------- */
   
    //cout<<"\n Hello, World!\n"<<endl;

    /* -------------- parameter initialization --------------- */	
	ttmin = 0;
    /* ------------------------------------------------------- */

	/* --------------- set free parameters ------------------- */

	//cout<<"Please enter epsilon, the deformation parameter"<<endl;
	//cin>>epsilon;        /* deformation parameter */  
	//cout<<"Please enter spin, the spin parameter"<<endl;
	//cin>>spin;      /* spin parameter */
	spin2 = spin*spin;
 	
    /* ----------- set parameters of the HOT SPOT ------------ */
    orbit_omega = 0;  /* the angular velocity of hot spot, will be 
				assigned further in function find_isco */
	isco = find_isco(spin, epsilon);
    r_orbit = r_multi_isco*isco;

    orbit_omega = r_omega(spin, epsilon, r_orbit);

	period = 2*Pi/orbit_omega;
	//period = 150;

    // cout << setprecision(2);
	cout << " ----------- PARAMETERS USED IN THIS STEP ----------- " << endl;
	cout << "\t" << "epsilon" << "\t   " << "spin" << "\t   " << "iobs_deg" 
	     << "\t" << "r_spot" << endl;
	cout << "\t   " << epsilon << "\t    " << spin << "\t      " << iobs_deg 
		 << "\t " << r_spot << endl;
	cout << " ---------------------------------------------------- " << endl;
	cout << "\t" << "r_isco" << "\t\t" << "r_orbit" << "\t\t" << "period" << endl;
	cout << "\t" << isco << "\t\t" << r_orbit << "\t\t" << period << endl;
	cout << " ---------------------------------------------------- " << endl;

	/*cout<<"\n"<<"Well, the ISCO and the PERIOD of the spot are: \n"
	<<" ISCO = "<<setprecision(12)<<isco
	<<"\n PERIOD = "<<period<<"\n\n"<<endl;*/
    
	//cout<<"Please enter iobs_deg, the inclination angle (degree)"<<endl;
	//cin>>iobs_deg;     /* inclination angle */
	//cout<<"And, please enter the beginning time ttmin"<<endl;
	//cin>>ttmin;
	//cout<<"Okay, please enter the R_ORBIT of Hot Spot"<<endl;
	//cin>>r_spot;
	  
    // cout<<"Very good! Everything is okay now, and let's go!\n"<<endl;
	// cout<<"=======================================================\n"<<endl;
	// cout<<"======= PLEASE FEEL FREE TO FIND SOME PLEASURE ========\n"<<endl; 

	iobs = Pi/180*iobs_deg;     /* inclination angle of the observer in rad */
	D  = 10.0;     /* distance Earth-binary system in kpc */
	D2 = D*D; 
			
	/* ----- Set model for the spectral line ----- */

	/* ------------ open files for data recording ------------ */
	stringstream docname;
	docname<<"data"; // <<"_e"<<epsilon<<"_a"<<spin<<"_i"<<iobs_deg;
	//_mkdir(docname.str().c_str());

	stringstream ofname_spectra, ofname_intspec, ofname_intspec_1, ofname_logfile;
	/* ofname_spectra<<docname.str()<<"/"<<"spectra"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<".dat"; */
	ofname_intspec<<docname.str()<<"/"<<"tot_intspec"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<"_Rorbit"<<r_multi_isco<<".dat";
	ofname_intspec_1<<docname.str()<<"/"<<"pri_intspec"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<"_Rorbit"<<r_multi_isco<<".dat";
	ofname_logfile<<docname.str()<<"/"<<"logfile"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<"_Rorbit"<<r_multi_isco<<".txt";

	// ofstream of_spectra(ofname_spectra.str().c_str());
	ofstream of_intspec(ofname_intspec.str().c_str());
	ofstream of_intspec_1(ofname_intspec_1.str().c_str());
	ofstream of_logfile(ofname_logfile.str().c_str());

	// docname<<"/"<<"gifdata"<<"_e"<<epsilon<<"_a"<<spin<<"_i"<<iobs_deg;
	// _mkdir(docname.str().c_str());
    of_intspec << setprecision(10);
	of_intspec << epsilon << "\t" << epsilon << endl;
	of_intspec << spin << "\t" << spin << endl;
	of_intspec << iobs_deg << "\t" << iobs_deg << endl;
	of_intspec << r_spot << "\t" << r_spot << endl;
	of_intspec << r_orbit << "\t" << r_orbit << endl;
	of_intspec << "\n" << endl;

	of_intspec_1 << setprecision(10);
	of_intspec_1 << epsilon << "\t" << epsilon << endl;
	of_intspec_1 << spin << "\t" << spin << endl;
	of_intspec_1 << iobs_deg << "\t" << iobs_deg << endl;
	of_intspec_1 << r_spot << "\t" << r_spot << endl;
	of_intspec_1 << r_orbit << "\t" << r_orbit << endl;
	of_intspec_1 << "\n" << endl;

	vector<string>  ofname_gif;
	/* for(int i = 0; i < N_gif; ++i){
		stringstream ofname_tmp;
		ofname_tmp<<docname.str()<<"/"<<"gifdata_"<<i<<".dat";
		ofname_gif.push_back(ofname_tmp.str());
	    // cout<<ofname_tmp.str()<<endl;
	} */

	// system("pause");

	double r_N = 100;
	double phi_N = 100;

	/* ------------ set computational parameters ------------- */	
	dobs = 1000000;    /* distance of the observer */
	
	robs_i = 0.5;
	robs_f = 20;
	p_i = 0;
	p_f = 2*Pi;	
	rstep  = 1.005;
	rstep2 = (rstep*rstep - 1)/2/rstep;
	pstep  = 2*Pi/phi_N;  

	aatol = 1.0e-6;
	rtol = 1.0e-6;   
	hstart = 100;  

	rstep_check = 1.1;
	pstep_check = pstep*10;
    robs_min = robs_i;
    robs_max = robs_f;

 	hit_tmp = 0;
	for (robs = robs_f; robs > robs_i; robs = robs/rstep_check) {
		pobs = p_i; 
		tt = ttmin;
		hit_sum = 0;
		for (pobs = p_i; pobs < p_f - 0.5*pstep; pobs = pobs + pstep_check) {
			raytracing(spin, epsilon, orbit_omega, r_orbit, 
	            r_spot, period, hstart, tt, robs, rstep2, 
				pobs, pstep, iobs, dobs, aatol, rtol, E_obs, 
				fphi, fphi_1, ofname_gif, hit_sum);
		}
		// cout << robs << "\t" << hit_sum <<endl;
        
        if (hit_sum > 0 && hit_tmp == 0) {
        	robs_max = robs*rstep_check;
        }
        if (hit_sum == 0 && hit_tmp > 0) {
        	robs_min = robs;
        	break;
        }
        hit_tmp = hit_sum;
	}
	// cout << robs_max << "\t" << robs_min <<endl;
    
    robs_i = robs_min;
    robs_f = robs_max;
	rstep  = exp(log(robs_f/robs_i)/r_N);
	rstep2 = (rstep*rstep - 1)/2/rstep;
    N_robs = (int)ceil(log(robs_f/robs_i)/log(rstep)); 

	// system("pause");


    /* --------- keep every parameter in a log file ---------- */
    of_logfile<<"\nHello, World!\n"<<endl;
    of_logfile<<"Welcome to this simple codes to calculate spectra "
		<<"emitted by a hot spot orbiting a Black hole and observed "
		<<"by distant observer. ";
	of_logfile<<"As follows you can find parameters used in this "
	    <<"specific calculation: \n"<<endl;
    of_logfile<<"\tDeformation parameter: e = "<<epsilon<<endl;
    of_logfile<<"\tSpin parameter: a = "<<spin<<endl;
    of_logfile<<"\nThe ISCO, r_orbit and the PERIOD of the spot are: \n"
		<<"\t ISCO = "<<setprecision(12)<<isco
		<<"\n\n\t r_orbit = "<<r_orbit<<"\n"
		<<"\t PERIOD = "<<period<<"\n"<<endl;
    of_logfile<<"\tInclination angle: i = "<<iobs_deg<<" degrees"<<endl;
    of_logfile<<"\tBeginning time: ttmin = "<<ttmin<<endl;
    of_logfile<<"\tSpot R_ORBIT: R_spot = "<<r_spot<<" M"<<endl;
    of_logfile<<"\tObserver distance: d_obs = "<<dobs<<endl;
    of_logfile<<"\nThe scanning region in observer's plane: \n";
    of_logfile<<"\tr_initial = "<<robs_i<<"\n\tr_end = "<<robs_f
	    <<"\n\tr_step = "<<rstep;
    of_logfile<<"\n\tphi_initial = "<<p_i<<"\n\tphi_end = "<<p_f
	    <<"\n\tphi_step = "<<pstep<<endl;
    of_logfile<<"\nThe total r loops scanned in observer plane is:\n\t"
	    <<N_robs<<endl;
    of_logfile<<"\nThe initial h step used in Runge-Kutta-Nystrom "
	    <<"methods is \n\t h_start = "<<hstart<<endl;
	of_logfile<<"\nThe tolerance value used in Runge-Kutta-Nystrom "
	    <<"methods is \n\t a_tol = "<<aatol<<" \n\t r_tol = "
	    <<rtol<<endl;   
	time_t writtentime = time(NULL);
//    tm* timenow = localtime(&writtentime);
	tm* timenow = NULL;
	localtime_s(timenow,&writtentime);
    of_logfile<<"\n------------------------  "<<timenow->tm_year + 1900
		<<"-"<<timenow->tm_mon + 1<<"-"<<timenow->tm_mday<<" "
		<<timenow->tm_hour<<":"<<timenow->tm_min<<":"<<timenow->tm_sec
		<<endl;

	of_logfile.close();
	 
	Estep = double(2)/N_E;    /* minimum photon energy detected
							   by the observer; in keV */		
	for (i = 0; i < N_E + 1; i++) {
		E_obs[i] = i*Estep;
		if(i > 0) {
			E_shift[i-1] = E_obs[i];
		}
	}
	
	for (i = 0; i < N_T; i++) {
		N_obs[i] = 0;
	}

	for (i = 0; i < N_T; i++) {
		/*T_obs[i] = i*period/N_T/period;*/
		T_obs[i] = (double) i*T_MAX/N_T;
	}

	//for (i = 0; i < N_T; i++) {
	//	N_tot[i]  = 0;
	//}

	for (i = 0; i < Ntot; i++) {
		*(fphi[0]+i) = 0;
		*(fphi_1[0]+i) = 0;
	}

	/* --------- assign photon position in the grid ---------- */
	robs = robs_i;
	for (i = 0; i < N_robs; robs = robs*rstep, ++i) {
		
		//timestart = clock();    /* calculate the time used in each loop */
		/*cout<<"=======================================================\n"<<endl;
		cout<<"The spectra at robs = "<<robs<<" is under calculating!\n"<<endl;*/
		
		hit_sum = 0;
		for (pobs = p_i; pobs < p_f - 0.5*pstep; pobs = pobs + pstep) {
			//cout<<"\t"<<robs<<"\t"<<pobs<<"\n\n"<<endl; 
			tt = ttmin;
			raytracing(spin, epsilon, orbit_omega, r_orbit, 
	            r_spot, period, hstart, tt, robs, rstep2, 
				pobs, pstep, iobs, dobs, aatol, rtol, E_obs, 
				fphi, fphi_1, ofname_gif, hit_sum);
		}
		/*timefinish = clock();
        
		timeduration = (double)(timefinish - timestart)/CLOCKS_PER_SEC;
		nettimeduration += timeduration;
		remaintime = timeduration*(N_robs-i);

		cout<<"\n\n------------- This LOOP is robs = "<<robs<<"----------------"<<endl;
		cout<<"\n\n------------- "<<i + 1<<" LOOPs have been calculated------------"<<endl;
		cout<<"\n\n---------- There are still "<<N_robs - i<<" LOOPs waiting ----------"<<endl;
		
	    timer(timeduration, nettimeduration, remaintime);*/
	}
	
	/* --------- data rearrangement and outputting: TOTAL----------- */
	//for (i = 0; i < N_T; i++) {
	//	for (j = 0; j < N_E; j++) {
	//		tempfphi = *(fphi[i] + j);
	//		N_tot[i] = N_tot[i] + tempfphi;	
	//		
	//		/* of_spectra<<setiosflags(ios::fixed)<<setprecision(3)
	//			<<T_obs[i]<<"\t"<<setprecision(3)<<E_shift[j]
	//		    <<"\t"<<setprecision(10)<<tempfphi<<"\n"; */
	//		}	
	//	}
	
	for (i = 0; i < N_T; i++) {
		intspectra[i] = 0;
		intspectra[i] = integral(sizeof(fphi[0])/sizeof(fphi[0][0]), E_shift, fphi[i]);
	}

	int i_zero, n_period;
	double tstep, n_tt;
	tstep = (double) T_MAX/N_T;
	n_period = int (T_MAX/period);
	n_tt = n_period*period/tstep;
	i_zero = floor(n_tt);

	intintensity = integral(N_T, T_obs, intspectra);
	for (i = 0; i < N_T; i++) {
    	intspectra[i] = n_period*period*intspectra[i]/intintensity;
	}

	vector <double> extend_spec;
	int extend_n_t = (N_T/(i_zero + 1) + 1)*(i_zero + 1);
	for (i = 0; i < extend_n_t; ++i) 
		extend_spec.push_back(intspectra[i % (i_zero + 1)]);

	vector <double> extend_spec_copy(extend_spec);

	double spec_peak;
	int peak_post, peak_shft;

	spec_peak = 0;
	peak_post = 0;
	for (i = 0; i < extend_n_t; i++) {
    	if (extend_spec[i] >= spec_peak) {
    		spec_peak = extend_spec[i];
    		peak_post = i;
    	}
	}

	// int ii_min, ii_max, ii_range;

	peak_shft = N_T/2 - peak_post;
	for (i = 0; i < extend_n_t; i++) {
		j = ( - peak_shft + i + extend_n_t) % extend_n_t;
    	extend_spec[i] = extend_spec_copy[j];
	}

    for (i = 0; i < N_T; i++) {
		of_intspec<<setiosflags(ios::fixed)<<setprecision(10)
			<<T_obs[i]<<"\t"<<setprecision(10)<<extend_spec[i]<<"\n";
	}

	// of_spectra.close();
    of_intspec.close();


	/* --------- data rearrangement and outputting: PRIMARY----------- */
	//for (i = 0; i < N_T; i++) {
	//	for (j = 0; j < N_E; j++) {
	//		tempfphi = *(fphi_1[i] + j);
	//		N_tot[i] = N_tot[i] + tempfphi;	
	//		
	//		/* of_spectra<<setiosflags(ios::fixed)<<setprecision(3)
	//			<<T_obs[i]<<"\t"<<setprecision(3)<<E_shift[j]
	//		    <<"\t"<<setprecision(10)<<tempfphi<<"\n"; */
	//		}	
	//	}
	
	for (i = 0; i < N_T; i++) {
		intspectra[i] = 0;
		intspectra[i] = integral(sizeof(fphi_1[0])/sizeof(fphi_1[0][0]), E_shift, fphi_1[i]);
	}

	tstep = (double) T_MAX/N_T;
	n_period = int (T_MAX/period);
	n_tt = n_period*period/tstep;
	i_zero = floor(n_tt);

	intintensity = integral(N_T, T_obs, intspectra);
	for (i = 0; i < N_T; i++) {
    	intspectra[i] = n_period*period*intspectra[i]/intintensity;
	}

	vector <double> extend_spec_1;
	extend_n_t = (N_T/(i_zero + 1) + 1)*(i_zero + 1);
	for (i = 0; i < extend_n_t; ++i) 
		extend_spec_1.push_back(intspectra[i % (i_zero + 1)]);

	vector <double> extend_spec_copy_1(extend_spec_1);

	spec_peak = 0;
	peak_post = 0;
	for (i = 0; i < extend_n_t; i++) {
    	if (extend_spec_1[i] >= spec_peak) {
    		spec_peak = extend_spec_1[i];
    		peak_post = i;
    	}
	}

	// int ii_min, ii_max, ii_range;

	peak_shft = N_T/2 - peak_post;
	for (i = 0; i < extend_n_t; i++) {
		j = ( - peak_shft + i + extend_n_t) % extend_n_t;
    	extend_spec_1[i] = extend_spec_copy_1[j];
	}

    for (i = 0; i < N_T; i++) {
		of_intspec_1<<setiosflags(ios::fixed)<<setprecision(10)
			<<T_obs[i]<<"\t"<<setprecision(10)<<extend_spec_1[i]<<"\n";
	}

	// of_spectra.close();
    of_intspec_1.close();
		
	/*cout<<"=======================================================\n"<<endl;
	cout<<"===== CONGRADULATIONS! The calculation is FINISHED! ====\n\n"<<endl;*/



	return 0;
}

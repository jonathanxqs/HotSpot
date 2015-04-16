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
	int size_tmp;

	double robs_min, robs_max;
	double rstep_check, pstep_check;
	double spin2;
	double D, D2;
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
	double intspec_tot[N_T], intspec_pri[N_T];
	double E_obs[N_E + 1], E_shift[N_E + 1];
	double N_obs[N_T];
	double T_obs[N_T];
	double fphi_tot[N_T][N_E];
	double fphi_pri[N_T][N_E];
	double outp_fphi_tot[N_T_PRT][N_E_PRT];
	double outp_fphi_pri[N_T_PRT][N_E_PRT];
	double aatol, rtol;
	double orbit_omega, r_orbit, period;
	double hstart;

	// clock_t timestart, timefinish;

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

	if (!(period < T_MAX - 2 || period > 0 ) )  {
		cout << epsilon << "\t " << spin << "\t " << isco 
			 << "\t " << orbit_omega << "\t " << period <<endl;
		return 0;
	}

	//period = 150;

    cout << std::fixed << setprecision(3) << endl;
	cout << " ----------- PARAMETERS USED IN THIS STEP ----------- " << endl;
	cout << "\t" << "epsilon" << "\t  " << "spin" << "\t  " << "iobs_deg" 
	     << "\t" << "r_spot" << endl;
	cout << "\t" << epsilon << "\t " << spin << " \t " << iobs_deg 
		 << " \t " << r_spot << endl;
	cout << std::fixed << setprecision(4) << endl;	
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

	stringstream ofname_spectra;
	stringstream ofname_intspec;
	stringstream ofname_logfile;
	/* ofname_spectra<<docname.str()<<"/"<<"spectra"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<".dat"; */
	ofname_spectra<<docname.str()<<"/"<<"spectra"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<"_Rorbit"<<r_multi_isco<<".dat";
	ofname_intspec<<docname.str()<<"/"<<"intspec"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<"_Rorbit"<<r_multi_isco<<".dat";
	ofname_logfile<<docname.str()<<"/"<<"logfile"<<"_e"<<epsilon<<"_a"
		<<spin<<"_i"<<iobs_deg<<"_Rspot"<<r_spot<<"_Rorbit"<<r_multi_isco<<".txt";

	// ofstream of_spectra(ofname_spectra.str().c_str());
	ofstream of_spectra(ofname_spectra.str().c_str());
	ofstream of_intspec(ofname_intspec.str().c_str());
	ofstream of_logfile(ofname_logfile.str().c_str());

	// docname<<"/"<<"gifdata"<<"_e"<<epsilon<<"_a"<<spin<<"_i"<<iobs_deg;
	// _mkdir(docname.str().c_str());
	of_spectra << setprecision(6)
		<< epsilon << "\t" << epsilon << "\t" << epsilon << "\t" << epsilon << "\n"
		<< spin << "\t" << spin << "\t" << spin << "\t" << spin << "\n"
		<< iobs_deg << "\t" << iobs_deg << "\t" << iobs_deg << "\t" << iobs_deg << "\n"
	    << r_spot << "\t" << r_spot << "\t" << r_spot << "\t" << r_spot << "\n"
		<< r_orbit << "\t" << r_orbit << "\t" << r_orbit << "\t" << r_orbit << "\n"
		<< isco << "\t" << isco << "\t" << isco << "\t" << isco << "\n"
		<< period << "\t" << period << "\t" << period << "\t" << period << "\n"
		<< endl;

	of_intspec << setprecision(6)
		<< epsilon << "\t" << epsilon << "\t" << epsilon << "\n"
		<< spin << "\t" << spin << "\t" << spin << "\n"
		<< iobs_deg << "\t" << iobs_deg << "\t" << iobs_deg << "\n"
	    << r_spot << "\t" << r_spot << "\t" << r_spot << "\n"
		<< r_orbit << "\t" << r_orbit << "\t" << r_orbit << "\n"
		<< isco << "\t" << isco << "\t" << isco << "\n"
		<< period << "\t" << period << "\t" << period << "\n"
		<< endl;

	vector<string>  ofname_gif;
	/* for(int i = 0; i < N_gif; ++i){
		stringstream ofname_tmp;
		ofname_tmp<<docname.str()<<"/"<<"gifdata_"<<i<<".dat";
		ofname_gif.push_back(ofname_tmp.str());
	    // cout<<ofname_tmp.str()<<endl;
	} */

	/* ------------ set computational parameters ------------- */	
	dobs = 1000000;    /* distance of the observer */
	
	robs_i = 0.5;
	robs_f = 20;
	p_i = 0;
	p_f = 2*Pi;	
	rstep  = 1.005;
	rstep2 = (rstep*rstep - 1)/2/rstep;
	pstep  = 2*Pi/PHI_N;  

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
				fphi_tot, fphi_pri, ofname_gif, hit_sum);
		}
        
        if (hit_sum > 0 && hit_tmp == 0) {
        	robs_max = robs*rstep_check;
        }
        if (hit_sum == 0 && hit_tmp > 0) {
        	robs_min = robs;
        	break;
        }
        hit_tmp = hit_sum;
	}
    
    robs_i = robs_min;
    robs_f = robs_max;
	rstep  = exp(log(robs_f/robs_i)/R_N);
	rstep2 = (rstep*rstep - 1)/2/rstep;
    N_robs = (int)ceil(log(robs_f/robs_i)/log(rstep)); 

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
    of_logfile<<"\tSpot radius: R_spot = "<<r_spot<<" M"<<endl;
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
    //tm* timenow = localtime(&writtentime);
	tm* timenow = NULL;
	localtime_s(timenow,&writtentime);
    of_logfile<<"\n------------------------  "<<timenow->tm_year + 1900
		<<"-"<<timenow->tm_mon + 1<<"-"<<timenow->tm_mday<<" "
		<<timenow->tm_hour<<":"<<timenow->tm_min<<":"<<timenow->tm_sec
		<<endl;
	of_logfile.close();
	 

	Estep = double(G_MAX)*1.0/N_E;    /* minimum photon energy detected
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
		T_obs[i] = (double) i*T_MAX/N_T;
	}

	for (i = 0; i < Ntot; i++) {
		*(fphi_tot[0]+i) = 0;
		*(fphi_pri[0]+i) = 0;
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
				fphi_tot, fphi_pri, ofname_gif, hit_sum);
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
	
	/* --------- data rearrangement and outputting ----------------- */

	for (i = 0; i < N_T; i++) {
		intspec_tot[i] = 0;
		size_tmp = sizeof(fphi_tot[0])/sizeof(fphi_tot[0][0]);
		intspec_tot[i] = integral(size_tmp, E_shift, fphi_tot[i]);
		intspec_pri[i] = 0;
		size_tmp = sizeof(fphi_pri[0])/sizeof(fphi_pri[0][0]);
		intspec_pri[i] = integral(size_tmp, E_shift, fphi_pri[i]);
	}

	shift_peak(period, T_obs, intspec_tot, fphi_tot);
	shift_peak(period, T_obs, intspec_pri, fphi_pri);

    for (i = 0; i < N_T; i++) {
		of_intspec << setiosflags(ios::fixed) 
			<< setprecision(4) << T_obs[i] << " \t " 
			<< setprecision(10) << intspec_pri[i] << " \t "
			<< setprecision(10) << intspec_tot[i] << "\n";
	}

	double* ptr_pri = outp_fphi_pri[0];
	double* ptr_tot = outp_fphi_tot[0];
	for (i = 0; i < N_T_PRT; i++) {
		for (j = 0; j < N_E_PRT; j++) {
			*ptr_pri++ = 0;
			*ptr_tot++ = 0;
		}
	}

	outp_spec(fphi_tot, outp_fphi_tot);
	outp_spec(fphi_pri, outp_fphi_pri);

	int E_stp_prt = N_E/N_E_PRT;
	int T_stp_prt = N_T/N_T_PRT;
	ptr_pri = outp_fphi_pri[0];
	ptr_tot = outp_fphi_tot[0];
	for (i = 0; i < N_T_PRT; i++) {
		for (j = 0; j < N_E_PRT; j++) {
			of_spectra << setiosflags(ios::fixed) 
				<< setprecision(4)
				<< T_obs[i*T_stp_prt + T_stp_prt/2] <<" \t " 
				<< setprecision(4) 
				<< E_shift[j*E_stp_prt + E_stp_prt/2] <<" \t " 
				<< setprecision(10) << *ptr_pri++ << " \t "
				<< setprecision(10) << *ptr_tot++ << "\n";
		}
	}

	of_spectra.close();
    of_intspec.close();
	
	/*cout<<"=======================================================\n"<<endl;
	cout<<"===== CONGRADULATIONS! The calculation is FINISHED! ====\n\n"<<endl;*/

	return 0;
}

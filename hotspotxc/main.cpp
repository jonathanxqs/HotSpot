#ifndef _DEF_H
#include "def.h"
#endif

int main(int argc, char **argv)
{
	double epsilon;
	double spin;
	double iobs_deg;
	double r_spot;
	double r_multi_isco;
	double timeduration = 0; 
	double nettimeduration = 0;
	double hurs, mins, secs;
	double a_min, a_max, a_step;
	double e_min, e_max, e_step;
	clock_t timestart, timefinish;

	epsilon = 0;      
	spin = 0;      
	iobs_deg = 60; 
	r_spot = 0.3;
	r_multi_isco = 1.0;

    /* default scan area in the plane spin-epsilon */
	a_min = 0;
	a_max = 1;
	a_step = 0.9;
	e_min = 0;
	e_max = 1;
	e_step = 2;

	/*a_min = atof(argv[1]);
    a_max = atof(argv[2]);
    a_step = atof(argv[3]);
    e_min = atof(argv[4]);
    e_max = atof(argv[5]);
    e_step = atof(argv[6]);*/
    
	double isco_0, isco;
	double omega_0, omega;

	isco_0 = find_isco(0, 0);
	omega_0 = r_omega(0, 0, isco_0);

	//std::cin >> epsilon >> spin >> iobs_deg >> r_spot;
	
	//for (spin = a_max; spin >= a_min; spin -= a_step)
		//for (epsilon = e_min; epsilon < e_max; epsilon += e_step)
			// for (iobs_deg = 20; iobs_deg < 90; iobs_deg += 25)
				// for (r_spot = 0.05; r_spot < 0.16; r_spot += 0.05)
					//for (r_multi_isco = 1.5; r_multi_isco < 1.6; r_multi_isco += 0.2)
		{
		    cout << "\n======================================================" << endl;

			timestart = clock(); 

			isco = find_isco(spin, epsilon);
			//// omega = r_omega(spin, epsilon, isco);
			r_multi_isco = r_orbit_of_omega(spin, epsilon, omega_0, isco);
			
			spec(epsilon, spin, iobs_deg, r_spot, r_multi_isco);

			timefinish = clock();
			cout << " ------ CONGRADULATIONS! This STEP is FINISHED! ----- "<<endl;
		    
		    timeduration = (double)(timefinish - timestart)/CLOCKS_PER_SEC;
			nettimeduration += timeduration;	    
			secs = floor(fmod(timeduration, 60)); 
		    mins = floor(timeduration/60);
			hurs = floor(mins/60);
		    mins = floor(fmod(mins, 60));
			
			cout << std::fixed << setprecision(0) << endl;
		    cout << " \n \t Elapsed Time:";
		    cout << "\t" << hurs << " hurs " << mins << " mins " 
		         << secs << " secs \t\n" << endl;
			cout << "======================================================\n" << endl;
			cout << "                 ********************                 \n" << endl;
		}

	cout << "\n======================================================\n" << endl;
	secs = floor(fmod(nettimeduration, 60));    
    mins = floor(nettimeduration/60);
	hurs = floor(mins/60);
    mins = floor(fmod(mins, 60));

	cout << std::fixed << setprecision(0) << endl;
    cout << " \t  Net Time:" ;
    cout << "\t" << hurs << " hurs " << mins << " mins " 
         << secs << " secs \t\n" << endl;
	cout << "======================================================\n" << endl;

	system("pause");
	return 0;

	/*ofstream of_logfile("periods.txt");
    double orbit_omega, r_orbit, isco, period;
	orbit_omega = 0;
	for (epsilon = e_min; epsilon < e_max; epsilon += e_step)
		for (spin = a_min; spin < a_max; spin += a_step)  {
			isco = find_isco(spin, epsilon, orbit_omega);
    		r_orbit = isco;
    		orbit_omega = r_omega(spin, epsilon, r_orbit);
    		period = 2*Pi/orbit_omega;
    		of_logfile << spin << "\t" << epsilon << "\t" 
    			       << orbit_omega << "\t" << period << endl;
    	    cout << spin << "\t" << epsilon << "\t" 
    			 << orbit_omega << "\t" << period << endl;
		}
	of_logfile.close();*/
}

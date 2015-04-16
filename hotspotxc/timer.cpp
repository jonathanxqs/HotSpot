#include <math.h>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;

int timer(double t_loop, double t_net, double t_remain)
{
	double hurs, mins, secs;
	cout << setprecision(0) << endl;
    
	secs = floor(fmod(t_loop, 60));    
    mins = floor(t_loop/60);
	hurs = floor(mins/60);
    mins = floor(fmod(mins, 60));
    cout << "Elapsed time of this single loop is: " << endl;
    cout << "\t" << hurs << " hurs " << mins << " mins " 
         << secs << " secs \n" << endl;
        
    secs = floor(fmod(t_net, 60));    
    mins = floor(t_net/60);
	hurs = floor(mins/60);
    mins = floor(fmod(mins, 60));
    cout << "Net elapsed time till now is: " << endl;
    cout << "\t" << hurs << " hurs " << mins << " mins " 
         << secs << " secs \n" << endl;       

    secs = floor(fmod(t_remain, 60));    
    mins = floor(t_remain/60);
	hurs = floor(mins/60);
    mins = floor(fmod(mins, 60));
    cout << "The estimated time remaining is: " << endl;
    cout << "\t" << hurs << " hurs " << mins << " mins " 
         << secs << " secs \n" << endl;   

	return 0;
}

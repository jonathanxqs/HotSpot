/* 
* This file will define a function named metric_JP, which can
* be used to calculate the no-vanishing conponents of JP metric.
* Johannsen-Psaltis (JP) metric can be seen as a metric 
* describing non-Kerr BHs in a putative alternative theory of 
* gravity [T. Johannsen and D. Psaltis, Phys. Rev. D 83, 124015
* (2011)]. This metric has an infinite number of deformation 
* parameters epsilon_i, and the Kerr solution is recovered 
* when all the deformation parameters are set to zero.In order
* to reproduce the correct Newtonian limit, we have to impose 
* epsilon_0 = epsilon_1 = 0, while epsilon_2 is strongly 
* constrained by Solar System experiments. In this function, 
* I will only examine the simplest cases where epsilon_3 != 0, 
* while all the other deformation parameters are set to zero.
*
* 2013-11-26 by Zilong Li
* zilongli@fudan.edu.cn
*/

#include <math.h>

void metric_JP(double spin, double epsilon, double r, 
	           double theta, double mn[4][4])
{
	double spin2 = spin*spin;
	double r2;
	double c2, s2;
	double Delta;
	double Sigma;
	double h;
	double fctr0, fctr1, fctr2;

	r2 = r*r;
	c2 = cos(theta)*cos(theta);
	s2 = 1 - c2;
	Delta = r2 - 2*r + spin2;
	Sigma  = r2 + spin2*c2;
	h = epsilon*r/Sigma/Sigma;
	fctr0 = 2*spin*r*s2/Sigma;
	fctr1 = r2 + spin2 + fctr0*spin;
	fctr2 = spin2*(Sigma + 2*r)*s2*s2/Sigma;
	
	mn[0][0] = - (1 + h)*(1 - 2*r/Sigma);
	mn[0][3] = - fctr0*(1 + h);
	mn[1][1] = (1 + h)*Sigma/(Delta + spin2*h*s2);
	mn[2][2] = Sigma;
	mn[3][0] = mn[0][3];
	mn[3][3] = s2*fctr1 + fctr2*h;

	return ;
}

double horizon_JP(double spin, double epsilon, double r,
                  double theta)
{
	double spin2 = spin*spin;
	double r2;
	double c2, s2;
	double Delta, Sigma;
	double h;
	double horizon;

	r2 = r*r;
	c2 = cos(theta)*cos(theta);
	s2 = 1 - c2;
	Delta = r2 - 2*r + spin2;
	Sigma  = r2 + spin2*c2;
	h = epsilon*r/Sigma/Sigma;
	horizon = Delta + spin2*s2*h;

	return horizon;
}

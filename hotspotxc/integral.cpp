#ifndef _DEF_H
#include "def.h"
#endif

double integral(int N, double x[], double y[])
{
	int i;
	double II;

	II = 0;
	for (i = 0; i < N - 1; i++) {
		II += double((x[i+1] - x[i])*(y[i+1] + y[i])/2);
	}

	return II;
}

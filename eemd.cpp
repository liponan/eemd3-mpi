// EEMD  coded by Po-Nan Li 2012
// requires: emd_core.cpp 

#include <cstdlib>
#include <cmath>
#include <ctime>
#include "emd_core.cpp"
#include <iostream>

double u = 0;
double v = 0;

void randn(double *W, int sz) {
	for (int i = 0; i < sz; i++) {
		u = (rand() *1.0) / RAND_MAX;
		v = (rand() *1.0) / RAND_MAX; 
		W[i] = sqrt( -2 * log(u) ) * cos( 2 * 3.1415926 * v );
	} // end of for-i
}

double Std(double *Y, int sz) {
	double mean = 0;
	for (int i = 0; i < sz; i++)
		mean = mean + Y[i];
	mean = mean / sz;
	double sigma = 0;
	for (int i = 0; i < sz; i++)
		sigma = sigma + pow( (Y[i] - mean) , 2);
	if (sigma > 0)
		sigma = sqrt( sigma / sz );
	else
		sigma = 0;
	return sigma;
} // end of Std()

void eemd(double *modes,
		double *Y, int sz, int goal, int ens, double nstd) {

	srand((int)time(NULL));
	int m, i, c, k, t;	
	int goal1 = goal + 1;

	/* Core function */
	double *m1 = new double[goal1*sz];
	double *m2 = new double[goal1*sz];
	double *tmp = new double[goal1*sz];
	double *wn = new double[sz];
	double *Y1 = new double[sz];
	double *Y2 = new double[sz];
	double sigma = Std(Y, sz);
    
    for (t = 0; t < sz*goal1; t++)
        tmp[t] = 0;
	for (k = 0; k < ens; k++) {
		randn(wn, sz);
		for (i = 0; i < sz; i++) {
			if (sigma > 0.000000001)
				Y1[i] = Y[i]+ sigma * wn[i] * nstd;
			else
				Y1[i] = Y[i] + 0.000000001 * wn[i] * nstd;
			if (nstd > 0) {
				if (sigma > 0.000000001)
					Y2[i] = Y[i] - sigma * wn[i] * nstd;
				else
					Y2[i] = Y[i] - 0.000000001 * wn[i] * nstd;
			} // end of if 
		} // end of for-i

		emd_core(m1, Y1, sz, goal);
		if (nstd > 0)
			emd_core(m2, Y2, sz, goal);
		if (nstd > 0)
			for (t = 0; t < sz*goal1; t++)
				tmp[t] = tmp[t] + m1[t] + m2[t];
		else
			for (t = 0; t < sz*goal1; t++)
				tmp[t] = tmp[t] + m1[t];
	} // end of for-k

	if (nstd > 0)
		for (t = 0; t < sz*goal1; t++)
			modes[t] = tmp[t] / (ens * 2);
	else 
		for (t = 0; t < sz*goal1; t++)
			modes[t] = tmp[t] / ens;

	// free dynamic arrays	
	delete[] m1;
	delete[] m2;
	delete[] tmp;
	delete[] wn;
	delete[] Y1;
	delete[] Y2;

} // end of eemd
	

/* 1-D EMD   coded by Po-Nan Li 2012 */ 

#include <iostream>
#include <fstream>
#include "find_extrema.cpp"
//#include "spline_alglib.cpp"
#include "spline_gsl.cpp"

	using namespace std;

void mean(double *mean, double *y1, double *y2, int sz) {
	for (int i = 0; i < sz; i++)
		mean[i] = (y1[i] + y2[i]) / 2;
}

void emd_core(double *modes, 
		double *Y, int sz, int goal) {
	int m, i, c;
	int MAX = sz;
	int goal1 = goal + 1;

	/* Core function */

	// variables for find_extrema
	double *vmax = new double[MAX];
	double *vmin = new double[MAX];
	int *pmax = new int[MAX];
	int *pmin = new int[MAX];
	int lmax = 0, lmin = 0;
	// variables for spline
	double *upper = new double[MAX];
	double *lower = new double[MAX];
	double *emean = new double[MAX];
	int itr = (int)ceil(pow(sz, 0.5));
	double *h = new double[MAX];
	double *r = new double[MAX];


	for (i = 0; i < sz; i++)
		emean[i] = 0;

	// r = Y
	for(i = 0; i < sz; i++)
		r[i] = Y[i];
	// for each mode
	for (m = 0; m < goal; m++) {
			// h = r
			for (i = 0; i < sz; i++) 
				h[i] = r[i];
				
			// solving mode
			for (c = 0; c < itr; c++) {
				// find extremas
				find_extrema(pmax, vmax, pmin, vmin, &lmax, &lmin, h, sz);
				// findupper envelope
				spline(upper, pmax, vmax, lmax, sz);
				// find lower envelope
				spline(lower, pmin, vmin, lmin, sz);
				// find mean curvs
				mean(emean, upper, lower, sz);
				for (i = 0; i < sz; i++) {
					if (emean[i] > 1000 || emean[i] < -1000)
						cout<< "Warning! Abnormal EMEAN value " << emean[i] 
							<< " at position " << i << " after " << c << " iterations..." << endl;
				} // end of for-i
				// h = h - mean
				for (i = 0; i < sz; i++)
					h[i] = h[i] - emean[i];
				for (i = 0; i < sz; i++) {
					if (h[i] > 1000 || h[i] < -1000)
						cout<< "Warning! Abnormal H value " << emean[i] 
							<< " at position " << i << " after " << c << " iterations..." << endl;
				} // end of for-i
			} // end of for-c
			// r = r - h
			for (i = 0; i < sz; i++) {
				r[i] = r[i] - h[i];
				modes[i + m*sz] = h[i];
			} // end of for-i

	} // end of for-m

	for (i = 0; i < sz; i++) // save trend
		modes[i + goal*sz] = r[i];
    
	delete [] vmax;
	delete [] pmax;
	delete [] vmin;
	delete [] pmin;
	delete [] upper;
	delete [] lower;
	delete [] emean;
	delete [] h, r;
	
	} // end of emd_core
	

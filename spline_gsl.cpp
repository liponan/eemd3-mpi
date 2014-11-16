// requires GNU Scientific Lib
// v0.3 @ 2014/11/12

#include <gsl/gsl_spline.h>

void spline(double *YY, 
		int *X, double *Y, int m1, int m2) {

	// m1: length of discrete extremas
	// m2: length of original data points


	
  	//const gsl_interp_type *t = gsl_interp_cspline; 
  	

	/* Core function */
  	double *xd = new double[m1];
  	for (int j = 0; j < m1; j++)
  		xd[j] = (int)X[j];

  	gsl_spline_init (spline, xd, Y, m1);

	double m;
	
	if (m1 > 3 ) { // use spline
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, m1);
		gsl_spline_init (spline, xd, Y, m1);
		for (int j = 0; j < m2; j++) {
			YY[j] = gsl_spline_eval (spline, j, acc);
		} // end of for-j
	} // end of if

	else {
		if (m1 > 2) { // use polynomial
			gsl_interp_accel *acc = gsl_interp_accel_alloc ();
			gsl_spline *spline = gsl_spline_alloc (gsl_interp_polynomial, m1);
			gsl_spline_init (spline, xd, Y, m1);
			for (int j = 0; j < m2; j++) {
			YY[j] = gsl_spline_eval (spline, j, acc);
		} // end of for-j

		} // end of if
		else {
			m = (Y[1] - Y[0]) / (m2 - 1);
			for (int j = 0; j < m2; j++) {
				YY[j] = Y[0] + m * j;
		} // end of else	
		} // end of for-j	
	} //end of else

	delete[] xd;
	if (m1 > 2) {
		gsl_interp_accel_free (acc);
		gsl_spline_free (spline);	
	}
	

}

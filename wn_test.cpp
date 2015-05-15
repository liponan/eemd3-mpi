// eemd3d_mpi_v2.cpp (open-mpi)
// v1:   written by Neil Po-Nan Li 2014/04/01 @ IPAS
// v2:   Now support up to N^2 cores parallel computing, 
//       where N is the length in one dimension
//       2014/04/02
// v3:   Interpolation engine changes to GNU Scientific Library
//       2014/04/30
// v3.1: bugs fixed 2014/05/05
// v3.2: now export to binary file
// v3.6  read binary file (2014/11/11)
// v3.7  fixed bad_alloc problem


#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
/* custom subfunctions */


	using namespace std;


int main(int argc, char *argv[])
{
	int size = 10000;
	int count = 0;
	if (argc > 1)
		size = atoi( argv[1] );

	double *wn = new double[size];

	srand((int)time(NULL));

	
	for (int i = 0; i++; i < size) {
		wn[i] = (rand() *1.0) / RAND_MAX;
		if (wn[i] == 0)
			count++;
	}

	cout << "0 occurence prob: " << 1.0 * count / size << endl;

	return 0;
} // end of main()


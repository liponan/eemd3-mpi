// now up to 4-D volume is supported

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>

	using namespace std;

void printArray(string filename, double *Y, int dim, int* lg) {
	int U = 1;
	int V = 1;
	int W = 1;
	int M = 1;
	ofstream fout;
	fout.open(filename.c_str());
	switch (dim) {
	case 1:
		U = lg[0];
		break;
	case 2:
		U = lg[0];
		V = lg[1];
		break;
	case 3:
		U = lg[0];
		V = lg[1];
		W = lg[2];
		break;
	case 4:
		U = lg[0];
		V = lg[1];
		W = lg[2];
		M = lg[3];
	default:
		break;
	}

	int count = 0;
	fout << "A = zeros(1, " << U*V*W*M << ");" << endl;
	fout << "A(1, :) = [";

	/* write data */
	for (int m = 0; m < M; m++)	{
		for (int k = 0; k < W; k++) {
			for(int j = 0; j < V; j++) {
				for (int i = 0; i < U; i++) {
					fout << Y[i + j*U + k*U*V + m*U*V*W] << " ";
					if (++count % 10 == 0)
						fout << "... " << endl;
				} // for-j
			} // for-i
		} // for-k
	} // for-m
	fout << "]; " << endl;
	fout << "A = reshape(A, [" << U << " " << V << " " << W << " " << M << "]);" << endl;
	fout.close();
}
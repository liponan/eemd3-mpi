// now up to 4-D volume is supported

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>

	using namespace std;

void print2bin(string filename, double *Y, int dim, int* lg) {
	FILE *file;
	char filename_char[20];
	strcpy(filename_char, filename.c_str()); 
	file = fopen(filename_char , "wb");


	// first byte: dimension number
	fwrite(&dim, sizeof(int), 1, file);

	// 2nd~5th bytes: size in each dimension
	int sz = 1;
	for (int i = 0; i < dim; i++) {
		sz *= lg[i];
		fwrite(&lg[i], sizeof(int), 1, file);
	}

	// remaing bytes: data (double)
	fwrite(Y, sizeof(double), sz, file);	

	// close file
	fclose(file);
}
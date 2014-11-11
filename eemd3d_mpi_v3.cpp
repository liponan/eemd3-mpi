// eemd3d_mpi_v2.cpp (open-mpi)
// v1:   written by Neil Po-Nan Li 2014/04/01 @ IPAS
// v2:   Now support up to N^2 cores parallel computing, 
//       where N is the length in one dimension
//       2014/04/02
// v3:   Interpolation engine changes to GNU Scientific Library
//       2014/04/30
// v3.1: bugs fixed 2014/05/05
// v3.2: now export to binary file

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <mpi.h>
/* custom subfunctions */
#include "eemd.cpp"
#include "binaryIO.cpp"

	using namespace std;

	int toDo(int, int, int);

int main(int argc, char *argv[])
{
	// initialize MPI
	MPI_Init(&argc, &argv);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// some variables
	int t1 = MPI_Wtime();
	int t2, t3, t4, t5, t6, t7;
	int dt, eta_time;
	int dim = 0;
	int d = 0;
	int lg[4] = {0};
	bool flag = true;
	int U = 1, V = 1, W = 1, SZ = 1;
	int i = 0, j = 0, k =0;
	int t = 0, m1 = 0, m2 = 0, m3 = 0, r = 0;

	ifstream fin;
	string tmpStr;

	if (world_rank == 0) {
		
		// read the input file		
		if (argc < 2) {
			cout << "No input file!!" << endl;
			return 0;
		} else {
			cout << "Loading " << argv[1] << endl;
			bin_flag1 = readBinaryHeader(&dim, lg, argv[1]);
			if (bin_flag1 == 0)
				flag = false;
			else {
				cout << "Can't open " << argv[1] << ". Please try again... " << endl;
				return 1;
			}
		}
			
		// read the header
		U = lg[1];
		V = lg[0];
		W = lg[2];
		SZ = U * V * W;
	} // end of if (world_rank == 0)

	MPI_Bcast(&SZ,  1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the data length SZ
	MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the dimension dim
	MPI_Bcast(&U,   1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the height U
	MPI_Bcast(&V,   1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the width V
	MPI_Bcast(&W,   1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the depth W


	/* print the dimensions */
	if (world_rank == 0) {
		cout << "# of dimensions: " << dim << endl;
		cout << "Size: [ ";
		for (d = 0; d < dim; d++)
			cout << lg[d] << " ";
		cout << "]" << endl;
	} // end of if (world_rank == 0)

	

	// Now we know the exact data size, so let's declare the memory block for the data
	double *img = NULL;
	double tmpNum;

	// load the image into memory
	if (world_rank == 0) {
		img = new double[SZ];
		bin_flag2 = readBinaryImage(img, argv[1]);
	} // end of if (world_rank == 0)

	//MPI_Bcast(img, SZ, MPI_DOUBLE, 0, MPI_COMM_WORLD); // broadcast the IMG data to all nodes

	/* find the parameters or use the default values */
	
	// init
	int goal = 3;
	int ens = 1;
	double nstd = 0;

	int timecode = (int)time(NULL) % 10000;

	// You all know the argument values, so I don't need to broadcast them
	switch (argc) {
		case 6: 
			timecode = atoi(argv[5]);
		case 5:
			nstd = atof(argv[4]);
		case 4:
			ens = atoi(argv[3]);	
		case 3:
			goal = atoi(argv[2]);
		default:
			break;
	}
	int goalt = goal + 1;

	char timecode_str[5];
	sprintf(timecode_str, "%d", timecode);

	double *modes1; 
	double *rootBuff1;
	if (world_rank == 0) {
		cout << "num of modes: " << goal << endl;
		cout << "num of ensembles: " << ens << endl;
		cout << "amp of white noise: " << nstd << endl;
		modes1 = new double[SZ*goalt]; 
		rootBuff1 = new double[SZ*goalt];
		cout << "EEMD3D starts!" << endl;
		cout << "=============================================" << endl;
	}

	/* now establish the plan for parallel programing */

	// new plan for storing image data
	int *uCnts0 = new int[world_size];
	int *vCnts0 = new int[world_size];
	int *wCnts0 = new int[world_size];
	int *uDisps0 = new int[world_size];
	int *vDisps0 = new int[world_size];
	int *wDisps0 = new int[world_size];
	uCnts0[0] = toDo(V*W, 0, world_size) * U;
	vCnts0[0] = toDo(U*W, 0, world_size) * V;
	wCnts0[0] = toDo(U*V, 0, world_size) * W;
	uDisps0[0] = 0;
	vDisps0[0] = 0;
	wDisps0[0] = 0;
	for (t = 1; t < world_size; t++) {
		uCnts0[t] = toDo(V*W, t, world_size) * U;
		vCnts0[t] = toDo(U*W, t, world_size) * V;
		wCnts0[t] = toDo(U*V, t, world_size) * W;
		uDisps0[t] = uDisps0[t-1] + uCnts0[t-1];
		vDisps0[t] = vDisps0[t-1] + vCnts0[t-1];
		wDisps0[t] = wDisps0[t-1] + wCnts0[t-1];
	} // end of for-t
	// new plan for gathering post-EEMD data
	int *uCnts1 = new int[world_size];
	int *vCnts1 = new int[world_size];
	int *wCnts1 = new int[world_size];
	int *uDisps1 = new int[world_size];
	int *vDisps1 = new int[world_size];
	int *wDisps1 = new int[world_size];
	uCnts1[0] = toDo(V*W, 0, world_size) * U*goalt;
	vCnts1[0] = toDo(U*W, 0, world_size) * V*goalt;
	wCnts1[0] = toDo(U*V, 0, world_size) * W*goalt;
	uDisps1[0] = 0;
	vDisps1[0] = 0;
	wDisps1[0] = 0;
	for (t = 1; t < world_size; t++) {
		uCnts1[t] = toDo(V*W, t, world_size) * U*goalt;
		vCnts1[t] = toDo(U*W, t, world_size) * V*goalt;
		wCnts1[t] = toDo(U*V, t, world_size) * W*goalt;
		uDisps1[t] = uDisps1[t-1] + uCnts1[t-1];
		vDisps1[t] = vDisps1[t-1] + vCnts1[t-1];
		wDisps1[t] = wDisps1[t-1] + wCnts1[t-1];
	} // end of for-t





	/* ############################################### */
	/* ############ 3D-EEMD core function ############ */
	/* ############################################### */
	




	/*****************************************/
	/* for 1st dimension: solve for each COL */
	/*****************************************/

	int myTodo = toDo(V*W, world_rank, world_size);
	double *myBuff1 = new double[uCnts0[world_rank]];
	double *myModes1 = new double[uCnts0[world_rank]*goalt];
	double *inTmp1 = new double[U];
	double *outTmp1 = new double[U * goalt];
	MPI_Scatterv(img, uCnts0, uDisps0, MPI_DOUBLE, 
		myBuff1, uCnts0[world_rank], MPI_DOUBLE,
		0, MPI_COMM_WORLD);
	for (t = 0; t < myTodo; t++) {
		for(i = 0; i < U; i++) 
			inTmp1[i] = myBuff1[i + t*U];
		eemd(outTmp1, inTmp1, U, goal, ens, nstd);
		for (i = 0; i < U; i++) {
			for (m1 = 0; m1 < goalt; m1++) {
				myModes1[i + m1*U + t*U*goalt]
				 = outTmp1[i + m1*U];
			} // end of for-m1
		} // end of for-i
	} // end of for-t
	delete[] myBuff1, inTmp1, outTmp1;

	// gather and re-sort the data: rootBuff1 => modes1
	MPI_Gatherv(myModes1, uCnts1[world_rank], MPI_DOUBLE, 
		rootBuff1, uCnts1, uDisps1, MPI_DOUBLE,
		0, MPI_COMM_WORLD);
	
	if (world_rank == 0) {
		for (m1 = 0; m1 < goalt; m1++)
			for (i = 0; i < U; i++)
				for (j = 0; j < V; j++)
					for (k = 0; k < W; k++)
						modes1[i + j*U + k*U*V + m1*SZ]
						 = rootBuff1[i + m1*U + j*U*goalt + k*U*goalt*V];				 
	} // end of if
	
	delete[] rootBuff1;
	delete[] img;
	delete[] myModes1;
	delete[] uCnts0, uCnts1, uDisps0, uDisps1;

	t2 = MPI_Wtime();
	dt = t2 - t1;
	if (world_rank == 0) {
		cout << "EEMD stage 1 done in " << dt << "s" << endl;
		cout << "==============================================" << endl;
	}
	
	



	/*****************************************/
	/* for 2nd dimension: solve for each ROW */
	/*****************************************/

	// MPI_Bcast(modes1, SZ*goalt, MPI_DOUBLE, 0, MPI_COMM_WORLD); // seed for 2nd EEMD phase
	myTodo = toDo(U*W, world_rank, world_size);
	// buffers for each layer with mutiple modes
	double *modeBuff2in = NULL; 
	double *myBuff2 = new double[vCnts0[world_rank]];
	double *myModes2 = new double[vCnts0[world_rank]*goalt];
	double *inTmp2 = new double[V];
	double *outTmp2 = new double[V*goalt];
	// buffer for inter-mode data collecting
	double *modeBuff2out = NULL;
	// buffers for post-EEMD data collecting
	double *modes2 = NULL;
	// memory allocation for root-only arrays
	if (world_rank == 0) {
		modeBuff2in = new double[SZ]; 
		modeBuff2out = new double[SZ * goalt];
		modes2 = new double[SZ*goalt*goalt];
	}

	// parallel in each mode
	for (m1 = 0; m1 < goalt; m1++) {
		
		t3 = MPI_Wtime();
		if (world_rank == 0) {
			for (i = 0; i < U; i++)
				for (j = 0; j < V; j++) 
					for (k = 0; k < W; k++)
						modeBuff2in[j + i*V + k*V*U]
						 = modes1[i + j*U + k*U*V + m1*SZ];
		} // end of if

		// scatter and EEMD
		MPI_Scatterv(modeBuff2in, vCnts0, vDisps0, MPI_DOUBLE, 
			myBuff2, vCnts0[world_rank], MPI_DOUBLE,
			0, MPI_COMM_WORLD);
		for (t = 0; t < myTodo; t++) {
			for (j = 0; j < V; j++) 
				inTmp2[j] = myBuff2[j + t*V];
		
			eemd(outTmp2, inTmp2, V, goal, ens, nstd);
			for (j = 0; j < V; j++) {
				for (m2 = 0; m2 < goalt; m2++) {
					myModes2[j + m2*V + t*V*goalt]
					 = outTmp2[j + m2*V];
				} // end of for-m2
			} // end of for-j			
		} // end of for-k

		// gather and re-sort
		MPI_Gatherv(myModes2, vCnts1[world_rank], MPI_DOUBLE, 
			modeBuff2out, vCnts1, vDisps1, MPI_DOUBLE,
			0, MPI_COMM_WORLD);
		if (world_rank == 0) {
			for (m2 = 0; m2 < goalt; m2++)
				for (i = 0; i < U; i++)
					for (j = 0; j < V; j++)
						for (k = 0; k < W; k++)
							modes2[i + j*U + k*U*V + m1*SZ + m2*SZ*goalt]
							 = modeBuff2out[j + m2*V + i*V*goalt + k*V*goalt*U];
		
		// ETA estimation
		t4 = MPI_Wtime();
		dt = t4 - t3;
		eta_time = (t4 - t2) * (goalt-m1-1) / (m1+1);
		cout << "Mode " << m1+1 << "/" << goalt 
			 << " solved in " << dt << "s.  ";
		cout << eta_time << "s to go..." << endl; 
		} // end of if (world_rank == 0)
	} // end of for-m1
	delete[] myBuff2, myModes2, inTmp2, outTmp2;
	delete[] vCnts0, vCnts1, vDisps0, vDisps1;


	if (world_rank == 0) {
		delete[] modes1, modeBuff2in, modeBuff2out;
		dt = t4 - t2;
		cout << "EEMD stage 2 done in " << dt << "s" << endl;
		cout << "==============================================" << endl;
	}





	/*****************************************/
	/* for 3rd dimension: solve for each STK */
	/*****************************************/

	//MPI_Bcast(modes2, SZ*goalt*goalt, MPI_DOUBLE, 0, MPI_COMM_WORLD); // seed for 3nd EEMD phase
	myTodo = toDo(U*V, world_rank, world_size);
	// buffers for each layer with mutiple modes
	double *modeBuff3in = NULL; 
	double *myBuff3 = new double[wCnts0[world_rank]];
	double *myModes3 = new double[wCnts0[world_rank]*goalt];
	double *inTmp3 = new double[W];
	double *outTmp3 = new double[W*goalt];
	// buffer for inter-mode data collecting
	double *modeBuff3out = NULL;
	// buffers for post-EEMD data collecting
	double *modes3 = NULL;
	// memory allocation for root-only arrays`
	if (world_rank == 0) {
		modeBuff3in = new double[SZ]; 
		modeBuff3out = new double[SZ * goalt];
		modes3 = new double[SZ*goalt*goalt*goalt];
	}
	// parallel in each m1 mode in each m2 mode
	for (m2 = 0; m2 < goalt; m2++) {
		for (m1 = 0; m1 < goalt; m1++) {
			
			t5 = MPI_Wtime();
			if (world_rank == 0) {
				for (i = 0; i < U; i++)
					for (j = 0; j < V; j++)
						for (k = 0; k < W; k++)
							modeBuff3in[k + i*W + j*W*U]
							 = modes2[i + j*U + k*U*V + m1*SZ + m2*SZ*goalt];
			} // end of if

			// scatter and EEMD
			MPI_Scatterv(modeBuff3in, wCnts0, wDisps0, MPI_DOUBLE, 
				myBuff3, wCnts0[world_rank], MPI_DOUBLE,
				0, MPI_COMM_WORLD);
			for (t = 0; t < myTodo; t++) {
				for (k= 0; k < W; k++)
					inTmp3[k] = myBuff3[k + t*W];
				eemd(outTmp3, inTmp3, W, goal, ens, nstd);
				for (k = 0; k < W; k++) {
					for (m3 = 0; m3 < goalt; m3++) {
						myModes3[k + m3*W + t*W*goalt]
						 = outTmp3[k + m3*W];
					} // end of for-m3
				} // end of for-k
			} // end of for-j

			// gather and re-sort
			MPI_Gatherv(myModes3, wCnts1[world_rank], MPI_DOUBLE, 
				modeBuff3out, wCnts1, wDisps1, MPI_DOUBLE,
				0, MPI_COMM_WORLD);
			if (world_rank == 0) {
				for (m3 = 0; m3 < goalt; m3++)
					for (i = 0; i < U; i++)
						for (j = 0; j < V; j++)
							for (k = 0; k < W; k++)
								modes3[i + j*U + k*U*V + m1*SZ + m2*SZ*goalt + m3*SZ*goalt*goalt]
								 = modeBuff3out[k + m3*W + i*W*goalt + j*W*goalt*U];
			// ETA estimation
			t6 = MPI_Wtime();
			dt = t6 - t5;
			eta_time
			 = (t6 - t4) * ( (goalt-m1-1) + (goalt-m2)*goalt )
			  / ((m1+1) + m2*goalt);
			cout << "Mode " << m1+1 << "/" << goalt
				 << " in mode " << m2+1 << "/" << goalt
				 << " solved in " << dt << "s.  ";
			cout << eta_time << "s to go..." << endl; 
			} // end of if (world_rank == 0)
		} // end of for-m1
	} // end of for-m2
	delete[] modeBuff3in, myBuff3, myModes3, inTmp3, outTmp3;
	delete[] wCnts0, wCnts1, wDisps0, wDisps1;

	if (world_rank == 0) {
		delete[] modes2, modeBuff3in, modeBuff3out;
		dt = t6 - t4;
		cout << "EEMD stage 3 done in " << dt << "s" << endl;
		cout << "==============================================" << endl;
	}
	




	/*****************************************/
	/*             Combine modes             */
	/*****************************************/

	if (world_rank == 0)
		cout << "Combining modes... " << endl;
	double *modes = NULL;
	if (world_rank == 0) {
		modes = new double[SZ * goalt];
		for (int t = 0; t < SZ*goalt; t++)
			modes[t] = 0;
		for (r = goal; r >= 0; r--) {
			for (k = 0; k < W; k++) {
				for (i = 0; i < U; i++) {
					for (j = 0; j < V; j++) {
						for (m1 = r; m1 < goalt; m1++) {
							for (m2 = r; m2 < goalt; m2++) {
								for (m3 = r; m3 < goalt; m3++) {
									modes[i + j*U + k*U*V + r*SZ]
									 = modes[i + j*U + k*U*V + r*SZ]
									  + modes3[i + j*U + k*U*V + m1*SZ + m2*SZ*goalt + m3*SZ*goalt*goalt];
								} // end of for-m3
			}}}}} // end of for-m2, -m1, -j, -i, -k
			if (r < goal) {
				for (k = 0; k < W; k++) {
					for (i = 0; i < U; i++) {
						for (j = 0; j < V; j++) {
							for (m1 = r+1; m1 < goalt; m1++) {
								modes[i + j*U + k*U*V + r*SZ]
								 = modes[i + j*U + k*U*V + r*SZ]
								  - modes[i + j*U + k*U*V + m1*SZ];
							} // end of for-m1
				}}} // end of for-j, -i, -k
			} // end of if (r < goal)
		} // end of for-r
		delete[] modes3;
	} // end of if (world_rank == 0) 
	

	/* ############ end of 3D-EEMD core function ############ */






	// export the result to file
	t7 = MPI_Wtime();
	if (world_rank == 0) {
		cout << "Done! Writing to file... " << endl;
		int dim2 = 4;
		int lg2[4] = {U, V, W, goalt};
		//int timecode = (int)time(NULL) % 10000;
		//char timecode_str[4];
		//sprintf(timecode_str, "%d", timecode);
		string filenameStr(argv[1]);
		string filename_export // v
		 = string(filenameStr,0,filenameStr.length()-0)+"_modes" + timecode_str + ".bin";
		string filename_log    // v
		 = string(filenameStr,0,filenameStr.length()-0)+"_log" + timecode_str + ".txt";

		// write output file
		writeBinary(filename_export, dim2, lg2, modes);

		cout << filename_export << " exported!" << endl;
		ofstream fout(filename_log.c_str());
		fout << "Input: " << filenameStr << endl;
		fout << "# of modes: " << goal << endl;
		fout << "# of ensembles: " << ens << endl;
		fout << "Amp. of white noise: " << nstd << endl;
		fout << "# of cores: " << world_size << endl;
		fout << "Elapsed time: " << t6 - t1 << "s" << endl;
		fout.close();
		cout << "Elapsed time: " << t6 - t1 << "s" << endl;
		delete[] modes;
	} // end of if (world_rank == 0)
	MPI_Finalize();
	return 0;
} // end of main()

/* toDo: calculate the partition size for parallel work */
int toDo(int N, int myrank, int world_size) {
	int num = floor((N*1.0)/world_size);
	int remaining = N - num * world_size;
	if (myrank >= world_size - remaining)
		num++;
	return num;
}

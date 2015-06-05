# MPI (parallel) version 3-D EEMD
by Po-Nan Li in Institute of Physics, Academia Sinica, Taiwan

# Compilation example
Require GSL (GNU Scientific Library) and MPI
Intel compiler:

```
mpicxx eemd3d_mpi.cpp -lgsl -lgslcblas -lm -O3 -o mpi-eemd3
```

# Running command example

```
mpirun -np 50 ./mpi-eemd3 INPUT.img 3 100 1
```

We will provide the data for test in the repo in the furture! 

## Arguments
- First argument: input file in CSV style with first number indicating the number of dimension (i.e. 3), followed by three integers indicating the size in three dimensions. Remaining values are the elements of the input data in 1-D array style.
- Second argument: numbers of modes to decompose. (Default: 3)
- Third argument: Number of ensembles. (Default: 1)
- Fourth argument: Sigma of the white noise. If 0 (default), the EEMD will degenerate to EMD.
- Fifth argument: (not shown) file name suffix for the output data. Can be any integers in 0~9999. By default a random number will be assigned. 

# File usage
- **eemd3d_mpi_v3.cpp**: main function
- **eemd.cpp**: 1-D EEMD (ensemble empirical mode decomposition) 
- **emd_core.cpp**: 1-D EMD (empirical mode decomposition)
- **find_extrema.cpp**: subfunction for emd_core.cpp
- **spline_gsl.cpp**: subfunction for emd_core.cpp. Powered by GSL.
- **print2bin.cpp**: subfunction for the main function.
- **README.md**: this file.

## Scripts for post-EEMD data processing
- **binsize.m**: show the size (in each dimension) of the exported data.
- **bin2m.m**: convert the binaray data to Matlab array.


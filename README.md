# rhOver

## Installation

1. Prerequisites

   Development files and libraries of LAPACK and BLAS need to be installed to build rhOver.

2. Compilation

   Select your Fortran compiler (gfortran or ifort) and the parallelization (serial, OpenMP) with the corresponding Makefile for your platform (only 64-bit supported):
   
   `ln -s Makefiles/TARGET Makefile`
    
    where TARGET refers to the following available targets:
	   gfortran-serial		GNU gfortran
	   gfortran-omp		GNU gfortran with OpenMP parallelization
	   ifort-serial		Intel Fortran Compiler
	   ifort-omp		Intel Fortran Compiler with OpenMP parallelization

    rhOver was tested with the following compilers:
    * gfortran 4.8.5
    * gfortran 5.2.0
    * ifort version 12.1.6

    The following command compiles rhOver and creates the executable './rhover' in the source directory:
    `make`

3. Run rhOver

    A variety of sample input files for rhOver is given in the Examples/ folder.
    The rhOver job can be started from the command line with the following command:

    `/path/to/rhover input.dy3 > output.log &`

    Note: 
    The stack limit should be set to 'unlimited' by running the command `ulimit -s unlimited` before running rhOver, otherwise it can happen that rhOver crashes during matrix algebra operations.

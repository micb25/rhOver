[![License](https://img.shields.io/github/license/micb25/rhOver.svg)](LICENSE)
[![Issues](https://img.shields.io/github/issues/micb25/rhOver.svg)](https://github.com/micb25/rhOver/issues)

# rhOver

## Installation

1. Prerequisites

   Development files and libraries of LAPACK and BLAS need to be installed to build rhOver succesfully.

2. Downloading

   The latest rhOver source code can be obtained by the following command:
   
   `git clone https://github.com/micb25/rhOver`
   
   Enter the newly created subdirectory "rhOver" by typing: 
   
   `cd rhOver`

   In addition, rhOver utilizes [libboys](https://github.com/micb25/libboys) to calculate molecular integrals. Hence, the corresponding files need to be obtained via the following command:
   
   `git submodule update --init --recursive`
   
3. Compilation

   Select your Fortran compiler (`gfortran` or `ifort`) and the parallelization (serial, OpenMP) with the corresponding Makefile for your platform (only 64-bit supported):
   
   `ln -s Makefiles/TARGET Makefile`
    
    where TARGET refers to the following available targets:
    
	   gfortran-serial		GNU gfortran
	   gfortran-omp		GNU gfortran with OpenMP parallelization
	   ifort-serial		Intel Fortran Compiler
	   ifort-omp		Intel Fortran Compiler with OpenMP parallelization

    rhOver was tested with the following compilers:
    * gfortran 4.4.6
    * gfortran 4.8.5
    * gfortran 5.2.0
    * gfortran 7.3.1
    * gfortran 7.4.1
    * ifort version 12.1.6

    The following command compiles rhOver and creates the executable './rhover' in the source directory:
    
    `make`

4. Run rhOver

    A variety of sample input files for rhOver is given in the Examples/ folder.
    The rhOver job can be started from the command line with the following command:

    `/path/to/rhover input.dy3 > output.log`

    Note: 
    The stack limit should be set to 'unlimited' by running the command `ulimit -s unlimited` before running rhOver, otherwise it can happen that rhOver crashes during matrix algebra operations.

## Citation

Please quote the usage of the rhOver program and any scientific results in any form obtained with it by the following reference:

*[rhOver: Determination of Magnetic Anisotropy and Related Properties for Dysprosium(III) Single-Ion Magnets by Semi-Empirical Approaches utilizing Hartree-Fock Wave Functions.](https://dx.doi.org/10.1002/jcc.25565)* Michael Böhme, Winfried Plass, *J. Comput. Chem.* **2018**, *39* (32), 2697–2712. *DOI: 10.1002/jcc.25565*

!
! rhOver - a FORTRAN program to determine magnetic anisotropy and related 
! properties for dysprosium(III) single-ion magnets by semi-empirical approaches
! Copyright (C) 2014-2019 Michael Böhme <boehme.mic@gmail.com>
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!

module global_c

     ! constants
     double precision, parameter :: au2pm = 5.2917721E-01
     double precision, parameter :: au2kJpmol = 2625.49953026d0
     double precision, parameter :: au2rcm = 219474.63633664d0
     double precision, parameter :: au2K = 315774.62979379d0
     double precision, parameter :: Pi = 3.1415926535897932d0
     double precision, parameter :: TwoPi = 2d0 * Pi
     double precision, parameter :: FourPi = 4d0 * Pi
     double precision, parameter :: SqrtPi = 1.7724538509055160
     integer, parameter :: sMaxBuffer = 256
     integer, parameter :: iMaxAtoms = 256
     integer, parameter :: iMaxCGTOs = 10000
     integer, parameter :: iMaxMOs = 10000
     double precision :: TotalNuclearCharge, LigTotalCharge, TotalElCharge
     
     character (len=sMaxBuffer) :: XFilenameMin, DEPPFile
     character (len=sMaxBuffer) :: OLFile
     character (len=sMaxBuffer) :: MOFile
     character (len=sMaxBuffer) :: JobTitle, SHostname, SWorkDir
     character (len=sMaxBuffer) :: GPFile
     
     double precision :: SphereRad = 1d0 / au2pm
     double precision :: DistMin = -0.01d0
     double precision :: QdivRThres = 0.0d0
     
     ! optimal parameters for Freeman/Watson
!      double precision :: Ln4fCutOff_I = 0.03d0
!      double precision :: Ln4fCutOff_O = 3.50d0
     ! optimal parameters for ANO-RCC
     double precision :: Ln4fCutOff_I = 0.02d0
     double precision :: Ln4fCutOff_O = 4.00d0
!      double precision :: Ln4fCutOff_O = 8.00d0
     double precision :: LFTCutOff_I = 0.00d0
     double precision :: LFTCutOff_O = 3.50d0
!      double precision :: LFTCutOff_O = 8.00d0
     double precision, parameter :: Ahlrichs_Alpha = 0.6d0
     integer :: nGrid
     
     ! mJ state: 
     !  1 = +/- 15/2, 
     !  2 = +/- 13/2, 
     !  3 = +/- 11/2,
     !  4 = +/-  9/2,
     !  5 = +/-  7/2,
     !  6 = +/-  5/2,
     !  7 = +/-  3/2,
     !  8 = +/-  1/2
     integer :: mJ = 1
     
     ! radial wave function: 
     !  1 = ANO-RCC 4f basis set for Dy(0) (as used in the original JCC article), 
     !  2 = Freeman-Watson 1962 for Dy(III)
     !  3 = CASSCF/DKH2/ANO-RCC-based Dy(III) (M. Böhme 2019, unpublished)
     integer :: iRadWF = 1
     
     ! file units
     integer, parameter :: uEner = 12
     integer, parameter :: uInp  = 13
     integer, parameter :: uXYZ  = 14
     integer, parameter :: uMolF = 15
     integer, parameter :: uPotF = 16
     
     integer :: NRad, NLigSpher, NLigGP, NGP
     integer :: NLigRadShells
     ! global program options
     logical :: OVerbose, OGenPotE, OMagVec, OGenXYZ, OAngstrom, OMax, OTraj, OPruning, ORandRot, ONoDEPPFile
     logical :: OOldRInts, OPrintBasis, OPrintMulliken, OCFOnly, OIESMode, ODelete4f, OExpert
     logical :: OSkipMom, OLDAX, OMixedMode, ODeleteAll
     logical :: OPCM, OSternheimer
     logical :: OExportGPs
     logical :: OUHF
     
     ! LAPACK
     logical :: LAPACKInit = .FALSE.
     integer :: NArraySize, LWORKSize, RWORKSize, NArraySizeLDA, lapackstatus
     complex*16, dimension(16, 16) :: LFMatrix
     complex*16, dimension(:), allocatable :: LWorkArr
     double precision, dimension(:), allocatable :: RWorkArr
     double precision, dimension(16) :: EigenValues
     character :: LAPACKJobType, LAPACKMatrixType
     
     ! Ligand-field theory
     logical :: CFInit = .FALSE.
     double precision, dimension(0:6) :: StevMultFac, OldRadInts, SternheimerShieldings
     double complex, dimension(-2:2,1:16, 1:16) :: StevOp2Mat, HermOp2O, HermOp2W, TensorOp2Mat
     double complex, dimension(-4:4,1:16, 1:16) :: StevOp4Mat, HermOp4O, HermOp4W, TensorOp4Mat
     double complex, dimension(-6:6,1:16, 1:16) :: StevOp6Mat, HermOp6O, HermOp6W, TensorOp6Mat
     double precision, dimension(6, 0:6) :: StevPLM, StevLam
     double precision, dimension(0:15, 16) :: WFDecomp
     double precision, dimension(0:6,-6:6) :: ARkq, BRkq
     double complex, dimension(0:6,-6:6) :: BWkq
     integer :: MaxKRank
     double precision, dimension(0:6) :: LFPScalingFactors
     ! mJ coefficients
     double precision, dimension(8,4) :: coeff
     
     ! basis set
     integer :: NPGTO = 0
     integer :: NumPGTO, NumCGTO, NumDelPGTO
     integer :: NCustomCharges = 0
     
     ! origin of Dy(III)
     double precision :: EDyX, EDyY, EDyZ    
     integer :: CAtomID = 0
     
     double precision :: Energy_SUM
     double precision :: TotDensMO, TotDensDS, TotNormWF
     double precision :: LinScanX1, LinScanY1, LinScanX2, LinScanY2, LinScanDX, LinScanDY
     integer :: LinScanS
     
     double precision :: NumberOf4fElectrons = 9d0
     
     ! Stevens-factors A_2, A_4, A_6 for Dy(III)
     ! J. Sievers, Z. Phys. B., 1982 (45), 289--296
     double precision, dimension(8,3) :: StevensFactors_Dy = reshape( &
     (/ & 
       (/ -0.33333333333E0, -0.200000000000E0, -0.08571428571E0 /), & ! 15/2
       (/  0.00952380952E0,  0.085714285710E0,  0.14285714286E0 /), & ! 13/2
       (/  0.18095238095E0,  0.200000000000E0, -0.12121212121E0 /), & ! 11/2
       (/  0.04040404040E0,  0.098124098120E0,  0.08924408924E0 /), & !  9/2
       (/  0.04484404484E0, -0.010212010210E0, -0.05727605728E0 /), & !  7/2
       (/ -0.08391608392E0,  0.058275058280E0, -0.10489510490E0 /), & !  5/2
       (/ -0.03496503497E0,  0.052895822130E0,  0.07799892415E0 /), & !  3/2
       (/  0.04034427111E0, -0.02241348395E0, -0.067240451856E0 /)  & !  1/2
     /), shape(StevensFactors_Dy), order=(/2, 1/) )
     
#ifdef _OPENMP
     integer :: omp_threads, omp_man_threads
#endif
     
end module global_c

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

module data_c
	use global_c
	implicit none

	double precision :: EnergyConv, AngConv
	double precision :: MagDist
	! 1 = single, 2 = screening, 3 = find min
	integer :: JobType, nStep, MaxIter	
	double precision :: incStep
	double precision :: cStates(8)

	double precision :: MagX, MagY, MagZ, fMRotX, fMRotY
	double precision :: rho
	double precision :: MagVecDev
	
	! rotation
	double precision :: fRotA, fRotB, fRotG
	double precision :: minRotX, minRotY, minValue
	double precision :: maxRotX, maxRotY, maxValue
	double precision :: minScanX, maxScanX, minScanY, maxScanY, minScanZ, maxScanZ, incScanX, incScanY, incScanZ

	character (len=sMaxBuffer) :: sLine, InpFile, SepLine, TRJFilename

end module

subroutine print_jobinfo
	use data_c
	use data_mo
#ifdef _OPENMP
	use omp_lib
#endif
	implicit none
!
	write(*,*) "> JOB INFO:"
	write(*,*) 

	write(*,'(A35,A)') "title: ", trim(JobTitle)
	call Hostnm(SHostname)
	write(*,'(A35,A)') "host name: ", trim(SHostname)
	call GetCWD(SWorkDir)
	write(*,'(A35,A)') "working directory: ", trim(SWorkDir)
	
	write(*,'(A35)', advance='no') "type of job: "
	select case ( JobType )
	      case ( 1 )
			write(*,'(A)') "Single-Point"
	      case ( 2 )
			write(*,'(A)') "SCAN"
	      case ( 3 )
			write(*,'(A)') "OPT"
	      case ( 4 )
			write(*,'(A)') "SCANOPT"
	      case ( 5 )
			write(*,'(A)') "LINSCAN"
	      case default
			write(*,'(A)') "unknown"
	end select  
	
	if ( JobType .eq. 4 ) then
		write(*,'(A35)',advance='no') "mixed ESC/LFT mode: "
		call write_option_bool(OMixedMode)
	end if
	
#ifdef _OPENMP
	write(*,'(A35,I14)') "parallel OpenMP threads: ", omp_get_max_threads()
#endif

	write(*,*)
	if ( CAtomID .eq. 0 ) then
		write(*,'(A35,F14.8,F14.8,F14.8)') "Dy(III) ion is located at (a.u.): ", EDyX, EDyY, EDyZ
		write(*,'(A35,F14.8,F14.8,F14.8)') "(angstrom): ", EDyX*au2pm, EDyY*au2pm, EDyZ*au2pm
	else
		write(*,'(A35,I14)') "Dy(III) ion has atom id: ", CAtomID
	end if

	write(*,*)
	write(*,'(A35,A8,I2,A)') "m_J state: ", "  | +/- ", int(15-(mJ-1)*2) , "/2 >"
	write(*,*)
	
	if ( OMagVec .eqv. .TRUE. ) then
	      write(*,'(A35,3F14.8)') "reference axis: ", MagX, MagY, MagZ
	end if

	write(*,*)
	write(*,*)

end subroutine

!---------------------------------------------------
!       write_xyz_file
!       generates structure file with the anisotropy axis (two He atoms)
!---------------------------------------------------
subroutine write_xyz_file(rotx, roty, trj, val)
	use data_c
	use data_mo
	use data_grid
	implicit none
	double precision, intent(in) :: rotx, roty
	logical, intent(in) :: trj
	double precision, intent(in) :: val
	double precision :: RDyX, RDyY, RDyZ    ! rotated coordinates
	integer :: i, iost
	character (len=sMaxBuffer) :: sBuffer
	
	if ( trj .eqv. .FALSE. ) then
		write(*,*)
		write(*,*)
		write(*,*) "> GENERATING STRUCTURE WITH QUANTIZATION AXIS ..."
		write(*,*)
		write(*,'(5X,A)') "writing XYZ file '" // trim(XFilenameMin) // "' ..."
		open(unit=uXYZ, file=XFilenameMin,action='write', status='replace',form='formatted')
	else
		open(unit=uXYZ, file=TRJFilename,action='write', status='old', form='formatted', position='append',iostat=iost)
		if ( iost .ne. 0 ) then
			open(unit=uXYZ, file=TRJFilename,action='write', status='new', form='formatted')
		end if 
	end if

	! file header
	if ( OPCM .eqv. .FALSE. ) then
		write(sBuffer, '(I8)') NA + 3
	else
		write(sBuffer, '(I8)') NCustomCharges + 3
	end if
	write(uXYZ, '(A)') trim(adjustl(sBuffer))
	
	if ( trj .eqv. .FALSE. ) then
		write(uXYZ, '(A/)', advance='no') "Generated with rhOver - (C) 2018 by Michael Böhme"
	else
		write(uXYZ, '(A,F12.6/)', advance='no') "Overlap ", val
	end if
	
	! Dy(III) 
	write(uXYZ,'(A,3F12.6)') adjustl(trim("Dy")), EDyX*au2pm, EDyY*au2pm, EDyZ*au2pm
	! add two He atoms to mark the main magnetic anisotropy axis
	call rotate_euler_3d(rotx, roty, 0d0, 0d0, 0d0, 0d0 + MagDist, RDyX, RDyY, RDyZ)
	write(uXYZ,'(A,3F12.6)') adjustl(trim("Xx")),  (EDyX+RDyX)*au2pm, (EDyY+RDyY)*au2pm, (EDyZ+RDyZ)*au2pm
	write(uXYZ,'(A,3F12.6)') adjustl(trim("Xx")),  (EDyX-RDyX)*au2pm, (EDyY-RDyY)*au2pm, (EDyZ-RDyZ)*au2pm

	! other atoms
	if ( OPCM .eqv. .FALSE. ) then
		do i = 1, NA

			if ( Atoms(i)%ghost .eqv. .FALSE. ) then
				call get_Element_by_Charge(Atoms(i)%charge, sBuffer)
			else
				call get_Element_by_Charge(0, sBuffer)
			end if
			
			RDyX = Atoms(i)%x-EDyX
			RDyY = Atoms(i)%y-EDyY
			RDyZ = Atoms(i)%z-EDyZ
			write(uXYZ,'(A,3F16.8)') adjustl(trim(sBuffer)), (RDyX+EDyX)*au2pm, (RDyY+EDyY)*au2pm, (RDyZ+EDyZ)*au2pm

		end do
	else
		do i = 1, NCustomCharges
			RDyX = PointCharges(i)%x-EDyX
			RDyY = PointCharges(i)%y-EDyY
			RDyZ = PointCharges(i)%z-EDyZ
			write(uXYZ,'(A,3F16.8)') adjustl(trim("Xx")), (RDyX+EDyX)*au2pm, (RDyY+EDyY)*au2pm, (RDyZ+EDyZ)*au2pm
		end do
	end if

	close(uXYZ)
	
	if ( trj .eqv. .FALSE. ) then
		call write_done
	end if
	

end subroutine

subroutine parse_input_file(filename)
	use data_c
	use data_mo 
	use data_grid
#ifdef _OPENMP
	use omp_lib
#endif	
	implicit none
!
	character (len=sMaxBuffer), intent(in) :: filename
	integer :: iost, i = 0, j
	logical :: EOF = .FALSE.
	logical :: FMagicStr = .FALSE.
	logical :: FMOFile   = .FALSE.
	logical :: FCoord    = .FALSE.
	logical :: FLinScan  = .FALSE.
	double precision :: tempVal, tempValX, tempValY, tempValZ
	character (len=sMaxBuffer) :: to_upper
!
	write(*,'(/X,A,A,A/)', advance='no') "PARSING INPUT FILE '", trim(filename), "' ..."
	
	open(unit=uInp, file=filename, action='read',status='old', iostat=iost)
	if ( iost .ne. 0 ) then
		write(*,*) "ERROR! Can't read input file!"
		stop
	end if

	do while ( EOF .eqv. .FALSE. )
	    
		i = i + 1
		read(uInp, '(A255)', iostat=iost) sLine
		    
		sLine = trim(adjustl(sLine))
		if ( iost .lt. 0 ) then
			EOF = .TRUE.
		end if	
		    
		! skip empty lines
		if ( len_trim(sLine) .ne. 0 ) then

			sLine = to_upper(sLine)
			  
			! parse input strings
			if ( ( sLine(1:1) == "*" ) .OR. ( sLine(1:1) == "!" ) .OR. ( sLine(1:1) == "#" ) ) then
			  
				! comment, do nothing
			  
			else if ( sLine == "RHOVER" ) then
			  
				FMagicStr = .TRUE.
				
			else if ( sLine == "TITLE" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				JobTitle = trim(adjustl(sLine))
				
			else if ( sLine == "JOB" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(to_upper(sLine)))
				if ( sLine == "SCANOPT" ) then
				    JobType = 4
				else if ( sLine == "OPT" ) then
				    JobType = 3
				else if ( sLine == "SCAN" ) then
				    JobType = 2
				else if ( sLine == "SP" ) then
				    JobType = 1
				else if ( sLine == "LINSCAN" ) then
				    JobType = 5
				else
					write(*,*) "ERROR! Unknown job type!"
					stop
				end if
				
			else if ( sLine == "MAXITER" ) then
			  
				i = i + 1
				read(uInp, *, iostat=iost) MaxIter
				if ( MaxIter .lt. 1 ) then
					write(*,*) "ERROR! MaxIter value is invalid!"
					stop
				end if
				
			else if ( sLine == "MOFILE" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				MOFile = trim(adjustl(sLine))
				FMOFile = .TRUE.
				
			else if ( sLine == "DYIII" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) EDyX, EDyY, EDyZ
				FCoord = .TRUE.
				
			else if ( sLine == "PRINTBASIS" ) then
			
				OPrintBasis = .TRUE.
				
			else if ( sLine == "MULLIKEN" ) then
			
				OPrintMulliken = .TRUE.
				
			else if ( ( sLine == "LDAX" ) .or. ( sLine == "SLATERX" ) ) then
			  
				OLDAX = .TRUE.
				
			else if ( ( sLine == "NOLFT" ) .or. ( sLine == "NOCFT" ) ) then
			  
				OIESMode = .TRUE.
				OMixedMode = .FALSE.
				OMax = .FALSE.
			
			else if ( ( sLine == "LFTONLY" ) .or. ( sLine == "CFTONLY" ) ) then
			  
				OIESMode = .FALSE.
				OMixedMode = .FALSE.
				OMax = .TRUE.
				
			else if ( sLine == "ATOMID" ) then
			
				i = i + 1
				read(uInp, *) CAtomID
				if ( ( CAtomID .lt. 1 ) .OR. ( CAtomID .gt. iMaxAtoms ) ) then
					write(*,*) "ERROR! Invalid atom id specified!"
					stop
				end if
				
			else if ( sLine == "REFERENCE" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) MagX, MagY, MagZ
				OMagVec = .TRUE.
				
				fMRotX = dasin( MagY )
				fMRotY = dasin( - MagX / dcos(fMRotX) )
				
				if ( fMRotX .lt. 0.0d0 ) then
					fMRotX = fMRotX + pi
				end if

			else if ( sLine == "SETROT" ) then
				
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) fRotA, fRotB
				fRotA = fRotA * pi/180d0
				fRotB = fRotB * pi/180d0
				fRotB = 0d0
				
			else if ( sLine == "SETROT3D" ) then
				
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) fRotA, fRotB, fRotG
				fRotA = fRotA * pi/180d0
				fRotB = fRotB * pi/180d0
				fRotG = fRotG * pi/180d0
				
			else if ( sLine == "DANGLE" ) then

				i = i + 1
				read(uInp, *, iostat=iost) incStep
				if ( ( incStep .lt. 0.001 ) .OR. ( incStep .gt. 180.0 ) ) then
					write(*,*) "ERROR! DANGLE option is invalid!"
					stop
				end if 
				incStep = incStep * 2.0 * pi / 360.0
				
			else if ( sLine == "GRID" ) then
			  
				i = i + 1
				read(uInp, *, iostat=iost) nGrid
				if ( ( nGrid .lt. 1 ) .OR. ( nGrid .gt. 5 ) )  then
					write(*,*) "ERROR! Invalid value for 'GRID'!"
					stop
				end if
			
			else if ( sLine == "SAVEGRIDPOINTS" ) then
			
				OExportGPs = .TRUE.
				
			else if ( sLine == "LFTCUTOFF" ) then
			
				i = i + 1
				read(uInp, *, iostat=iost) LFTCutOff_O
				if ( LFTCutOff_O .le. LFTCutOff_I ) then
					write(*,*) "ERROR! Invalid value for 'LFTCUTOFF'!"
					stop
				end if
				
			else if ( sLine == "PARA" ) then
			
#ifdef _OPENMP
				i = i + 1
				read(uInp, *, iostat=iost) omp_man_threads
				if ( omp_man_threads .lt. 1 )  then
					write(*,*) "ERROR! Invalid value for 'PARA'!"
					stop
				end if
#else
				write(*,*) "ERROR! Keyword 'PARA' is not available in serial mode!"
				stop
#endif
			else if ( sLine == "SCANRANGE" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) minScanX, minScanY, maxScanX, maxScanY
				if ( minScanX .gt. maxScanX ) then
					tempVal = maxScanX
					maxScanX = minScanX
					minScanX = tempVal
				end if
				if ( minScanY .gt. maxScanY ) then
					tempVal = maxScanY
					maxScanY = minScanY
					minScanY = tempVal
				end if
				
			else if ( sLine == "LINSCANRANGE" ) then
				
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) LinScanS, LinScanX1, LinScanY1, LinScanX2, LinScanY2
				
				if ( LinScanS .lt. 2 ) then
					write(*,*) "ERROR! Number of points is invalid!"
					stop
				end if
				
				LinScanDX = (LinScanX2-LinScanX1) / (LinScanS-1)
				LinScanDY = (LinScanY2-LinScanY1) / (LinScanS-1)
				FLinScan = .TRUE.
				
			else if ( sLine == "SCANINC" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				sLine = trim(adjustl(sLine))
				read(sLine, *) incScanX, incScanY
				if ( incScanX .lt. 0.0 ) then
					incScanX = - incScanX
				end if
				if ( incScanY .lt. 0.0 ) then
					incScanY = - incScanY
				end if
				if ( ( incScanX .eq. 0.0 ) .OR. ( incScanY .eq. 0.0 ) ) then
					write(*,*) "ERROR! Value of SCANINC is invalid!"
					stop
				end if

			else if ( sLine == "NOGRIDPRUNING" ) then
			  
				OPruning = .FALSE.
				
			else if ( sLine == "GRIDPRUNING" ) then
			  
				OPruning = .TRUE.
				
			else if ( sLine == "NORANDOMROT" ) then
			  
				ORandRot = .FALSE.
				
			else if ( sLine == "RANDOMROT" ) then
			  
				ORandRot = .TRUE.
				
			else if ( sLine == "EXPERT" ) then
			  
				OExpert = .TRUE.
			
			else if ( sLine == "POTCORR" ) then
			
                                        OPotCorr = .TRUE.
			
			else if ( sLine == "NODEPPFILE" ) then
			  
				ONoDEPPFile = .TRUE.

			else if ( sLine == "NOSHIELDING" ) then
			  
				OSternheimer = .FALSE.
				
			else if ( sLine == "SHIELDING" ) then
			  
				OSternheimer = .TRUE.
				
			else if ( sLine == "PCM" ) then
			
				i = i + 1
				
				if ( OPCM .eqv. .TRUE. ) then
					write(*,*) "ERROR! Point charges were already defined!"
					stop
				end if
				
				OPCM = .TRUE.
                                        
				read(uInp, *, iostat=iost) NCustomCharges
				if ( ( NCustomCharges .lt. 1 ) .or. ( NCustomCharges .gt. iMaxAtoms ) ) then
					write(*,*) "ERROR! Too many charges (PCM)!"
					stop
				end if
				
				allocate(PointCharges(NCustomCharges))
				
				do j = 1, NCustomCharges
				
					i = i + 1
					read(uInp,*, iostat=iost) tempValX, tempValY, tempValZ, tempVal
					
					PointCharges(j)%id = j
					PointCharges(j)%x = tempValX
					PointCharges(j)%y = tempValY
					PointCharges(j)%z = tempValZ
					PointCharges(j)%Charge = tempVal
					
				end do
				
			else if ( sLine == "SCALING" ) then
			
				i = i + 1
				
				i = i + 1
				read(uInp,*, iostat=iost) LFPScalingFactors(2), LFPScalingFactors(4), LFPScalingFactors(6)
			
			else if ( sLine == "MJ" ) then
			  
				i = i + 1
				read(uInp, *, iostat=iost) mJ
				if ( ( mJ .lt. 1 ) .OR. ( mJ .gt. 8 ) )  then
					write(*,*) "ERROR! mJ value is invalid!"
					stop
				end if			
				
			else if ( sLine == "MAXGP" ) then
			  
				i = i + 1
				read(uInp, *, iostat=iost) NGP
				if ( ( NGP .lt. 1000 ) .or. ( NGP .gt. 10000000 ) ) then
					write(*,*) "ERROR! MaxGP value is invalid!"
					stop
				end if
				
			else if ( sLine == "TRJFILE" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				TRJFilename = trim(adjustl(sLine))
				OTraj = .TRUE.			
				
			  else if ( sLine == "OLFILE" ) then
			  
				i = i + 1
				read(uInp, '(A)', iostat=iost) sLine
				OLFile = trim(adjustl(sLine))
				OGenPotE = .TRUE.
				
			  else if ( sLine == "ANGSTROM" ) then
			  
				OAngstrom = .TRUE.
				
			  else if ( sLine == "DELETE4F" ) then
			  
				ODelete4f = .TRUE.
				
			  else if ( sLine == "DELETEALL" ) then
			  
				ODeleteAll = .TRUE. 
				
			  else if ( sLine == "NOXYZ" ) then
			  
				OGenXYZ = .FALSE.
							
			  else if ( sLine == "VERBOSE" ) then
			  
				OVerbose = .TRUE.
				
			  else if ( sLine == "END" ) then
			  
				EOF = .TRUE.
				
			  else
			  
				write(*,*)
				write(*,'(X,A,A,A,I4,A)') "ERROR! Unknown keyword '",trim(sLine),"' found on line ", i, " !"
				stop
				
			  end if

		end if
		    
		if ( iost .lt. 0 ) then
			EOF = .TRUE.
		end if
	    
	end do
	
	close(uINP)     
	
	XFilenameMin = trim(adjustl(InpFile))//".xyz"
	OLFile = trim(adjustl(InpFile))//".dat"  
	DEPPFile = trim(adjustl(InpFile))//".depp"
	GPFile = trim(adjustl(InpFile))//".grid.dat"
	
	! convert angstrom to a.u.
	if ( OAngstrom .eqv. .TRUE. ) then
		EDyX = EDyX / au2pm
		EDyY = EDyY / au2pm
		EDyZ = EDyZ / au2pm
		
		if ( OPCM .eqv. .TRUE. ) then
			do i = 1, NCustomCharges
				PointCharges(i)%x = PointCharges(i)%x / au2pm
				PointCharges(i)%y = PointCharges(i)%y / au2pm
				PointCharges(i)%z = PointCharges(i)%z / au2pm
			end do
		end if
		
		OAngstrom = .FALSE.
	end if
	
	! prepare point charges, if necessary
	if ( OPCM .eqv. .TRUE. ) then
		do i = 1, NCustomCharges
			PointCharges(i)%x = PointCharges(i)%x - EDyX
			PointCharges(i)%y = PointCharges(i)%y - EDyY
			PointCharges(i)%z = PointCharges(i)%z - EDyZ
			PointCharges(i)%dist = dsqrt(PointCharges(i)%x**2 + PointCharges(i)%y**2 + PointCharges(i)%z**2)
			
			write(*,'(5X,I3,5F18.8)') i, PointCharges(i)%x, PointCharges(i)%y, PointCharges(i)%z, PointCharges(i)%dist, PointCharges(i)%Charge
		end do
	end if
	
	if      ( FMagicStr .eqv. .FALSE. ) then
		write(*,*) "ERROR! Input file is invalid! (KEYWORD 'RHOVER' IS MISSING!)"
		stop   
	else if ( ( FMOFile .eqv. .FALSE. ) .AND. ( OPCM .eqv. .FALSE. ) ) then
		write(*,*) "ERROR! Neither MOFile nor PCM keyword not found!"
		stop
	else if ( ( OLDAX .eqv. .TRUE. ) .AND. ( OPCM .eqv. .TRUE. ) ) then
		write(*,*) "ERROR! Keywords 'LDAX' and 'PCM' are incompatible!"
		stop
	else if ( ( FMOFile .eqv. .TRUE. ) .AND. ( OPCM .eqv. .TRUE. ) ) then
		write(*,*) "ERROR! Keywords 'MOFile' and 'PCM' are incompatible!"
		stop
	else if ( ( FCoord .eqv. .FALSE. ) .AND. ( CAtomID .eq. 0 ) ) then
		write(*,*) "ERROR! Neither DyIII nor AtomID keyword found!"
		stop   
	else if ( ( JobType .eq. 5 ) .and. ( FLinScan .eqv. .FALSE. ) ) then
		write(*,*) "ERROR! LINSCAN mode activated, but no parameters specified!"
		stop   
	end if
	
#ifdef _OPENMP
	if ( omp_man_threads .gt. 0 ) then
		call omp_set_num_threads(omp_man_threads)
	end if
#endif
	
	call write_done
	write(*,*)

end subroutine

subroutine set_default_values
	use data_c
	use data_mo
	use global_c
#ifdef _OPENMP  
	use omp_lib
#endif
	implicit none
!
	! default values for OPT
	incStep       = 10d0 * pi / 180d0
	fRotA         = 90d0 * pi / 180d0
	fRotB         = 90d0 * pi / 180d0
	fRotG         = 90d0 * pi / 180d0
	EnergyConv    = 1.0E-03
	AngConv       = 5.0E-02
	MaxIter       = 999
	
	! default values for SCAN
	minScanX      =     0d0
	minScanY      =   -90d0
	minScanZ      =     0d0
	maxScanX      =   180d0
	maxScanY      =    90d0
	maxScanZ      =   360d0
	incScanX      =    10d0
	incScanY      =    10d0
	incScanZ      =    10d0
	
	! misc
	minValue      = 1.0E+09 
	minRotX       =     0d0
	minRotY       =     0d0
	maxValue      =-1.0E+09
	maxRotX       =     0d0
	maxRotY       =     0d0
	MagDist       =     5d0
	JobType       =       4
	NLigRadShells =      55
	NGP           =  100000
	nGrid         =       3
	MaxKRank      =       6
	
	LFPScalingFactors = 0d0
	LFPScalingFactors(2) = 1d0
	LFPScalingFactors(4) = 1d0
	LFPScalingFactors(6) = 1d0

	OVerbose      = .FALSE.
	OMagVec       = .FALSE.
	OAngstrom     = .FALSE.
	OMax          = .FALSE.
	OTraj         = .FALSE.
	OGenPotE      = .TRUE.
	OGenXYZ       = .TRUE.
	OPruning      = .TRUE.
	ORandRot      = .FALSE.
	ONoDEPPFile   = .FALSE.
	OPrintBasis   = .FALSE.
	OPrintMulliken= .FALSE.
	OIESMode      = .TRUE.
	ODelete4f     = .FALSE.
	OExpert       = .FALSE.
	OSkipMom      = .TRUE.
	OLDAX         = .FALSE.
	OMixedMode    = .TRUE.
	OPotCorr      = .FALSE.
	OPCM          = .FALSE.
	OSternheimer  = .FALSE.
	OExportGPs    = .FALSE.
	
#ifdef _OPENMP
	omp_threads   = omp_get_max_threads()
	omp_man_threads = 0
#endif
	
end subroutine

subroutine do_scan
	use data_c
	use data_mo
	use data_grid
	implicit none
!
	double precision :: angX, angY, angZ
	double precision :: rhoVal, incScanZOld
!
	angX = minScanX
	angY = minScanY
	angZ = minScanZ

	minScanZ = 0d0
	incScanZOld = incScanZ 
	incScanZ = 360d0
	maxScanZ = 360d0
	write(*,'(/X,A/)') "PERFORMING SCAN OF EULER ANGLES ALPHA AND BETA ..."
	write(*,'(A33,A/)') "energies determined by: ", "electrostatic calculation"
	
	write(*,'(A33,F7.2,A,F7.2,A,F6.2,A)') "alpha scan range: ", minScanX, " deg -- ", maxScanX, " deg (inc. ", incScanX, " deg)"
	write(*,'(A33,F7.2,A,F7.2,A,F6.2,A)') " beta scan range: ", minScanY, " deg -- ", maxScanY, " deg (inc. ", incScanY, " deg)"
	write(*,*)

	if ( OGenPotE .eqv. .TRUE. ) then
		open(unit=uEner, file=OLFile, action='write',status='replace')
	end if
	
	if ( OIESMode .eqv. .FALSE. ) then
		call write_rline(51)
		write(*,"(30X,2A12,A25)") "alpha", "beta", "m_J contribution"
		write(*,"(30X,2A12,A25)") "/ degree", "/ degree", "      "
		call write_rline(51)
	else
		call write_rline(51)
		write(*,"(30X,2A12,A25)") "alpha", "beta", "relative energy"
		write(*,"(30X,2A12,A25)") "/ degree", "/ degree", "/ a.u."
		call write_rline(51)
	end if

	do while ( angX .le. maxScanX )
		fRotA = angX * pi / 180d0
		if ( angX .ne. minScanX ) then
			write(*,*)
		end if
		do while ( angY .le. maxScanY )
			fRotB = angY * pi / 180d0
			do while ( angZ .lt. maxScanZ ) 
				fRotG = angZ * pi / 180d0
				
				if ( OIESMode .eqv. .FALSE. ) then
					call calc_cf_energy(fRotA, fRotB, fRotG, rhoVal)
					write(*,"(30X,2F12.2,F25.8)") angX, angY, rhoVal
				else
					call cpot_calc_energy(fRotA, fRotB, 0d0, rhoVal)
					write(*,"(30X,2F12.2,F25.8)") angX, angY, rhoVal
				end if
				
				if ( OTraj .eqv. .TRUE. ) then
					call write_xyz_file(fRotA, fRotB, .TRUE., rhoVal) 
				end if 
				
				if ( rhoVal .lt. minValue ) then
					minValue = rhoVal
					minRotX = fRotA		  
					minRotY = fRotB
				end if
				if ( rhoVal .gt. maxValue ) then
					maxValue = rhoVal
					maxRotX = fRotA		  
					maxRotY = fRotB
				end if

				if ( OGenPotE .eqv. .TRUE. ) then
					if ( OIESMode .eqv. .FALSE. ) then
						write(uEner,"(3X,2F8.2,3F18.10)") angX, angY, rhoVal
					else 
						write(uEner,"(3X,2F8.2,F18.10)") angX, angY, rhoVal
					end if
				end if
				angZ = angZ + incScanZ
			end do
			angY = angY + incScanY
			angZ = minScanZ
		  end do
		if ( OGenPotE .eqv. .TRUE. ) then
			write(uEner,*) 
		end if  
		

		angX = angX + incScanX
		angY = minScanY
	end do
	
	if ( OIESMode .eqv. .FALSE. ) then
		call write_rline(63)
	else
		call write_rline(51)
		incScanZ = incScanZOld
	end if
	
	write(*,'(/5X,A)') " > energetic minimum:"
	if ( OIESMode .eqv. .FALSE. ) then
		write(*,'(/A35,F16.8/)') " KD1 |+/-15/2> contribution = ", maxValue
		write(*,'(A35,F10.2,A)') " alpha = ", maxRotX/pi*180d0, " degree"
		write(*,'(A35,F10.2,A)') "  beta = ", maxRotY/pi*180d0, " degree"
	else
		write(*,'(/A35,F16.8,A/)') " relative energy = ", minValue, " a.u."
		write(*,'(A35,F10.2,A)') " alpha = ", minRotX/pi*180d0, " degree"
		write(*,'(A35,F10.2,A)') "  beta = ", minRotY/pi*180d0, " degree"
	end if
	
	write(*,*)

	if ( OGenPotE .eqv. .TRUE. ) then
		close(uEner)
	end if  

end subroutine

subroutine do_scan_ab
	use data_c
	use data_mo
	use data_grid
	implicit none
!
	double precision :: angX, angY, angZ
	double precision :: rhoVal, incScanZOld
!
	angX = minScanX
	angY = minScanY
	angZ = minScanZ
	
	write(*,'(/X,A/)') "PERFORMING SCAN OF EULER ANGLES ALPHA AND BETA ..."
	if ( OIESMode .eqv. .TRUE. ) then
		write(*,'(A33,A/)') "energies determined by: ", "electrostatic calculation"
	else
		write(*,'(A33,A/)') "energies determined by: ", "crystal-field theory"
	end if
	
	write(*,'(A33,F7.2,A,F7.2,A,F6.2,A)') "alpha scan range: ", minScanX, " deg -- ", maxScanX, " deg (inc. ", incScanX, " deg)"
	write(*,'(A33,F7.2,A,F7.2,A,F6.2,A)') " beta scan range: ", minScanY, " deg -- ", maxScanY, " deg (inc. ", incScanY, " deg)"
	write(*,*)

	if ( OGenPotE .eqv. .TRUE. ) then
		open(unit=uEner, file=OLFile, action='write',status='replace')
	end if
	
	call write_rline(63)
	write(*,"(18X,3A12,A25)") "alpha", "beta", "gamma", "relative energy"
	write(*,"(18X,3A12,A25)") "/ degree", "/ degree", "/ degree", "/ cm-1"
	call write_rline(63)

	do while ( angX .lt. maxScanX )
		fRotA = angX * pi / 180d0
		if ( angX .ne. minScanX ) then
			write(*,*)
		end if
		do while ( angY .le. maxScanY )
			fRotB = angY * pi / 180d0
				
				if ( OIESMode .eqv. .FALSE. ) then
					call calc_cf_energy(fRotA, fRotB, fRotG, rhoVal)
					write(*,"(18X,3F12.2,F25.8)") angX, angY, angZ, rhoVal
				else
					call cpot_calc_energy(fRotA, fRotB, fRotG, rhoVal)
					write(*,"(30X,2F12.2,F25.8)") angX, angY, rhoVal
				end if
				
				if ( OTraj .eqv. .TRUE. ) then
					call write_xyz_file(fRotA, fRotB, .TRUE., rhoVal) 
				end if 
				
				if ( rhoVal .lt. minValue ) then
					minValue = rhoVal
					minRotX = fRotA
					minRotY = fRotB
				end if
				if ( rhoVal .gt. maxValue ) then
					maxValue = rhoVal
					maxRotX = fRotA
					maxRotY = fRotB
				end if

				if ( OGenPotE .eqv. .TRUE. ) then
					if ( OIESMode .eqv. .TRUE. ) then
						write(uEner,"(3X,2F8.2,F18.10)") angX, angY, rhoVal
					else
						write(uEner,"(3X,3F8.2,F18.10)") angX, angY, angZ, rhoVal
					end if
				end if	

			angY = angY + incScanY
			angZ = minScanZ
		end do
		if ( OGenPotE .eqv. .TRUE. ) then
			write(uEner,*) 
		end if  

		angX = angX + incScanX
		angY = minScanY
	end do
	
	if ( OIESMode .eqv. .FALSE. ) then
		call write_rline(63)
	else
		call write_rline(51)
		incScanZ = incScanZOld
	end if
	
	write(*,'(/5X,A)') " > energetic minimum:"
	write(*,'(/A35,F16.8,A/)') " relative energy = ", minValue, " 1/cm"
	write(*,'(A35,F10.2,A)') " alpha = ", minRotX/pi*180d0, " degree"
	write(*,'(A35,F10.2,A)') "  beta = ", minRotY/pi*180d0, " degree"
	write(*,*)

	if ( OGenPotE .eqv. .TRUE. ) then
		close(uEner)
	end if  

end subroutine


subroutine do_opt(bequiet)
	use data_c
	use data_grid
	implicit none
!
	logical, intent(in) :: bequiet
	double precision :: GDyX, GDyY, GDyZ    ! current grid coordinates	
	double precision :: RDyX, RDyY, RDyZ    ! rotated coordinates   	
	logical :: OConverged = .FALSE.
	double precision :: rhoMin(9), SGridX(9), SGridY(9), calc_vec_angle
	logical :: found
	integer :: i, j, minI
!
	j = 1
	do i = 1, 9
		if ( i .ne. 5 ) then
			SGridX(i) = dcos( 0.25d0*pi*dble(j-1) )
			SGridY(i) = dsin( 0.25d0*pi*dble(j-1) )
			j = j + 1
		else
			SGridX(i) = 0d0
			SGridY(i) = 0d0
		end if
	end do
	
	if ( bequiet .eqv. .FALSE. ) then
	
		write(*,*)
		if ( OMax .eqv. .FALSE. ) then
			write(*,*) "PERFORMING A MINIMIZATION OF ELECTROSTATIC ENERGY ..."
		else
			write(*,*) "PERFORMING A MAXIMIZATION OF KD1 |+/-15/2> CONTRIBUTION ..."
		end if
		write(*,*)
		write(*,'(A35)') "Convergence criteria: "
		write(*,'(A35,F18.5,A)') "dE: ", EnergyConv, " cm-1"
		write(*,'(A35,F18.5,A)') "dRot: ", AngConv, " deg"
		if ( OIESMode .eqv. .FALSE. ) then
			write(*,'(A35,F18.2,A)') "gamma: ", fRotG*180d0/pi, " deg"
		end if
		write(*,'(A35,I14)'    ) "max. iterations: ", MaxIter
		
		write(*,*)

		if ( OMagVec .eqv. .TRUE. ) then
			SepLine = "-----------------------------------------------------------------------------"
			write(*,'(X,A)') trim(SepLine)
			write(*,"(66X,A12)") "dev. from"
			write(*,"(2A3,2A8,2A18,A8,A12)") "#", " SI", "alpha", "beta", "rel. energy", "D[dE]/D[rot]", "rot.", "ab initio"
			write(*,"(2A3,2A8,2A18,A8,A12)") "", "", "deg", "deg", "cm-1", "cm-1 / deg", "deg", "deg"

			GDyX = 0.0d0
			GDyY = 0.0d0
			GDyZ = 1.0d0
		else
			SepLine = "--------------------------------------------------------------------------"			
			write(*,'(A)') trim(SepLine)
			write(*,"(2X,2A3,2A8,2A15,A12)") "#", " SI", "alpha", "beta", "E", "D[dE]/D[rot]", "angle inc."
			write(*,"(2X,2A3,2A8,2A15,A12)") "", "", "deg", "deg", "cm-1", "cm-1 / deg", "deg"
		end if
		write(*,'(X,A)') trim(SepLine)
		write(*,*)

	end if
	
	nStep = 0
	minRotX = fRotA
	minRotY = fRotB

	if ( OIESMode .eqv. .FALSE. ) then
		call calc_cf_energy(fRotA, fRotB, fRotG, rho)
	else
		call cpot_calc_energy(fRotA, fRotB, 0d0, rho)
	end if

	rhoMin(5) = rho
	found = .TRUE.

	do while ( nStep .lt. MaxIter )

		minValue = rho

		if ( found .eqv. .TRUE. ) then
			nStep = nStep + 1
			if ( bequiet .eqv. .FALSE. ) then
				write(*,'(2I3,2F8.2,2F18.8,F8.2)',advance='no') & 
					nStep, 0, fRotA*180d0/pi, fRotB*180d0/pi, rho, ( rho - rhoMin(5) ) / incStep, incStep*180d0/pi
			end if
			if ( OTraj .eqv. .TRUE. ) then

				call write_xyz_file(fRotA, fRotB, .TRUE., rho) 

			end if 
			
			if ( rho .gt. maxValue ) then
				maxValue = rho
				maxRotX = fRotA		  
				maxRotY = fRotB
			end if

			if ( bequiet .eqv. .FALSE. ) then
				if ( OMagVec .eqv. .TRUE. ) then       
					call rotate_euler_3d(fRotA, fRotB, 0d0, 0d0, 0d0, 1d0, RDyX, RDyY, RDyZ)
					MagVecDev = calc_vec_angle(RDyX, RDyY, RDyZ, MagX, MagY, MagZ)
					write(*,"(F12.2)") MagVecDev
				else
					write(*,*)
				end if
			end if
		end if

		found = .FALSE.

		! convergence check
		if ( ( abs(rho - rhoMin(5)) .lt. EnergyConv ) .AND. &
		      ( abs(rho - rhoMin(5)) .ne. 0.0 ) .AND. &
		      ( incStep*180d0/pi .lt. AngConv ) ) then
			minRotX = fRotA
			minRotY = fRotB
			OConverged = .TRUE.
			exit
		end if

		rhoMin(5) = rho

		do i = 1, 9

			fRotA = minRotX + SGridX(i) * incStep
			fRotB = minRotY + SGridY(i) * incStep

			if ( i .ne. 5 ) then
				if ( OIESMode .eqv. .FALSE. ) then
					call calc_cf_energy(fRotA, fRotB, fRotG, rhoMin(i))
				else
					call cpot_calc_energy(fRotA, fRotB, 0d0, rhoMin(i))
				end if
			end if
			
			if ( bequiet .eqv. .FALSE. ) then
				write(*,'(3X,I3,2F8.2,2F18.8)') & 
					i, fRotA*180d0/pi, fRotB*180d0/pi, rhoMin(i), (rhoMin(i) - rhoMin(5)) / incStep
			end if
			
			if ( OMax .eqv. .FALSE. ) then

				if ( minValue .gt. rhoMin(i) ) then
					minValue = rhoMin(i)
					found = .TRUE.
					minI = i
				end if

			else

				if ( minValue .lt. rhoMin(i) ) then
					minValue = rhoMin(i)
					found = .TRUE.
					minI = i
				end if
			end if

		end do
		
		if ( bequiet .eqv. .FALSE. ) then
			write(*,*)
		end if

		if ( found .eqv. .TRUE. ) then
			minRotX = minRotX + SGridX(minI) * incStep
			minRotY = minRotY + SGridY(minI) * incStep
			
			fRotA = minRotX
			fRotB = minRotY
			rho = minValue
			
			if ( minI .eq. 5 ) then
				incStep = incStep / 2d0
			end if
		else
			fRotA = minRotX
			fRotB = minRotY
			incStep = incStep / 2d0
			
			if ( incStep*180d0/pi .lt. AngConv ) then
				minRotX = fRotA
				minRotY = fRotB
				OConverged = .TRUE.
				exit
			end if
			
		end if
	
	end do

	if ( bequiet .eqv. .FALSE. ) then
		write(*,*)
		write(*,'(2X,A)') trim(SepLine)
	end if

	if ( OConverged .eqv. .FALSE. ) then
		write(*,"(25X,A)") "! OPTIMIZATION FAILED TO CONVERGE !"
		write(*,*)
		stop
	else
		if ( bequiet .eqv. .FALSE. ) then
			write(*,'(/26X,A,I4,A)') "CONVERGENCE AFTER ", nStep, " ITERATIONS"
			write(*,*)
		end if
	end if

end subroutine

subroutine check_angles
	use data_c
	use global_c
	implicit none
!
	if ( fRotA .lt. 0d0 ) then
		fRotA = fRotA + TwoPi
	end if
	if ( fRotB .lt. -TwoPi ) then
		fRotB = fRotB + TwoPi
	end if
	if ( fRotG .lt. 0d0 ) then
		fRotG = fRotG + TwoPi
	end if
	
end subroutine

subroutine print_results
	use data_c
	use data_mo
	use data_grid
	implicit none
!
	double precision :: RDyX, RDyY, RDyZ, calc_vec_angle
!
	call rotate_euler_3d(minRotX, minRotY, 0d0, 0d0, 0d0, 1d0, RDyX, RDyY, RDyZ)
	
	write(*,'(/X,A/)') "=============================== FINAL RESULTS ==============================="
	write(*,'(A35)') "results of the "
	if ( OIESMode .eqv. .FALSE. ) then
		
		write(*,'(A35)') "LFT approach calculation: "
		write(*,*)
		write(*,'(A35,F16.8)') "KD1 |+/-15/2> contribution: ", minValue
	else
		write(*,'(A35)') "improved electrostatic calculation: "
		write(*,*)
		write(*,'(A35,F16.8,A)') "minimal energy: ", minValue, " a.u."	    
	end if
	write(*,'(A35,F12.2,A)') "alpha: ", minRotX*180.0/pi, " deg"
	write(*,'(A35,F12.2,A)') "beta: ", minRotY*180.0/pi, " deg"
	write(*,'(A35,3F14.8)')  "quantization axis: ", RDyX, RDyY, RDyZ

	if ( OMagVec .eqv. .TRUE. ) then

		MagVecDev = calc_vec_angle(RDyX, RDyY, RDyZ, MagX, MagY, MagZ)
		write(*,'(A35,F12.2,A)') "deviation from reference: ", MagVecDev, " deg"
		
	end if
	
	write(*,'(/X,A)') "=============================== FINAL RESULTS ==============================="
	
end subroutine

subroutine do_scan_and_opt
	use data_c
	implicit none
!
	call do_scan()

	if ( OMax .eqv. .FALSE. ) then
		fRotA = minRotX
		fRotB = minRotY
		fRotG = 0d0
	else
		fRotA = maxRotX
		fRotB = maxRotY
		fRotG = 0d0
	end if
	
	if ( OMixedMode .eqv. .TRUE. ) then
		OIESMode = .FALSE.
	end if
	
	if (OIESMode .eqv. .TRUE. ) then
		OMax = .FALSE.
	else
		OMax = .TRUE.
	end if

	incStep    =  0.5d0 * dsqrt(incScanX**2 + incScanY**2 ) * pi/180d0
	call do_opt(.FALSE.)

end subroutine

subroutine do_lin_scan
	use data_c
	implicit none
!
	integer :: i
	double precision :: minEner, maxEner, cx, cy, ce
!
	write(*,*) "PERFORMING A LINEAR SCAN ..."
	write(*,*)
	
	if ( ( OIESMode .eqv. .FALSE. ) .and. ( OMixedMode .eqv. .FALSE. ) ) then
		write(*,*) "ERROR! LINEAR SCAN ONLY AVAILABLE FOR ELECTROSTATIC CALCULATIONS!"
		stop
	end if
	
	minEner = 1d8
	maxEner = 0d0
	
	write(*,'(5X,A)') "> FIRST COORDINATE:"
	write(*,'(10X,A35,F8.2,A)') "alpha = ", LinScanX1, " deg."
	write(*,'(10X,A35,F8.2,A)') " beta = ", LinScanX2, " deg."
	write(*,*)
	write(*,'(5X,A)') "> SECOND COORDINATE:"
	write(*,'(10X,A35,F8.2,A)') "alpha = ", LinScanX2, " deg."
	write(*,'(10X,A35,F8.2,A)') " beta = ", LinScanY2, " deg."
	write(*,*)
	write(*,'(10X,A35,I8)') "number of steps: ", LinScanS
	write(*,'(10X,A35,F8.2,A)') "increment alpha = ", LinScanDX, " deg."
	write(*,'(10X,A35,F8.2,A)') " increment beta = ", LinScanDY, " deg."
	write(*,*)
	write(*,*)
	call write_rline(70)
	
	do i = 1, LinScanS
	
		cx = ( LinScanX1 + dble(i-1) * LinScanDX ) / 180d0 * pi
		cy = ( LinScanY1 + dble(i-1) * LinScanDY ) / 180d0 * pi
		
		call cpot_calc_energy(cx, cy, 0d0, ce)
		
		if ( ce .lt. minEner ) then
			minEner = ce
		end if
		if ( ce .gt. maxEner ) then
			maxEner = ce
		end if
		
		write(*,'(16X,I5,2F10.2,F20.8)') i, cx*180d0/pi, cy*180d0/pi, ce
		
	end do
	
	call write_rline(70)
	write(*,*)
	write(*,'(10X,A35,F16.8,A)') "maximum: ", maxEner, " a.u."
	write(*,'(10X,A35,F16.8,A)') "minimum: ", minEner, " a.u."
	write(*,'(10X,A35,F16.8,A)') "difference: ", (maxEner-minEner), " a.u."
	write(*,'(10X,A35,F16.8,A)') "difference: ", (maxEner-minEner)*au2rcm, " 1/cm"
	write(*,*)
	write(*,*)
	
	stop
	
end subroutine

subroutine check_metal_center
	use data_mo
	implicit none
!
	integer :: i
	double precision :: distm, dist
!
	distm = 1d3
	
	do i = 1, NA
		dist = dsqrt((Atoms(i)%x-EDyX)**2 + (Atoms(i)%y-EDyY)**2 + (Atoms(i)%z-EDyZ)**2)
		if ( ( dist .lt. distm ) .AND. ( Atoms(i)%ghost .eqv. .FALSE. ) ) then
			distm = dist
		end if
	end do
	
	if ( distm .gt. 3d0 / au2pm ) then
		write(*,*) "*****************************************************************************"
		write(*,*) "*                                   ERROR!                                  *"
		write(*,*) "*                                                                           *"
		write(*,*) "*             Please check the position of the lanthanide ion!              *"
		write(*,*) "*****************************************************************************"
		write(*,*)
		if ( OExpert .eqv. .FALSE. ) then
			stop
		end if
	end if
	
end subroutine

! *************************************************************************
! * program starts here
! *************************************************************************
program rhover
	use data_c
	use data_mo
	use data_grid
	implicit none
!
	integer :: load_DEPP
!       
	! prints program header and sets default values
	call print_header	
	call set_default_values
	call generate_coefficients
	
	! checks for input file
	if ( iargc() .eq. 0 ) then
		write(*,*) "ERROR! No input file given!"
		write(*,*)
		stop
	else
		call getarg(1, InpFile)
	end if

	! reads the input file
	call parse_input_file(InpFile)
	call print_jobinfo
	call print_stevens_factors(OVerbose)
	
	if ( OPCM .eqv. .FALSE. ) then
		! reads data from ab initio calculation
		write(*,*) "> PREPARING JOB ..."
		call read_molden_file(MOFile)
		call reorder_cgtos
		call print_bs_info
		
		if ( OPrintBasis .eqv. .TRUE. ) then
			call print_cgto_info
			call print_pgto_info
		end if
		
		! calculates properties
		call check_metal_center
		call allocate_matrices
		call calc_overlap_matrix(.FALSE.)
		call calc_mulliken_charges(.TRUE.)
		
		! further checks
		call check_ghost_atoms
		if ( ( ODelete4f .eqv. .TRUE. ) .or. ( ODeleteAll .eqv. .TRUE. ) ) then 
			call calc_overlap_matrix(.TRUE.)
		end if
		call check_for_ecp
		
	else
	
		write(*,*) "> rhOver IS IN PCM MODE! ..."
		write(*,*)
		
	end if	
	
	! generates 4f electron grid
	call lebedev_init
	call spherical_dy_grid_init
	
	if ( OPCM .eqv. .FALSE. ) then
		
		! tries to load diamagnetic exchange pseudo-potential
		if ( load_DEPP(DEPPFile) .ne. 2 ) then
			! calculates the electrostatic potential
			call calc_potentials_nuc
			call calc_potentials_anal_elec
		
			if ( ONoDEPPFile .eqv. .FALSE. ) then
				call save_DEPP(DEPPFile)
			end if
		end if
		
		if ( OExportGPs .eqv. .TRUE. ) then
			call save_GPs(GPFile)
		end if
		
		! calculates ligand density and Slater-X, if necessary
		if ( OLDAX .eqv. .TRUE. ) then	
			call calc_ligand_dens
			call calc_ligand_ldax
		end if
		
	else
		
		call calc_potentials_pcm
		
	end if
	
	if ( JobType .eq. 1 ) then
		
		write(*,'(/X,A/)') "> CALCULATING SINGLE EULER ROTATION ..."
		call write_rline(75)
		if ( OIESMode .eqv. .FALSE. ) then
			write(*,"(10X,3A15,A24)") "alpha / deg", "beta / deg", "gamma / deg", "|+/-15/2> contr. / 1"
		else
			write(*,"(10X,3A15,A24)") "alpha / deg", "beta / deg", "gamma / deg", "rel. energy / cm-1"
		end if
		call write_rline(75)
		write(*,*)
		
		call check_angles
		
		fRotG = 0d0
		if ( OIESMode .eqv. .FALSE. ) then
			call calc_cf_energy(fRotA, fRotB, fRotG, minValue)
		else
			call cpot_calc_energy(fRotA, fRotB, fRotG, minValue)
		end if
		
		minRotX = fRotA
		minRotY = fRotB
		
		write(*,"(10X,3F15.2,F24.8)") fRotA*180d0/pi, fRotB*180d0/pi, fRotG*180d0/pi, minValue
		write(*,*)
		call write_rline(75)
		write(*,*)
		
	else if ( JobType .eq. 2 ) then
		
		call do_scan()
		
	else if ( JobType .eq. 3 ) then
		
		call do_opt(.FALSE.)
		
	else if ( JobType .eq. 4 ) then
		
		call do_scan_and_opt()
		
	else if ( JobType .eq. 5 ) then
		
		call do_lin_scan()
		
	else 
		write(*,*) "ERROR! Unknown job type!"
	end if
	
	! prints the final results
	call print_results
	
	if ( ( OGenXYZ .eqv. .TRUE. ) .and. ( OPCM .eqv. .FALSE. ) ) then
		call write_xyz_file(minRotX, minRotY, .FALSE., 0d0)
	end if

	if ( ( OIESMode .eqv. .TRUE. ) .and. ( OMixedMode .eqv. .FALSE. ) ) then
		write(*,*) "> LIGAND-FIELD PART IS SKIPPED!"
		write(*,*)
	else
		call calc_cf_energies_final(minRotX, minRotY, 0.0d0)
	end if
	
	write(*,'(/X,A/)') "END."
	
end program

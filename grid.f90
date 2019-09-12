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

module data_grid

	use global_c

	integer :: NLGP, NGridPoints
	
	type t_gridpoint
		integer*4 :: id
		integer*4 :: aid
		double precision :: x
		double precision :: y
		double precision :: z
		double precision :: dV
		double precision :: RadWFdV
		double precision :: RadWF
		double precision :: dist
		double precision :: NucPotVal
		double precision :: ElecPotVal
		double precision :: LDens
	end type t_gridpoint
	
	type t_pointcharge
		integer*4 :: id
		double precision :: x
		double precision :: y
		double precision :: z
		double precision :: dist
		double precision :: Charge
	end type t_pointcharge
	
	type t_smallgp
		double precision :: x
		double precision :: y
		double precision :: z
		double precision :: ElecPotVal
	end type t_smallgp
	
	type(t_gridpoint), dimension(:), allocatable :: GridPoints
	type(t_pointcharge), dimension(:), allocatable :: PointCharges
	
end module

subroutine spherical_dy_grid_init
	use global_c
	use data_grid
	use data_mo
	use data_lebedev
	implicit none
!
	double precision :: RadWF, xx, yy, zz, totaldens
	double precision :: yrot, xrot, zrot, xi, ri, wi, pxi, shelldens, dVSpher
	integer :: i, k, l, Ltype
	integer :: NRealRadShells, get_radial_shells
!
	NLigGP = 0
	pxi = 0.90d0
	
	allocate(GridPoints(NGP))
	
	write(*,*) "> CONSTRUCTING NUMERICAL 4f ELECTRON GRID ..."
	write(*,*)
#ifdef DEBUG
	write(*,'(A35,I14)') "max. grid points: ", NGP
#endif
	write(*,'(A35,I14)') "grid size: ", nGrid
	write(*,*)
	
	NLigRadShells = NRS(nGrid)
	NRealRadShells = get_radial_shells(NLigRadShells, pxi)
	
	write(*,'(A35,I14)') "number of radial shells: ", NLigRadShells
#ifdef DEBUG
	write(*,'(A35,I14)') "number of real radial shells: ", NRealRadShells 
#endif
	write(*,'(A35,F14.8,A)') "inner radial cut-off: ", Ln4fCutOff_I, " a0"
	write(*,'(A35,F14.8,A)') "outer radial cut-off: ", Ln4fCutOff_O, " a0"
	write(*,'(A35,F14.8,A)') "xi value: ", pxi, " a0"
	write(*,*)
	
	write(*,'(A35)') "spherical grid: "
	write(*,'(A35)',advance='no') "pruning: "
	call write_option_bool(OPruning)
	
	write(*,'(A35)',advance='no') "random rotation: "
	call write_option_bool(ORandRot)
	write(*,*)
	
	if ( OVerbose .eqv. .TRUE. ) then
		write(*,'(5X,4A)') ('------------------', i=1, 4)
		write(*,'(3X,2A11,A26,A14,A12)') "radius", "radius", "total 4f", "density",  "4f grid"
		write(*,'(3X,2A11,A26,A14,A12)') "/ a0", "/ angstrom", "charge density", "per shell",  "points"
		write(*,'(5X,4A)') ('------------------', i=1, 4)
	end if
	
	i = 0
	totaldens = 0d0
	
	do k = NRealRadShells, 1, -1
	
		! M4 radial scheme by Treutler et al. with Chebyshev quadrature of second kind (T2) 
		xi = dcos(real(k) * Pi/(NRealRadShells+1.0d0))
		ri = pxi / dlog(2d0) * (1d0+xi)**Ahlrichs_Alpha * log(2d0/(1d0-xi))

		wi = Pi/(NRealRadShells+1.0d0) * pxi**3 * ( ( xi+1d0 )**Ahlrichs_Alpha / dlog(2d0) )**3 * ( &
			 dsqrt( ( xi+1d0 ) / ( 1d0-xi ) ) * (dlog( ( 1d0-xi ) / 2d0 ))**2 - &
			 Ahlrichs_Alpha * dsqrt( ( 1d0 - xi ) / ( 1d0+xi ) ) * (log( ( 1d0-xi ) / ( 2d0 ) ))**3 )
		
		! check for radial cut-off
		if ( ( ri .gt. Ln4fCutOff_I ) .AND. ( ri .lt. Ln4fCutOff_O ) ) then
		
			shelldens = 0.0d0
			
			call Radial_4f_WF(ri, RadWF)
			RadWF = RadWF**2
		 
			if ( ORandRot .eqv. .TRUE. ) then
				call random_number(xrot)
				call random_number(yrot)
				call random_number(zrot)

				xrot = xrot*2d0*pi - pi
				yrot = yrot*1d0*pi - pi/2d0
				zrot = zrot*2d0*pi - pi
			end if
			
			call get_lebedev_grid_type(k, NRealRadShells, Ltype)
			call get_lebedev_grid_size(Ltype, NRad)
			
			do l = 1, NRad
					  
				i = i + 1
				NLigGP = NLigGP + 1
				  
				if ( i .ge. NGP ) then
					write(*,*) "ERROR! Please increase the number of grid points!"
					write(*,*)
					stop
				end if
				  
				GridPoints(i)%id = i
				GridPoints(i)%dist = ri

				call get_lebedev_grid_point(Ltype, l, GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z, dVSpher)
				GridPoints(i)%x = GridPoints(i)%x * ri
				GridPoints(i)%y = GridPoints(i)%y * ri
				GridPoints(i)%z = GridPoints(i)%z * ri
				GridPoints(i)%dV = dVSpher * wi
				GridPoints(i)%RadWF = RadWF 
				GridPoints(i)%RadWFdV = GridPoints(i)%dV * GridPoints(i)%RadWF

				if ( ORandRot .eqv. .TRUE. ) then
					call rotate_euler_3d(xrot, yrot, zrot, GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z, xx, yy, zz)
					GridPoints(i)%x = xx
					GridPoints(i)%y = yy
					GridPoints(i)%z = zz
				end if

				shelldens = shelldens + NumberOf4fElectrons * GridPoints(i)%RadWFdV 
				
			  end do
			  
			  totaldens = totaldens + shelldens

			  if ( OVerbose .eqv. .TRUE. ) then
				write(*,'(3X,2F11.4,F16.6,F9.2,A,F14.6,I12)') ri, ri*au2pm, totaldens, totaldens/NumberOf4fElectrons*100d0,"%", shelldens, i
!				 write(*,'(3X,3F22.16)') 0d0, 0d0, ri
			  endif
			  
		end if
		
	end do

	NGP = NLigGP
	NLGP = NGP
	
	if ( OVerbose .eqv. .TRUE. ) then
		write(*,'(5X,4A)') ('------------------', i=1, 4)
		write(*,*)
		write(*,*)
	end if
	
end subroutine


function get_radial_shells(n, pxi) result(res)
	use global_c
	implicit none
!
	integer, intent(in) :: n
	double precision, intent(in) :: pxi
	integer :: res, k, nshells, nrs
	double precision :: ri, xi
	logical :: found
!

	found = .FALSE.
	nrs = n

	do while ( found .eqv. .FALSE. )
	
		nshells = 0
		do k = nrs, 1, -1
		
			xi = dcos(real(k) * Pi/(nrs+1.0d0))
			ri = pxi / dlog(2d0) * (1d0+xi)**Ahlrichs_Alpha * log(2d0/(1d0-xi))
			
			if ( ( ri .gt. Ln4fCutOff_I ) .AND. ( ri .lt. Ln4fCutOff_O ) ) then
				nshells = nshells + 1
			end if
		end do
		
		if ( nshells .ge. n) then
			res = nrs
			found = .TRUE.
		end if
		nrs = nrs + 1
	end do
		
end function

subroutine save_DEPP(filename)
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	character (len=sMaxBuffer), intent(in) :: filename
	integer :: i
!	
	write(*,*) "   SAVING DEPP ..."
	write(*,*)
	
	open(unit=uPotF, file=filename,action='write', status='replace',form='formatted')
	
	write(uPotF,'(A)') "RHOVER GRID DATA"
	write(uPotF,'(2I6,L6,I6,2F6.3)') nGrid, NLGP, ODelete4f, NLigRadShells, Ln4fCutOff_I, Ln4fCutOff_O
	write(uPotF,'(6X,3E23.15)') EDyX, EDyY, EDyZ
	
	do i = 1, NLGP
		write(uPotF, '(I6,6E23.15)') GridPoints(i)%id, GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z, GridPoints(i)%dV, GridPoints(i)%NucPotVal, GridPoints(i)%ElecPotVal
	end do
	
	close(uPotF)
	
	write(*,'(A35,A)') "data written to file: ", filename
	write(*,*)
	call write_done

end subroutine

subroutine save_GPs(filename)
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	character (len=sMaxBuffer), intent(in) :: filename
	integer :: i
!	
	write(*,*) "   SAVING DEPP GRID POINTS ..."
	write(*,*)
	
	open(unit=uPotF, file=filename,action='write', status='replace',form='formatted')
	
	write(uPotF, *) NLGP
	
	do i = 1, NLGP
		write(uPotF, '(3E23.15)') GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z
	end do
	
	close(uPotF)
	
	write(*,'(A35,A)') "data written to file: ", filename
	write(*,*)
	call write_done

end subroutine

function load_DEPP(filename) result(res)
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	character (len=sMaxBuffer), intent(in) :: filename
	integer :: res
	integer :: i, FNLGP, iost, fnGrid, fNLigRadShells
	logical :: fexist, fDelete4f
	double precision :: fLn4fCutOff_I, fLn4fCutOff_O, fEDyX, fEDyY, fEDyZ
	character (len=sMaxBuffer) :: fStr
!	
	res = 0
	inquire(file=filename,exist=fexist)
	
	if ( fexist .eqv. .TRUE. ) then
		res = 1
		write(*,*) "DEPP FILE FOUND, TRYING TO LOAD ..."
		write(*,*)
		open(unit=uPotF, file=filename, action='read', status='old', iostat=iost)

		read(uPotF, '(A)') fStr
		
		if ( trim(fStr) .eq. "RHOVER GRID DATA" ) then
			read(uPotF, '(2I6,L6,I6,2F6.3)') fnGrid, FNLGP, fDelete4f, fNLigRadShells, fLn4fCutOff_I, fLn4fCutOff_O
			read(uPotF, '(6X,3E23.15)') fEDyX, fEDyY, fEDyZ
			
			! compare grid data settings
			if ( ( nGrid .eq. fnGrid ) .and. ( FNLGP .eq. NLGP ) .and. ( fDelete4f .eqv. ODelete4f ) .and. ( abs(fLn4fCutOff_I-Ln4fCutOff_I) .lt. 1d-4 ) &
				 .and. ( abs(fLn4fCutOff_O-Ln4fCutOff_O) .lt. 1d-4 ) .and. (fNLigRadShells .eq. NLigRadShells ) &
				 .and. ( abs(fEDyX-EDyX) .lt. 1d-4 ) .and. ( abs(fEDyY-EDyY) .lt. 1d-4 ) .and. ( abs(fEDyZ-EDyZ) .lt. 1d-4 ) ) then
			
				do i = 1, FNLGP
					read(uPotF, '(I6,6E23.15)') GridPoints(i)%id, GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z, GridPoints(i)%dV, GridPoints(i)%NucPotVal, GridPoints(i)%ElecPotVal
				
					GridPoints(i)%dist = dsqrt(GridPoints(i)%x**2 + GridPoints(i)%y**2 + GridPoints(i)%z**2)
					
				end do
				res = 2
				call write_done
			else
				write(*,'(65X,A)') " ... FAILED!"
				write(*,*)
			end if
		else
			write(*,'(65X,A)') " ... FAILED!"
			write(*,*)
		end if
		close (uPotF)
	end if

end function

subroutine calc_potentials_nuc
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	integer :: i, j
	double precision :: dist
!
	write(*,*) "CALCULATING ANALYTICAL NUCLEAR POTENTIAL ..."

	do i = 1, NLGP
		GridPoints(i)%NucPotVal = 0d0
		do j = 1, NA
			if ( ( Atoms(j)%ghost .eqv. .FALSE. ) .AND. ( Atoms(j)%skip .eqv. .FALSE. ) ) then
				dist = dsqrt( (GridPoints(i)%x-Atoms(j)%x+EDyX )**2 + (GridPoints(i)%y-Atoms(j)%y+EDyY)**2 + (GridPoints(i)%z-Atoms(j)%z+EDyZ )**2 )
				if ( dist .gt. 1.0d-5 ) then
					GridPoints(i)%NucPotVal = GridPoints(i)%NucPotVal + dble(Atoms(j)%charge)/dist
				end if
			end if
		end do
	end do
	
	call write_done
	
end subroutine

subroutine calc_potentials_pcm
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	integer :: i, j
	double precision :: dist
!
	write(*,*) "CALCULATING PCM-based POTENTIAL ..."

	do i = 1, NLGP
		GridPoints(i)%NucPotVal = 0d0
		GridPoints(i)%ElecPotVal = 0d0
		do j = 1, NCustomCharges
			dist = dsqrt( (GridPoints(i)%x-PointCharges(j)%x+EDyX )**2 + (GridPoints(i)%y-PointCharges(j)%y+EDyY)**2 + (GridPoints(i)%z-PointCharges(j)%z+EDyZ )**2 )
			if ( dist .gt. 1.0d-6 ) then
				GridPoints(i)%ElecPotVal = GridPoints(i)%ElecPotVal - ( dble(PointCharges(j)%Charge) / dist )
			end if
		end do
	end do
	
	call write_done
	
end subroutine

subroutine calc_potentials_anal_elec
	use data_mo
	use global_c
	use data_grid
	implicit none
!
	double precision :: sssval, ssssval, gamma, cpr2
	double precision, dimension(3) :: PAB, PA, PB, C
	double precision, dimension(:,:), allocatable :: ExpArr
	integer :: k, dxyz, vi, vj, vk, pgto_get_l
	double precision, dimension(9, 3) :: GVal
	double precision :: G_binom_factors, boys_func
	integer*1, dimension(:,:,:), allocatable :: CValMaxDim
	integer*4, dimension(:,:), allocatable :: L1Arr
	type(t_smallgp), dimension(:), allocatable :: TempElecPot
	integer*4 :: i, j, ic, jc, ii
	integer :: t1, t2, cmax, progress, cprogress
#ifdef __GFORTRAN__
	integer :: crate
#else
	double precision :: crate
#endif
!
	! delete old values
	do i = 1, NLGP
		GridPoints(i)%ElecPotVal = 0d0
	end do
	
	progress = NumCGTO / 10
	cprogress = 0
	
	write(*,*) "CALCULATING ANALYTICAL ELECTRONIC POTENTIAL ..."
	write(*,*)
	
	! calculates constant stuff before
	allocate(ExpArr(NumPGTO, NumPGTO))
	allocate(CValMaxDim(3, NumPGTO, NumPGTO))
	allocate(L1Arr(3, NumPGTO))
	
	!$OMP PARALLEL PRIVATE(TempElecPot)
	
	!$OMP DO PRIVATE(dxyz)
	do i = 1, NumPGTO
		do dxyz = 1, 3
			L1Arr(dxyz, i) = pgto_get_l(PGTOs(i)%shelltype, PGTOs(i)%subtype, dxyz)
		end do
	end do
	!$OMP END DO
	
	!$OMP BARRIER
	
	!$OMP DO PRIVATE(j, dxyz)
	do i = 1, NumPGTO
		do j = 1, i
			do dxyz = 1, 3
				CValMaxDim(dxyz, i, j) = L1Arr(dxyz, i) + L1Arr(dxyz, j)
				if ( i .ne. j ) then
					CValMaxDim(dxyz, j, i) = CValMaxDim(dxyz, i, j)
				end if
			end do
		end do
	end do
	!$OMP END DO
	
	!$OMP DO PRIVATE(j, gamma)
	do i = 1, NumPGTO
		do j = 1, i
			gamma = PGTOs(i)%coeff_alpha + PGTOs(j)%coeff_alpha
			ExpArr(i, j) = D_CGTO(PGTOs(i)%cgtoid, PGTOs(j)%cgtoid) * FourPi * PGTOs(i)%coeff_d * PGTOs(j)%coeff_d / gamma * dexp(-PGTOs(i)%coeff_alpha*PGTOs(j)%coeff_alpha * AtomsDistMatSq(PGTOs(i)%atomid, PGTOs(j)%atomid) / gamma )
			if ( i .ne. j ) then
				ExpArr(i, j) = 2d0 * ExpArr(i, j)
				ExpArr(j, i) = ExpArr(i, j) 
			end if
			if ( dabs(ExpArr(i, j)) .lt. 1d-10 ) then
				ExpArr(i, j) = 0d0
				ExpArr(j, i) = 0d0
			end if
			! TODO: potential calculation for g-orbitals
!			 if ( PGTOs(i)%shelltype .ge. 5 ) then
!				 ExpArr(i, j) = 0d0
!				 ExpArr(j, i) = 0d0
!			 end if
		end do
	end do
	!$OMP END DO NOWAIT
		
	allocate(TempElecPot(NLGP))
	do i = 1, NLGP
		TempElecPot(i)%x = GridPoints(i)%x + EDyX
		TempElecPot(i)%y = GridPoints(i)%y + EDyY
		TempElecPot(i)%z = GridPoints(i)%z + EDyZ
		TempElecPot(i)%ElecPotVal = 0d0
	end do
	
	!$OMP BARRIER
	
	!$OMP SINGLE
	call system_clock(t1, crate, cmax)
	!$OMP END SINGLE
	
	 !$OMP DO PRIVATE(ic, jc, i, j, PA, PB, PAB, C, sssval, ssssval, cpr2, gamma, dxyz, ii, vi, vj, vk, GVal) SCHEDULE(DYNAMIC)
	do ic = 1, NumCGTO
		PA = (/ CGTOs(ic)%ox, CGTOs(ic)%oy, CGTOs(ic)%oz /)
		do jc = 1, ic
			! check for orbital overlap
			if ( abs(D_CGTO(jc, ic)) .gt. 1d-10 ) then
				PB = (/ CGTOs(jc)%ox, CGTOs(jc)%oy, CGTOs(jc)%oz /)
				do i = CGTOs(ic)%npgto_beg, CGTOs(ic)%npgto_beg + CGTOs(ic)%npgto - 1
					do j = CGTOs(jc)%npgto_beg, CGTOs(jc)%npgto_beg + CGTOs(jc)%npgto - 1
						if ( j .le. i ) then
							if ( ExpArr(j, i) .ne. 0d0 ) then

								gamma = PGTOs(i)%coeff_alpha + PGTOs(j)%coeff_alpha
								PAB(1) = ( PGTOs(i)%coeff_alpha * PA(1) + PGTOs(j)%coeff_alpha * PB(1) ) / gamma
								PAB(2) = ( PGTOs(i)%coeff_alpha * PA(2) + PGTOs(j)%coeff_alpha * PB(2) ) / gamma
								PAB(3) = ( PGTOs(i)%coeff_alpha * PA(3) + PGTOs(j)%coeff_alpha * PB(3) ) / gamma
								
								! check for angular momentum > 0 
								if ( CValMaxDim(1, j, i) + CValMaxDim(2, j, i) + CValMaxDim(3, j, i) .eq. 0 ) then
								
									! loop over all grid points
									do ii = 1, NLGP
										
										cpr2 = (PAB(1)-TempElecPot(ii)%x)**2 + (PAB(2)-TempElecPot(ii)%y)**2 + (PAB(3)-TempElecPot(ii)%z)**2
										if ( dabs(cpr2) .le. 1d-10 ) then
										    cpr2 = 1d-8
										end if
										sssval = ExpArr(j, i) * dsqrt(pi/(4d0*gamma*cpr2))*derf(dsqrt(gamma*cpr2))
#if ( __GFORTRAN__ && __GNUC__ < 5 )
										!$OMP CRITICAL
#else
										!$OMP ATOMIC
#endif
										TempElecPot(ii)%ElecPotVal = TempElecPot(ii)%ElecPotVal + sssval
#if ( __GFORTRAN__ && __GNUC__ < 5 )
										!$OMP END CRITICAL
#else
										!$OMP END ATOMIC
#endif
									end do
									
								else
									! loop over all grid points
									do ii = 1, NLGP
										C = (/ TempElecPot(ii)%x, TempElecPot(ii)%y, TempElecPot(ii)%z /)
										cpr2 = (PAB(1)-C(1))**2 + (PAB(2)-C(2))**2 + (PAB(3)-C(3))**2
										sssval = 0d0
										GVal = 0d0
										
										! most time consuming step
										do dxyz = 1, 3
											do k = 0, CValMaxDim(dxyz, j, i)
												GVal(k+1, dxyz) = G_binom_factors(k, int(L1Arr(dxyz, i), 8), int(L1Arr(dxyz, j), 8), PA(dxyz), PB(dxyz), C(dxyz), PAB(dxyz), gamma)
											end do
										end do
					
										do vi = 0, CValMaxDim(1, j, i)
											do vj = 0, CValMaxDim(2, j, i)
												do vk = 0, CValMaxDim(3, j, i)
													ssssval = GVal(vi+1, 1) * GVal(vj+1, 2) * GVal(vk+1, 3)
													if ( abs(ssssval) .gt. 1d-09 ) then
														sssval = sssval + ssssval * boys_func(vi+vj+vk, gamma*cpr2)
													end if
												end do
											end do
										end do
										ssssval = sssval * ExpArr(j, i)
										
#if ( __GFORTRAN__ && __GNUC__ < 5 )
										!$OMP CRITICAL
#else
										!$OMP ATOMIC
#endif
										TempElecPot(ii)%ElecPotVal = TempElecPot(ii)%ElecPotVal + ssssval
#if ( __GFORTRAN__ && __GNUC__ < 5 )
										!$OMP END CRITICAL
#else
										!$OMP END ATOMIC
#endif
									end do
								end if
							end if
						end if
					end do
				end do
			end if
		end do

#if ( __GFORTRAN__ && __GNUC__ < 5 )
		!$OMP CRITICAL
#else
		!$OMP ATOMIC
#endif
		cprogress = cprogress + 1		
#if ( __GFORTRAN__ && __GNUC__ < 5 )
		!$OMP END CRITICAL
#else
		!$OMP END ATOMIC
#endif
		
		if ( modulo(cprogress, progress) .eq. 0 ) then
			write(*,'(41X,I6,A)') nint(10d0*dble(cprogress)/dble(NumCGTO))*10, " %"
		end if
		
	end do
	!$OMP END DO NOWAIT	
	!$OMP BARRIER
	
	! copy the values back to the main GP array
	do ii = 1, NLGP
		!$OMP CRITICAL
		GridPoints(ii)%ElecPotVal = GridPoints(ii)%ElecPotVal + TempElecPot(ii)%ElecPotVal
		!$OMP END CRITICAL
	end do
	
	deallocate(TempElecPot)
	!$OMP END PARALLEL
	
	call system_clock(t2, crate, cmax)
	
	write(*,*)
	write(*,'(A55,F14.2,A)') " time spent: ", dble((t2-t1)/crate), " s"
	write(*,*)
	call write_done

	deallocate(L1Arr)
	deallocate(CValMaxDim)
	deallocate(ExpArr)

end subroutine

function calc_anal_elec_pot_point(ox, oy, oz) result(res)
	use data_mo
	use global_c
	use data_grid
	implicit none
!
	double precision, intent(in) :: ox, oy, oz
	double precision :: ssval, sssval, ssssval, gamma, r2, cpr2, prefac, pprefac, ppprefac, pppprefac
	double precision, dimension(3) :: PAB, PA, PB, C
	integer :: k, dxyz, l1, l2, vi, vj, vk, pgto_get_l
	double precision, dimension(7, 3) :: GVal
	double precision :: G_binom_factors, boys_func
	integer, dimension(3) :: CValMaxDim
	integer :: i, j, ic, jc
	double precision :: res
!
	res = 0d0

	do ic = 1, NumCGTO
		do jc = 1, ic
			if ( abs(D_CGTO(ic, jc)) .gt. 1d-8 ) then
				
				r2 = AtomsDistMatSq(CGTOs(ic)%atomid, CGTOs(jc)%atomid)
				PA = (/ CGTOs(ic)%ox, CGTOs(ic)%oy, CGTOs(ic)%oz /)
				PB = (/ CGTOs(jc)%ox, CGTOs(jc)%oy, CGTOs(jc)%oz /)
				
				if ( ic .ne. jc ) then
					prefac = 4d0 * D_CGTO(ic, jc) 
				else
					prefac = 2d0 * D_CGTO(ic, jc)
				end if
				
				do i = CGTOs(ic)%npgto_beg, CGTOs(ic)%npgto_beg + CGTOs(ic)%npgto - 1
				
					pprefac = prefac * PGTOs(i)%coeff_d						 
					if ( abs(pprefac) .gt. 1d-8 ) then
					
						do j = CGTOs(jc)%npgto_beg, CGTOs(jc)%npgto_beg + CGTOs(jc)%npgto - 1
						
								ppprefac = pprefac * PGTOs(j)%coeff_d 
						
								gamma = PGTOs(i)%coeff_alpha + PGTOs(j)%coeff_alpha
								ssval = 2d0*pi/gamma * dexp( -PGTOs(i)%coeff_alpha*PGTOs(j)%coeff_alpha * r2 / gamma )
								pppprefac = ppprefac * ssval
								
								if ( abs(pppprefac) .gt. 1d-10 ) then
								
									do dxyz = 1, 3
										PAB(dxyz) = ( PGTOs(i)%coeff_alpha * PA(dxyz) + PGTOs(j)%coeff_alpha * PB(dxyz) ) / gamma
									end do
									
									
										C  = (/ ox+EDyX, oy+EDyY, oz+EDyZ /)
										CValMaxDim = 0
										cpr2 = 0d0
										GVal = 0d0
										do dxyz = 1, 3
											l1 = pgto_get_l(PGTOs(i)%shelltype, PGTOs(i)%subtype, dxyz)
											l2 = pgto_get_l(PGTOs(j)%shelltype, PGTOs(j)%subtype, dxyz)
											cpr2 = cpr2 + ( PAB(dxyz) - C(dxyz) )**2
											
											do k = 0, l1+l2
												GVal(k+1, dxyz) = G_binom_factors(k, l1, l2, PA(dxyz), PB(dxyz), C(dxyz), PAB(dxyz), gamma)
												if ( CValMaxDim(dxyz) .lt. k ) then
													CValMaxDim(dxyz) = k
												end if
											end do
											
										end do
										
										sssval = 0d0
					
										do vi = 0, CValMaxDim(1)
											do vj = 0, CValMaxDim(2)
												do vk = 0, CValMaxDim(3)
													ssssval = GVal(vi+1, 1) * GVal(vj+1, 2) * GVal(vk+1, 3)
													if ( abs(ssssval) .gt. 1d-20 ) then
														sssval = sssval + ssssval * boys_func(vi+vj+vk, gamma * cpr2)
													end if
												end do
											end do
										end do
										
										res = res + pppprefac * sssval
										
								end if
							
!							 end if
							
						end do
					
					end if
					
				end do
			end if
		end do
	end do	
	

end function

function calc_anal_nuc_pot_point(ox, oy, oz) result(res)
	use data_mo
	use global_c
	use data_grid
	implicit none
!
	double precision, intent(in) :: ox, oy, oz
	integer :: j
	double precision :: res, dist
!
	res = 0d0

	do j = 1, NA
		
		if ( ( Atoms(j)%ghost .eqv. .FALSE. ) .AND. ( Atoms(j)%skip .eqv. .FALSE. ) ) then

			dist = dsqrt( (ox - Atoms(j)%x + EDyX )**2 + (oy - Atoms(j)%y + EDyY )**2 + (oz - Atoms(j)%z + EDyZ )**2 )
				
			if ( dist .gt. 1.0d-5 ) then
				res = res + real(Atoms(j)%charge) / dist				
			end if
			
		end if
		
	end do
	
end function

subroutine calc_ligand_dens
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	integer :: i
!
	write(*,*) "> CALCULATING ELECTRON DENSITY ON GRID POINTS ..."
	write(*,*)
	
	!$OMP PARALLEL DO SCHEDULE(dynamic)
	do i = 1, NLGP		
		call calc_density_point_dmat4(GridPoints(i)%x+EDyX, GridPoints(i)%y+EDyY, GridPoints(i)%z+EDyZ, GridPoints(i)%LDens)		
	end do
	!$OMP END PARALLEL DO 
	
	write(*,*)
	
end subroutine

subroutine calc_ligand_ldax
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	integer :: i
	double precision :: XVal 
!
	write(*,*) "> CALCULATING SLATER-DIRAC EXCHANGE CONTRIBUTION ..."
	write(*,*)
	
	do i = 1, NLGP
		XVal = - (3d0/Pi)**(1d0/3d0) * (GridPoints(i)%LDens)**(1d0/3d0)
		GridPoints(i)%ElecPotVal = GridPoints(i)%ElecPotVal + (-XVal)
	end do
	
	call write_done
	write(*,*)
	
end subroutine

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

! *************************************************************************
! * calculates the overlap integral of a PGTO
! *************************************************************************
function calc_olint_pgto(i, j, c) result(val)
	use data_mo
	use global_c
	implicit none
!
	integer, intent(in) :: i, j
	integer, dimension(3), intent(in) :: c
	double precision :: val, gamma, r2, sval
	integer :: k, l1, l2, dxyz
	integer :: pgto_get_l
	integer (8) :: n_fac2
	double precision, dimension(3) :: PAB, PA, PB
	double precision :: f_binom_factors
!
	val = 0d0
	
	if ( ( i .gt. NumPGTO ) .OR. ( j .gt. NumPGTO ) .OR. ( i .lt. 1 ) .OR. ( j .lt. 1 ) ) then
		write(*,*) "Internal Error!"
		stop 
	end if
	
	gamma = PGTOs(i)%coeff_alpha + PGTOs(j)%coeff_alpha
	r2 = AtomsDistMatSq(PGTOs(i)%atomid, PGTOs(j)%atomid)
	PAB(1) = ( PGTOs(i)%coeff_alpha * PGTOs(i)%ox + PGTOs(j)%coeff_alpha * PGTOs(j)%ox ) / gamma
	PAB(2) = ( PGTOs(i)%coeff_alpha * PGTOs(i)%oy + PGTOs(j)%coeff_alpha * PGTOs(j)%oy ) / gamma
	PAB(3) = ( PGTOs(i)%coeff_alpha * PGTOs(i)%oz + PGTOs(j)%coeff_alpha * PGTOs(j)%oz ) / gamma
	
	PA = (/ PGTOs(i)%ox, PGTOs(i)%oy, PGTOs(i)%oz /)
	PB = (/ PGTOs(j)%ox, PGTOs(j)%oy, PGTOs(j)%oz /)
	
	val = (Pi/gamma)**(3d0/2d0) * dexp(-PGTOs(i)%coeff_alpha*PGTOs(j)%coeff_alpha*r2/gamma)
	
	do dxyz = 1, 3
	
		l1 = pgto_get_l(PGTOs(i)%shelltype, PGTOs(i)%subtype, dxyz)
		l2 = pgto_get_l(PGTOs(j)%shelltype, PGTOs(j)%subtype, dxyz) + c(dxyz)
		
		sval = 0d0
		do k = 0, (l1+l2)/2
			sval = sval + f_binom_factors(2*k, l1, l2, PAB(dxyz)-PA(dxyz), PAB(dxyz)-PB(dxyz)) * dble(n_fac2(2*k-1)) / (2d0*gamma)**(k)
		end do
		val = val * sval
	
	end do
	
end function


! *************************************************************************
! * calculates the overlap matrix
! *************************************************************************
subroutine calc_overlap_matrix(silent)
	use data_mo
	implicit none
!
	logical, intent(in) :: silent
	integer :: i, j, k, k1, k2
	double precision :: calc_olint_pgto
!
	S_PGTO = 0d0
	S_CGTO = 0d0
	D_CGTO = 0d0
	CSThresD = 0d0
	TotDensDS = 0d0
	
	if ( silent .eqv. .FALSE. ) then
		write(*,'(5X,A)') "> CALCULATING OVERLAP MATRIX..."
		write(*,*)
	end if
	
	!$OMP PARALLEL
	
	!$OMP DO PRIVATE(j)
	do i = 1, NumPGTO
		S_PGTO(i, i) = calc_olint_pgto(i, i, (/ 0, 0, 0 /))
		do j = i + 1, NumPGTO
			S_PGTO(i, j) = calc_olint_pgto(i, j, (/ 0, 0, 0 /))
			S_PGTO(j, i) = S_PGTO(i, j)
		end do
	end do
	!$OMP END DO NOWAIT
	
	!$OMP BARRIER
	
	!$OMP DO PRIVATE(i, k2, j)
	do k1 = 1, NumPGTO
		i = PGTOs(k1)%cgtoid
		do k2 = 1, NumPGTO
			j = PGTOs(k2)%cgtoid
#if ( __GFORTRAN__ && __GNUC__ < 5 )
			!$OMP CRITICAL
#else
			!$OMP ATOMIC
#endif
			S_CGTO(i, j) = S_CGTO(i, j) + PGTOs(k1)%coeff_d * PGTOs(k2)%coeff_d * S_PGTO(k1, k2)
#if ( __GFORTRAN__ && __GNUC__ < 5 )
			!$OMP END CRITICAL
#else
			!$OMP END ATOMIC
#endif
		end do
	end do
	!$OMP END DO NOWAIT
	
	!$OMP BARRIER
	
	!$OMP SINGLE
	if ( silent .eqv. .FALSE. ) then
		write(*,'(5X,A)') "> CALCULATING DENSITY MATRIX..."
		write(*,*)
	end if
	!$OMP END SINGLE

	!$OMP BARRIER
	
	!$OMP DO PRIVATE(i, j, k, k1, k2) SCHEDULE(dynamic)
	do i = 1, NumPGTO
		if ( PGTOs(i)%deleted .eqv. .FALSE. ) then
			k1 = PGTOs(i)%cgtoid
			do j = 1, NumPGTO
				if ( PGTOs(j)%deleted .eqv. .FALSE. ) then
					k2 = PGTOs(j)%cgtoid
					do k = 1, NMO
						if ( MOs(k)%occup .gt. 0d0 ) then
#if ( __GFORTRAN__ && __GNUC__ < 5 )
							!$OMP CRITICAL
#else
							!$OMP ATOMIC
#endif
							D_CGTO(k1, k2) = D_CGTO(k1, k2) + MOs(k)%occup * 0.5 * MOs(k)%mocoeff(k1) * MOs(k)%mocoeff(k2) / &
								 ( (dble(PGTOs(i)%cgtonpgto)) * ( dble(PGTOs(j)%cgtonpgto)) )
#if ( __GFORTRAN__ && __GNUC__ < 5 )
							!$OMP END CRITICAL
#else
							!$OMP END ATOMIC
#endif
						end if
					end do
				end if
			end do
		end if
	end do
	!$OMP END DO NOWAIT

	!$OMP BARRIER

	!$OMP DO PRIVATE(j)
	do i = 1, NumCGTO
		do j = 1, i
			if ( dabs(D_CGTO(i, j)) .gt. CSThresD ) then
#if ( __GFORTRAN__ && __GNUC__ < 5 )
				!$OMP CRITICAL
#else
				!$OMP ATOMIC WRITE
#endif
				CSThresD = dabs(D_CGTO(i, j))
#if ( __GFORTRAN__ && __GNUC__ < 5 )
				!$OMP END CRITICAL
#else
				!$OMP END ATOMIC 
#endif
			end if
		end do
	end do
	!$OMP END DO
	
	!$OMP SINGLE
	CSThresI = CSThresC * CSThresD
	!$OMP END SINGLE
	
	!$OMP DO REDUCTION(+:TotDensDS) PRIVATE(j)
	do i = 1, NumCGTO
		do j = 1, NumCGTO
			TotDensDS = TotDensDS + S_CGTO(i, j) * D_CGTO(i, j) 
		end do
	end do
	!$OMP END DO NOWAIT
	
	!$OMP BARRIER
	
	!$OMP END PARALLEL
	
	! needs "-heap-arrays" or "ulimit -s unlimited"
	DSMat = matmul(D_CGTO, S_CGTO)
	
	write(*,*)
	
end subroutine

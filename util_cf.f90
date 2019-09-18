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

subroutine init_lapack
	use global_c
	implicit none
!
	! initialization of LAPACK
	NArraySize = 16
	NArraySizeLDA = 16
	LAPACKJobType = 'V'
	LAPACKMatrixType = 'L'
	LFMatrix = cmplx(0d0, 0d0)
	
	! get the optimal work array size
	RWORKSize = max(1, 3*NArraySize-2)
	allocate(RWorkArr(RWORKSize))
	LWORKSize = max(1, 2*NArraySize-1)
	allocate(LWorkArr(LWORKSize))
	LWORKSize = -1
	
	call zheev(LAPACKJobType, LAPACKMatrixType, NArraySize, LFMatrix, NArraySizeLDA, EigenValues, LWorkArr, LWORKSize,  RWorkArr, lapackstatus)
	if ( lapackstatus .ne. 0 ) then
		if ( lapackstatus .ne. 0 ) then
			write(*,*) "LAPACK Error! Error Code = ", lapackstatus
			stop
		end if
	end if
	
	LWORKSize = int(real(LWorkArr(1)))
	deallocate(LWorkArr)
	allocate(LWorkArr(LWORKSize))
	
	LAPACKInit = .TRUE.
	
end subroutine

subroutine init_cf_matrix_values
	use data_mo
	use global_c
	use data_grid
	implicit none
!
	integer :: i, k
	double precision, dimension(1:16) :: Jz
	double precision, dimension(1:6, 1:16) :: JP, JM
	double precision :: JJPO, J
	double precision :: Radial_Expectation_Value_4f
!
	! values for dysprosium
	J = 15d0/2d0
	JJPO = J*(J+1d0)
	
	! computes Jz values
	Jz = 0d0
	do i = 1, 16
		Jz(i) = ( dble(i) ) - J - 1d0
	end do
	
	! computes J+ and J- values
	JP = 0d0
	JM = 0d0
	
	do k = 1, 6
		do i = 1, 16
			if ( i .le. 16-k ) then
				JP(k, i) = dsqrt((dble(J) - Jz(i) - dble(k-1)) * (dble(J)+Jz(i)+1d0 + dble(k-1)))
				if ( k .ne. 1 ) then
					JP(k, i) = JP(k, i) * JP(k-1, i) 
				end if
			end if
			if ( i .ge. 1+k ) then
				JM(k, i) = dsqrt((dble(J) + Jz(i) - dble(k-1)) * (dble(J)-Jz(i)+1d0 + dble(k-1)))
				if ( k .ne. 1 ) then
					JM(k, i) = JM(k, i) * JM(k-1, i) 
				end if
			end if
		end do
	end do
	
	! computes ESO matrix elements
	StevOp2Mat = cmplx(0.0, 0.0)
	StevOp4Mat = cmplx(0.0, 0.0)
	StevOp6Mat = cmplx(0.0, 0.0)
	
	do i = 1, 16
		! k = 2, q = 0
		StevOp2Mat(0, i, i) = 3d0 * Jz(i)**2 - JJPO
		! k = 2, q = 1, -1
		if ( i + 1 .le. 16 ) then
			StevOp2Mat( 1, i+1, i) = 0.25d0 * ( Jz(i)*JP(1,i) + JP(1,i)*Jz(i+1) ) 
			StevOp2Mat(-1, i+1, i) = - cmplx(0.0d0, 1.0d0) * 0.25d0 * ( Jz(i)*JP(1,i) + JP(1,i)*Jz(i+1) ) 
		end if
		if ( i - 1 .ge. 1 ) then
			StevOp2Mat( 1, i-1, i) = 0.25d0 * ( Jz(i)*JM(1,i) + JM(1,i)*Jz(i-1) ) 
			StevOp2Mat(-1, i-1, i) = - cmplx(0.0d0, 1.0d0) * 0.25d0 * ( -Jz(i)*JM(1,i) - JM(1,i)*Jz(i-1) ) 
		end if
		! k = 2, q = 2, -2
		if ( i + 2 .le. 16 ) then
			StevOp2Mat( 2, i+2, i) = 0.50d0 * JP(2,i)
			StevOp2Mat(-2, i+2, i) = cmplx(0.0d0, -1.0d0) / 2d0 * (JP(2,i))
		
		end if
		if ( i - 2 .ge. 1 ) then
			StevOp2Mat( 2, i-2, i) = 0.50d0 * JM(2,i)
			StevOp2Mat(-2, i-2, i) = cmplx(0.0d0, 1.0d0) / 2d0 * (JM(2,i))
		end if
		
		! k = 4, q = 0
		StevOp4Mat(0, i, i) = 35d0 * Jz(i)**4 - (30d0*JJPO - 25d0) * Jz(i)**2 + 3d0*JJPO**2 - 6d0*JJPO
		! k = 4, q = 1, -1
		if ( i + 1 .le. 16 ) then
			StevOp4Mat( 1, i+1, i) = 0.25d0 * ( (7d0*Jz(i)**3 - (3d0*JJPO + 1d0)*Jz(i))*JP(1,i) + &
							JP(1,i)*(7d0*Jz(i+1)**3 - (3d0*JJPO + 1d0)*Jz(i+1) ) )
			StevOp4Mat(-1, i+1, i) = - cmplx(0.0d0, 1.d0) / 4d0 * ( (7d0*Jz(i)**3 - (3d0*JJPO + 1d0)*Jz(i))*JP(1,i) + &
							JP(1,i)*(7d0*Jz(i+1)**3 - (3d0*JJPO + 1d0)*Jz(i+1) ) )
		end if
		if ( i - 1 .ge. 1 ) then
			StevOp4Mat( 1, i-1, i) = 0.25d0 * ( (7d0*Jz(i)**3 - (3d0*JJPO + 1d0)*Jz(i))*JM(1,i) + &
							JM(1,i)*(7d0*Jz(i-1)**3 - (3d0*JJPO + 1d0)*Jz(i-1) ) )
			StevOp4Mat(-1, i-1, i) = cmplx(0.0d0, 1.d0) / 4d0 * ( (7d0*Jz(i)**3 - (3d0*JJPO + 1d0)*Jz(i))*JM(1,i) + &
							JM(1,i)*(7d0*Jz(i-1)**3 - (3d0*JJPO + 1d0)*Jz(i-1) ) )
		end if
		! k = 4, q = 2, -2
		if ( i + 2 .le. 16 ) then
			StevOp4Mat( 2, i+2, i) = 0.25d0 * ( JP(2,i) * ( 7d0*Jz(i+2)**2 - JJPO - 5d0 ) + &
							( 7d0*Jz(i)**2 - JJPO - 5d0 ) * JP(2,i) )
			StevOp4Mat(-2, i+2, i) = - cmplx(0.0d0, 1.d0) / 4d0 * ( JP(2,i) * ( 7d0*Jz(i+2)**2 - JJPO - 5d0 ) + &
							( 7d0*Jz(i)**2 - JJPO - 5d0 ) * JP(2,i) )
		end if
		if ( i - 2 .ge. 1 ) then
			StevOp4Mat( 2, i-2, i) = 0.25d0 * ( JM(2,i) * ( 7d0*Jz(i-2)**2 - JJPO - 5d0 ) + &
							( 7d0*Jz(i)**2 - JJPO - 5d0 ) * JM(2,i) )
			StevOp4Mat(-2, i-2, i) = cmplx(0.0d0, 1.d0) / 4d0 * ( JM(2,i) * ( 7d0*Jz(i-2)**2 - JJPO - 5d0 ) + &
							( 7d0*Jz(i)**2 - JJPO - 5d0 ) * JM(2,i) )
		end if
		! k = 4, q = 3, -3
		if ( i + 3 .le. 16 ) then
			StevOp4Mat( 3, i+3, i) = 0.25d0 * ( JP(3,i) * Jz(i+3) + Jz(i) * JP(3,i) )
			StevOp4Mat(-3, i+3, i) = - cmplx(0.0d0, 1.d0) / 4d0 * ( JP(3,i) * Jz(i+3) + Jz(i)*JP(3,i) )
		end if
		if ( i - 3 .ge. 1 ) then
			StevOp4Mat( 3, i-3, i) = 0.25d0 * ( JM(3,i) * Jz(i-3) + Jz(i) * JM(3,i) )
			StevOp4Mat(-3, i-3, i) = cmplx(0.0d0, 1.d0) / 4d0 * ( JM(3,i) * Jz(i-3) + Jz(i)*JM(3,i) )
		end if
		! k = 4, q = 4, -4
		if ( i + 4 .le. 16 ) then
			StevOp4Mat( 4, i+4, i) = 0.50d0 * JP(4,i)
			StevOp4Mat(-4, i+4, i) = - cmplx(0.0d0, 1.d0) / 2d0 * JP(4,i)
		end if
		if ( i - 4 .ge. 1 ) then
			StevOp4Mat( 4, i-4, i) = 0.50d0 * JM(4,i)
			StevOp4Mat(-4, i-4, i) = cmplx(0.0d0, 1.d0) / 2d0 * JM(4,i)
		end if
		
		! k = 6, q = 0
		StevOp6Mat(0, i, i) = 231d0 * Jz(i)**6 - (315d0*JJPO - 735d0) * Jz(i)**4 + (105d0*JJPO**2 - 525d0*JJPO + 294d0)*Jz(i)**2 - 5d0*JJPO**3 + 40d0*JJPO**2 - 60d0*JJPO
		! k = 6, q = 1, -1
		if ( i + 1 .le. 16 ) then
			StevOp6Mat( 1, i+1, i) = 0.25d0 * ( JP(1,i)*(33d0*Jz(i+1)**5-(30d0*JJPO-15d0)*Jz(i+1)**3+(5d0*JJPO**2-10d0*JJPO+12d0)*Jz(i+1)) + (33d0*Jz(i)**5-(30d0*JJPO-15d0)*Jz(i)**3+(5d0*JJPO**2-10d0*JJPO+12d0)*Jz(i))*JP(1,i))
			StevOp6Mat(-1, i+1, i) = -cmplx(0.0d0, 1.d0) * StevOp6Mat( 1, i+1, i)
		end if
		if ( i - 1 .ge. 1 ) then
			StevOp6Mat( 1, i-1, i) = 0.25d0 * ( JM(1,i)*(33d0*Jz(i-1)**5-(30d0*JJPO-15d0)*Jz(i-1)**3+(5d0*JJPO**2-10d0*JJPO+12d0)*Jz(i-1)) + (33d0*Jz(i)**5-(30d0*JJPO-15d0)*Jz(i)**3+(5d0*JJPO**2-10d0*JJPO+12d0)*Jz(i))*JM(1,i))
			StevOp6Mat(-1, i-1, i) = cmplx(0.0d0, 1.d0) * StevOp6Mat( 1, i-1, i)
		end if
		! k = 6, q = 2, -2
		if ( i + 2 .le. 16 ) then
			StevOp6Mat( 2, i+2, i) = 0.25d0 * ( JP(2,i)*(33d0*Jz(i+2)**4-(18d0*JJPO+123d0)*Jz(i+2)**2+JJPO**2+10d0*JJPO+102d0) + (33d0*Jz(i)**4-(18d0*JJPO+123d0)*Jz(i)**2+JJPO**2+10d0*JJPO+102d0)*JP(2,i))
			StevOp6Mat(-2, i+2, i) = -cmplx(0.0d0, 1.d0)/4d0 * ( JP(2,i)*(33d0*Jz(i+2)**4-(18d0*JJPO+123d0)*Jz(i+2)**2+JJPO**2+10d0*JJPO+102d0) + (33d0*Jz(i)**4-(18d0*JJPO+123d0)*Jz(i)**2+JJPO**2+10d0*JJPO+102d0)*JP(2,i))
		end if
		if ( i - 2 .ge. 1 ) then
			StevOp6Mat( 2, i-2, i) = 0.25d0 * ( JM(2,i)*(33d0*Jz(i-2)**4-(18d0*JJPO+123d0)*Jz(i-2)**2+JJPO**2+10d0*JJPO+102d0) + (33d0*Jz(i)**4-(18d0*JJPO+123d0)*Jz(i)**2+JJPO**2+10d0*JJPO+102d0)*JM(2,i))
			StevOp6Mat(-2, i-2, i) = cmplx(0.0d0, 1.d0)/4d0 * ( JM(2,i)*(33d0*Jz(i-2)**4-(18d0*JJPO+123d0)*Jz(i-2)**2+JJPO**2+10d0*JJPO+102d0) + (33d0*Jz(i)**4-(18d0*JJPO+123d0)*Jz(i)**2+JJPO**2+10d0*JJPO+102d0)*JM(2,i))
		end if
		! k = 6, q = 3, -3
		if ( i + 3 .le. 16 ) then
			StevOp6Mat( 3, i+3, i) = 0.25d0 * ( JP(3,i)*(11*Jz(i+3)**3-(3d0*JJPO+59d0)*Jz(i+3)) + (11d0*Jz(i)**3-(3d0*JJPO+59d0)*Jz(i))*JP(3,i))
			StevOp6Mat(-3, i+3, i) = -cmplx(0.0d0, 1.d0)/4d0 * ( JP(3,i)*(11*Jz(i+3)**3-(3d0*JJPO+59d0)*Jz(i+3)) + (11d0*Jz(i)**3-(3d0*JJPO+59d0)*Jz(i))*JP(3,i))
		end if
		if ( i - 3 .ge. 1 ) then
			StevOp6Mat( 3, i-3, i) = 0.25d0 * ( JM(3,i)*(11d0*Jz(i-3)**3-(3d0*JJPO+59d0)*Jz(i-3)) + (11d0*Jz(i)**3-(3d0*JJPO+59d0)*Jz(i))*JM(3,i))
			StevOp6Mat(-3, i-3, i) = cmplx(0.0d0, 1.d0)/4d0 * ( JM(3,i)*(11*Jz(i-3)**3-(3d0*JJPO+59d0)*Jz(i-3)) + (11d0*Jz(i)**3-(3d0*JJPO+59d0)*Jz(i))*JM(3,i))
		end if
		! k = 6, q = 4, -4
		if ( i + 4 .le. 16 ) then
			StevOp6Mat( 4, i+4, i) = 0.25d0 * (JP(4,i)*(11d0*Jz(i+4)**2-JJPO-38d0) + (11d0*Jz(i)**2-JJPO-38d0)*JP(4,i))
			StevOp6Mat(-4, i+4, i) = -cmplx(0.0d0, 1.d0)/4d0 * (JP(4,i)*(11*Jz(i+4)**2-JJPO-38d0) + (11d0*Jz(i)**2-JJPO-38d0)*JP(4,i))
		end if
		if ( i - 4 .ge. 1 ) then
			StevOp6Mat( 4, i-4, i) = 0.25d0 * (JM(4,i)*(11d0*Jz(i-4)**2-JJPO-38d0) + (11d0*Jz(i)**2-JJPO-38d0)*JM(4,i))
			StevOp6Mat(-4, i-4, i) = cmplx(0.0d0, 1.d0) / 4d0 * (JM(4,i)*(11*Jz(i-4)**2-JJPO-38d0) + (11d0*Jz(i)**2-JJPO-38d0)*JM(4,i))
		end if
		! k = 6, q = 5, -5
		if ( i + 5 .le. 16 ) then
			StevOp6Mat( 5, i+5, i) = 0.25d0 * (JP(5,i)*Jz(i+5) + Jz(i)*JP(5,i))
			StevOp6Mat(-5, i+5, i) = - cmplx(0.0d0, 1.d0) / 4d0 * (JP(5,i)*Jz(i+5) + Jz(i)*JP(5,i))
		end if
		if ( i - 5 .ge. 1 ) then
			StevOp6Mat( 5, i-5, i) = 0.25d0 * (JM(5,i)*Jz(i-5) + Jz(i)*JM(5,i))
			StevOp6Mat(-5, i-5, i) = cmplx(0.0d0, 1.d0) / 4d0 * (JM(5,i)*Jz(i-5) + Jz(i)*JM(5,i))
		end if
		! k = 6, q = 6, -6
		if ( i + 6 .le. 16 ) then
			StevOp6Mat( 6, i+6, i) = 0.50d0 * JP(6,i)
			StevOp6Mat(-6, i+6, i) = - cmplx(0.0d0, 1.d0) / 2d0 * JP(6,i)
		end if
		if ( i - 6 .ge. 1 ) then
			StevOp6Mat( 6, i-6, i) = 0.50d0 * JM(6,i)
			StevOp6Mat(-6, i-6, i) = cmplx(0.0d0, 1.d0) / 2d0 * JM(6,i)
		end if
		
	end do
	
	! prefactors from tesseral harmonics
	StevPLM = 0d0
	
	StevPLM(2, 2) = 0.25d0 * dsqrt(15d0/pi)
	StevPLM(2, 1) = 0.50d0 * dsqrt(15d0/pi)
	StevPLM(2, 0) = 0.25d0 * dsqrt(5d0/pi)
	
	StevPLM(4, 4) = 0.1875d0 * dsqrt(35d0/pi)
	StevPLM(4, 3) = 0.3750d0 * dsqrt(70d0/pi)
	StevPLM(4, 2) = 0.3750d0 * dsqrt(5d0/pi)
	StevPLM(4, 1) = 0.7500d0 * dsqrt(5d0/(2d0*pi))
	StevPLM(4, 0) = 0.1875d0 * dsqrt(1d0/pi)
	
	StevPLM(6, 6) = 3.609375d0 * dsqrt(26d0/(231d0*pi))
	StevPLM(6, 5) = dsqrt(9009d0/(512d0*pi))
	StevPLM(6, 4) = 0.656250d0 * dsqrt(13d0/(7d0*pi))
	StevPLM(6, 3) = 0.031250d0 * dsqrt(2730d0/pi)
	StevPLM(6, 2) = 0.015625d0 * dsqrt(2730d0/pi)
	StevPLM(6, 1) = 0.125000d0 * dsqrt(273d0/(4d0*pi))
	StevPLM(6, 0) = 0.031250d0 * dsqrt(13d0/pi)
	
	! conversion factors between Stevens CFPs and Wybourne CFPs
	! M. Karbowiak, C. Rudowicz, J. Comput. Chem. 2014, 35, 1935--1941
	StevLam = 0d0

	StevLam(2, 2) = 2d0 / dsqrt(6d0)
	StevLam(2, 1) = 1d0 / dsqrt(6d0)
	StevLam(2, 0) = 2d0
! 
	StevLam(4, 4) = 8d0 / dsqrt(70d0)
	StevLam(4, 3) = 2d0 / dsqrt(35d0)
	StevLam(4, 2) = 4d0 / dsqrt(10d0)
	StevLam(4, 1) = 2d0 / dsqrt(5d0)
	StevLam(4, 0) = 8d0
! 
	StevLam(6, 6) = 16d0 / dsqrt(231d0)
	StevLam(6, 5) = 8d0 / (3d0 * dsqrt(77d0))
	StevLam(6, 4) = 16d0 / (3d0 * dsqrt(14d0))
	StevLam(6, 3) = 8d0 / dsqrt(105d0)
	StevLam(6, 2) = 16d0 / dsqrt(105d0)
	StevLam(6, 1) = 8d0 / dsqrt(42d0)
	StevLam(6, 0) = 16d0
	
	! theta_k values, 
	! R.J. Elliott, K.W.H. Stevens, Proc. Roy. Soc. A218, 553 (1953).
	StevMultFac = 1d0
	StevMultFac(2) = -2d0/(9d0*35d0)
	StevMultFac(4) = -8d0/(11d0*45d0*273d0)
	StevMultFac(6) =  4d0/(13d0*13d0*11d0*11d0*27d0*7d0)

	! radial integrals
	OldRadInts = 0d0
	OldRadInts(2) = Radial_Expectation_Value_4f(2)
	OldRadInts(4) = Radial_Expectation_Value_4f(4)
	OldRadInts(6) = Radial_Expectation_Value_4f(6)
	
	! Sternheimer shielding parameters for Dy(III)
	SternheimerShieldings = 0d0
	SternheimerShieldings(2) =  0.5270d0
	SternheimerShieldings(4) = -0.0199d0
	SternheimerShieldings(6) = -0.0316d0
	
	CFInit = .TRUE.

end subroutine

function ckq(k, q) result(val)
	implicit none
!
	integer, intent(in) :: k, q
	double precision :: val
!
	val = 0d0
	select case ( k ) 
		case ( 2 ) 
			select case ( q ) 
				case ( 2 )
					val = 2d0 * dsqrt(2d0/3d0) / 1d0
					return
				case ( 1 )
					val = - 2d0 * dsqrt(2d0/3d0) / 2d0
					return
				case ( 0 )
					val = 8d0 / 4d0
					return
			end select
		case ( 4 )
			select case ( q )
				case ( 4 )
					val = (8d0*dsqrt(2d0/35d0)) / 1d0
					return
				case ( 3 )
					val = (-16d0/dsqrt(35d0)) / 4d0
					return
				case ( 2 )
					val = (32d0*dsqrt(2d0/5d0)) / 8d0
					return
				case ( 1 )
					val = (-96d0/dsqrt(5d0)) / 24d0
					return
				case ( 0 )
					val = 384d0 / 48d0
					return
			end select
		case ( 6 )
			select case ( q )
				case ( 6 )
					val = (32d0/dsqrt(231d0)) / 1d0
					return
				case ( 5 )
					val = (-32d0/dsqrt(77d0)) / 6d0
					return
				case ( 4 )
					val = (64d0*dsqrt(2d0/7d0)) / 12d0
					return
				case ( 3 )
					val = (-64d0*dsqrt(15d0/7d0)) / 60d0
					return
				case ( 2 )
					val = ((14625d0-dsqrt(213310669d0))/26d0) / 360d0
					return
				case ( 1 )
					val = ((-10659d0-dsqrt(113889345d0))/6d0) / 1440d0
					return
				case ( 0 )
					val = 46080d0 / 2880d0
					return
			end select
	end select
	
	write(*,*) "ERROR! value not implemented!"
	stop
	
end function

subroutine calc_cf_energy(lrota, lrotb, lrotg, val)
	use global_c
	implicit none
!
	double precision, intent(in) :: lrota, lrotb, lrotg
	double precision, intent(out) :: val
!
	if ( OPCM .eqv. .FALSE. ) then
		call calc_cf_energy_grid(lrota, lrotb, lrotg, val)
	else
		call calc_cf_energy_pcm(lrota, lrotb, lrotg, val)
	end if

end subroutine

subroutine calc_cf_energy_grid(lrota, lrotb, lrotg, val)
	use data_mo
	use data_grid
	use global_c
	implicit none
!
	double precision, intent(in) :: lrota, lrotb, lrotg
	double precision, intent(out) :: val
	integer :: i, k, q, j
	double precision :: srcX, srcY, srcZ, ttheta, tphi, ZdivR
	double complex :: scval, SphericalHarmonicY
!
	if ( LAPACKInit .eqv. .FALSE. ) then
		call init_lapack
	end if
	
	if ( CFInit .eqv. .FALSE. ) then
		call init_cf_matrix_values
	end if
	
	ARkq = 0d0
	BRkq = 0d0
	BWkq = dcmplx(0.0d0, 0.0d0)
	
	! determine Waybourne LFPs
	do k = 2, MaxKRank, 2
		do q = 0, k, 1
			scval = dcmplx(0d0, 0d0)
			do i = 1, NLGP
				if ( ( GridPoints(i)%dist .le. LFTCutOff_O ) .and. ( GridPoints(i)%dist .gt. LFTCutOff_I ) ) then
			
					call rotate_euler_3d_inv(lrota, lrotb, lrotg, GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z, srcX, srcY, srcZ)
					
					! avoid numerical issues
					ZdivR = srcZ / GridPoints(i)%dist
					if ( ZdivR .gt. 1d0 ) then
						ZdivR = 1d0
					else if ( ZdivR .lt. -1d0 ) then
						ZdivR = -1d0
					end if
					
					ttheta = dacos(ZdivR)
					tphi = datan2(srcY, srcX)
					
					scval = scval + dconjg(SphericalHarmonicY(k, q, ttheta, tphi)) * (GridPoints(i)%RadWFdV * (GridPoints(i)%ElecPotVal - GridPoints(i)%NucPotval))
				end if
			end do
			
			BWkq(k, q) = dsqrt((2d0*dble(k)+1d0)/FourPi) * NumberOf4fElectrons * StevMultFac(k) * au2rcm * scval * EnergyScalingFactors(k)
			if ( q .gt. 0 ) then
				BWkq(k,-q) = (-1)**q * dconjg(BWkq(k, q))
			end if
		end do
	end do
		
	! conversion to Stevens LFPs 
	BRkq = 0d0
	do k = 2, MaxKRank, 2
		BRkq(k, 0) = dble(BWkq(k, 0)) / StevLam(k, 0)
		do q = 0, k, 1
			BRkq(k, q) = dble(BWkq(k, q)) / StevLam(k, abs(q))
		end do
		do q = -k, -1, 1
			BRkq(k, q) = dimag(BWkq(k, abs(q))) / StevLam(k, abs(q))
		end do
		do q = -k, k
						ARkq(k, q) = BRkq(k, q) / StevMultFac(k)
		end do
	end do
	
	if ( OSternheimer .eqv. .TRUE. ) then
		do k = 2, MaxKRank, 2
			do q = -k, k, 1
				BRkq(k, q) = BRkq(k, q) * ( 1d0 - SternheimerShieldings(k) )
			end do
		end do
	end if
	
	! build the ligand-field matrix
	LFMatrix = cmplx(0.0d0, 0.0d0)
	do i = 1, 16
		do j = 1, 16
			! k = 2
			do q = -2, 2
				LFMatrix(i, j) = LFMatrix(i, j) + BRkq(2, q) * StevOp2Mat(q, i, j)
			end do
			! k = 4
			do q = -4, 4, 1
				LFMatrix(i, j) = LFMatrix(i, j) + BRkq(4, q) * StevOp4Mat(q, i, j)
			end do
			! k = 6
			do q = -6, 6, 1
				LFMatrix(i, j) = LFMatrix(i, j) + BRkq(6, q) * StevOp6Mat(q, i, j)
			end do
		end do
	end do

	EigenValues = 0d0
	
	! solve eigensystem
	call zheev(LAPACKJobType, LAPACKMatrixType, NArraySize, LFMatrix, NArraySizeLDA, EigenValues, LWorkArr, LWORKSize,  RWorkArr, lapackstatus)
	if ( lapackstatus .ne. 0 ) then
		write(*,*) "LAPACK Error!"
		stop
	end if

	! get the 15/2 contribution to the first KD
	val = max(dimag(LFMatrix(1,1))**2+dble(LFMatrix(1,1))**2, dimag(LFMatrix(16,1))**2+dble(LFMatrix(16,1))**2)

end subroutine

subroutine calc_cf_energy_pcm(lrota, lrotb, lrotg, val)
	use data_mo
	use data_grid
	use global_c
	implicit none
!
	double precision, intent(in) :: lrota, lrotb, lrotg
	double precision, intent(out) :: val
	integer :: i, k, q, j
	double precision :: srcX, srcY, srcZ, ttheta, tphi, ZdivR 
	double complex :: scval, SphericalHarmonicY
!
	if ( LAPACKInit .eqv. .FALSE. ) then
		call init_lapack
	end if
	
	if ( CFInit .eqv. .FALSE. ) then
		call init_cf_matrix_values
	end if
	
	ARkq = 0d0
	BRkq = 0d0
	BWkq = dcmplx(0.0d0, 0.0d0)
	
	! determine Waybourne LFPs
	do k = 2, MaxKRank, 2
		do q = 0, k, 1
			scval = dcmplx(0d0, 0d0)
			do i = 1, NCustomCharges
				if ( PointCharges(i)%dist .gt. 1.0d-6 ) then
			
					call rotate_euler_3d_inv(lrota, lrotb, lrotg, PointCharges(i)%x, PointCharges(i)%y, PointCharges(i)%z, srcX, srcY, srcZ)
					
					! avoid numerical issues
					ZdivR = srcZ / PointCharges(i)%dist
					if ( ZdivR .gt. 1d0 ) then
						ZdivR = 1d0
					else if ( ZdivR .lt. -1d0 ) then
						ZdivR = -1d0
					end if
					
					ttheta = dacos(ZdivR)
					tphi = datan2(srcY, srcX)
					
					scval = scval + dconjg(SphericalHarmonicY(k, q, ttheta, tphi)) * ( - PointCharges(i)%Charge ) / PointCharges(i)%dist**(k+1)
				end if
			end do
			
			if ( OSternheimer .eqv. .TRUE.  ) then
				BWkq(k, q) = scval * StevMultFac(k) * au2rcm * dsqrt((2d0*dble(k)+1d0)/FourPi)**2 * OldRadInts(k) * ( 1.0 - SternheimerShieldings(k) )
			else
				BWkq(k, q) = scval * StevMultFac(k) * au2rcm * dsqrt((2d0*dble(k)+1d0)/FourPi)**2 * OldRadInts(k) 
			end if
			
			if ( q .gt. 0 ) then
				BWkq(k,-q) = (-1)**q * dconjg(BWkq(k, q))
			end if
		end do
	end do
		
	! conversion to Stevens LFPs 
	BRkq = 0d0
	do k = 2, MaxKRank, 2
		BRkq(k, 0) = dble(BWkq(k, 0)) / StevLam(k, 0)
		do q = 0, k, 1
			BRkq(k, q) = dble(BWkq(k, q)) / StevLam(k, abs(q))
		end do
		do q = -k, -1, 1
			BRkq(k, q) = dimag(BWkq(k, abs(q))) / StevLam(k, abs(q))
		end do
		do q = -k, k
			ARkq(k, q) = BRkq(k, q) / StevMultFac(k)
		end do
	end do
	
	! build the ligand-field matrix
	LFMatrix = cmplx(0.0d0, 0.0d0)
	do i = 1, 16
		do j = 1, 16
			! k = 2
			do q = -2, 2
				LFMatrix(i, j) = LFMatrix(i, j) + BRkq(2, q) * StevOp2Mat(q, i, j)
			end do
			! k = 4
			do q = -4, 4, 1
				LFMatrix(i, j) = LFMatrix(i, j) + BRkq(4, q) * StevOp4Mat(q, i, j)
			end do
			! k = 6
			do q = -6, 6, 1
				LFMatrix(i, j) = LFMatrix(i, j) + BRkq(6, q) * StevOp6Mat(q, i, j)
			end do
		end do
	end do

	EigenValues = 0d0
	
	! solve eigensystem
	call zheev(LAPACKJobType, LAPACKMatrixType, NArraySize, LFMatrix, NArraySizeLDA, EigenValues, LWorkArr, LWORKSize,  RWorkArr, lapackstatus)
	if ( lapackstatus .ne. 0 ) then
		write(*,*) "LAPACK Error!"
		stop
	end if

	! get the 15/2 contribution to the first KD
	val = max(dimag(LFMatrix(1,1))**2+dble(LFMatrix(1,1))**2, dimag(LFMatrix(16,1))**2+dble(LFMatrix(16,1))**2)

end subroutine

subroutine calc_cf_energies_final(lrota, lrotb, lrotg)
	use data_mo
	use global_c
	use data_grid
	implicit none
!
	double precision, intent(inout) :: lrota, lrotb, lrotg
	integer :: k, l, q, i
	double precision :: zrot, zrotmin, zenergy, zenergymin, Radial_Expectation_Value_4f, dens
	double precision, dimension(8) :: KDEnergies, KDEnergiesMin, KDEnergiesMax
	integer :: KDSteps
!
	KDSteps = 0
	KDEnergies = 0d0
	WFDecomp = 0d0
	write(*,*)
	write(*,*) "> ENTERING LIGAND-FIELD PART ..."
	write(*,*)
	
	if ( LAPACKInit .eqv. .FALSE. ) then
		call init_lapack
	end if
	
	if ( CFInit .eqv. .FALSE. ) then
		call init_cf_matrix_values
	end if
	
	OSkipMom = .FALSE.
	
	if ( OVerbose .eqv. .TRUE. ) then
		write(*,'(5x,A)') "> PARAMETERS: "
		write(*,*)
		write(*,'(A35)') "multiplicative Stevens factors: "
		do l = 2, 6, 2
			write(*,'(A31,I1,A,F14.8)') "theta_", l," = " , StevMultFac(l)
		end do
		write(*,*)
		
		write(*,'(A35)') "Scaling factors: "
		do l = 2, 6, 2
			write(*,'(A31,I1,A,F14.8)') "f_", l," = " , EnergyScalingFactors(l)
		end do
		write(*,*)

		write(*,'(A35)') "radial expectation values: "
		do l = 2, 6, 2
			write(*,'(A30,I1,A,F14.8)') "<r^", l, "> = " , OldRadInts(l)
		end do
		write(*,*)
		write(*,*)
	end if
	
	write(*,'(5X,A)') "> ENERGIES OF KRAMERS DOUBLETS (IN 1/cm)"
	write(*,*)
	write(*,*) "	  ... utilizing Hartree-Potential to determine energies"
	write(*,*)
	if ( OLDAX .eqv. .TRUE. ) then
		write(*,*) "	  ... with LDA exchange for ligand charge density"
		write(*,*)
	end if
	write(*,*) "	  ... performing rotation about Euler angles:"
	write(*,*)
	write(*,'(A35,F14.2,A)') " alpha = ", lrota/pi*180d0, " degree"
	write(*,'(A35,F14.2,A)') "  beta = ", lrotb/pi*180d0, " degree"
	write(*,'(A35,F14.2,A)') " gamma = ", lrotg/pi*180d0, " degree"
	write(*,*)
	
	write(*,'(A35,F14.2,A)') " inner cutoff= ", LFTCutOff_I, " a0"
	write(*,'(A35,F14.2,A)') " outer cutoff= ", LFTCutOff_O, " a0"
	dens = 0d0
	do i = 1, NLGP
		if ( ( GridPoints(i)%dist .le. LFTCutOff_O ) .and. ( GridPoints(i)%dist .gt. LFTCutOff_I ) ) then
			dens = dens + NumberOf4fElectrons * GridPoints(i)%RadWFdV
		end if
	end do
	write(*,'(A35,F14.8)') " 4f charge density = ", dens
	write(*,*)
	
	zenergymin = 0d0
	KDEnergiesMin =  1d6
	KDEnergiesMax = -1d6
	
	call calc_cf_energy(lrota, lrotb, lrotg, zenergy)
		
	do l = 2, 8
		KDEnergies(l) = KDEnergies(l) + (EigenValues(l*2-1)-EigenValues(1))
	end do
	
	if ( EigenValues(1) .lt. KDEnergiesMin(1) ) then
		KDEnergiesMin(1) = EigenValues(1)
	end if
	if ( EigenValues(1) .gt. KDEnergiesMax(1) ) then
		KDEnergiesMax(1) = EigenValues(1)
	end if
	
	do l = 2, 8
		if ( (EigenValues(l*2-1)-EigenValues(1)) .lt. KDEnergiesMin(l) ) then
			KDEnergiesMin(l) = (EigenValues(l*2-1)-EigenValues(1))
		end if
		if ( (EigenValues(l*2-1)-EigenValues(1)) .gt. KDEnergiesMax(l) ) then
			KDEnergiesMax(l) = (EigenValues(l*2-1)-EigenValues(1))
		end if
	end do
		
	zenergymin = zenergy
	zrotmin = zrot
	do k = 1, 16
		do l = 0, 15
			WFDecomp(l, k) = 100d0 * ( dble(LFMatrix(1+l, k))**2 + dimag(LFMatrix(1+l, k))**2 )
		end do
	end do
		
	call write_rline(72)
	write(*,'(7X,8A9)') "KD1", "KD2", "KD3", "KD4", "KD5", "KD6", "KD7", "KD8"
	call write_rline(72)
	write(*,*)
	write(*,'(7X,8F9.2)') KDEnergies(1), KDEnergies(2), KDEnergies(3), KDEnergies(4), KDEnergies(5), KDEnergies(6), KDEnergies(7), KDEnergies(8)
	write(*,*)
	call write_rline(72)
	write(*,*)
	write(*,*)
	
	write(*,'(5X,A)') " > DETERMINED CRYSTAL-FIELD PARAMETERS (STEVENS NOTATION):"
	write(*,*)
	call write_rline(54)
	write(*,'(25X,2A6,A22,A20)') "k", "q", "A(k,q) / 1/cm", "B(k,q) / 1/cm"
	call write_rline(54)
		do k = 2, 6, 2
			do q = -k, k, 1
				write(*,'(25X,2I6,F22.8,E20.8)') k, q, BRkq(k, q) / StevMultFac(k) / Radial_Expectation_Value_4f(k), BRkq(k, q)
			end do
			if ( k .ne. 6 ) then
				write(*,*)
			end if
		end do
	call write_rline(54)
	write(*,*)
	
	if ( OVerbose .eqv. .TRUE. ) then
		write(*,*)
		write(*,'(5X,A)') "> MAXIMUM ESO RANK k AND RESULTING ENERGIES (IN 1/cm):"
		write(*,*)
		call write_rline(72)
		write(*,'(7X,A2,A7,7A9)') " k", "KD1", "KD2", "KD3", "KD4", "KD5", "KD6", "KD7", "KD8"
		call write_rline(72)
		write(*,*)
		do k = 2, 6, 2
			MaxKRank = k
			KDEnergies = 0d0
			call calc_cf_energy(lrota, lrotb, lrotg, zenergy)
			do l = 2, 8
				KDEnergies(l) = KDEnergies(l) + (EigenValues(l*2-1)-EigenValues(1))
			end do
			write(*,'(7X,I2,F7.2,7F9.2)') k, KDEnergies(1), KDEnergies(2), KDEnergies(3), KDEnergies(4), KDEnergies(5), KDEnergies(6), KDEnergies(7), KDEnergies(8)
		end do
		write(*,*)
		call write_rline(72)
		write(*,*)
		write(*,*)

		write(*,'(5X,A)') "> DECOMPOSITION OF KRAMERS DOUBLETS:"
		write(*,*)
		write(*,'(2X,A/)',advance='no') "---------------------------------------------------------------------------------------------------------------------"
		write(*,'(7X,16A7)') "+15/2", "+13/2", "+11/2", "+9/2", "+7/2", "+5/2", "+3/2", "+1/2", "-1/2", "-3/2", "-5/2", "-7/2", "-9/2", "-11/2", "-13/2", "-15/2"
		write(*,'(2X,A/)',advance='no') "---------------------------------------------------------------------------------------------------------------------"
		do k = 1, 8
			write(*,'(3X,A2,I1,X,16F7.1)') "KD", k, ( WFDecomp(l, k*2-1), l=0, 15)
			write(*,'(3X,A2,I1,X,16F7.1)') "KD", k, ( WFDecomp(l, k*2), l=0, 15)
			if ( k .ne. 8 ) then
				write(*,*)
			end if
		end do
		write(*,'(2X,A/)',advance='no') "---------------------------------------------------------------------------------------------------------------------"
		write(*,*)
		
	end if
	
end subroutine

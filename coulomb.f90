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

subroutine cpot_calc_energy(alpha, beta, gamma, value)
	use global_c
	use data_mo
	use data_grid
	implicit none
!
	double precision, intent(in) :: alpha, beta, gamma
	double precision, intent(out) :: value
	integer :: i, k
	double precision :: theta, RDyX, RDyY, RDyZ, GDyVal
	double complex :: SphericalHarmonicY
!
	Energy_SUM = 0d0

	!$OMP PARALLEL 
	!$OMP DO REDUCTION(+:Energy_SUM) PRIVATE(RDyX, RDyY, RDyZ, theta, GDyVal, k) SCHEDULE(dynamic)
	do i = 1, NLGP

		call rotate_euler_3d_inv(alpha, beta, gamma, GridPoints(i)%x, GridPoints(i)%y, GridPoints(i)%z, RDyX, RDyY, RDyZ)

		! avoid NaN due to numerical noise
		if ( RDyZ / GridPoints(i)%dist .gt. 1d0 ) then
			theta = dacos(  1d0 )
		else if ( RDyZ / GridPoints(i)%dist .lt. -1d0 ) then
			theta = dacos( -1d0 )
		else
			theta = dacos( RDyZ / GridPoints(i)%dist )
		end if

		GDyVal = coeff(mJ, 1) * real(SphericalHarmonicY(0, 0, theta, 0d0))
		do k = 2, 4
			GDyVal = GDyVal + coeff(mJ, k) * real(SphericalHarmonicY(2*(k-1), 0, theta, 0d0))
		end do

		Energy_SUM = Energy_SUM + GDyVal * FourPi * GridPoints(i)%RadWFdV * ( GridPoints(i)%ElecPotVal - GridPoints(i)%NucPotVal )
        
	end do
	!$OMP END DO NOWAIT
	!$OMP BARRIER
	!$OMP END PARALLEL

	value = Energy_SUM

end subroutine

!
! rhOver - a FORTRAN program to determine magnetic anisotropy and related 
! properties for dysprosium(III) single-ion magnets by semi-empirical approaches
! Copyright (C) 2014-2018 Michael Böhme <boehme.mic@gmail.com>
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
! * generates c_k parameters from Stevens factors A2, A4, A6
! *************************************************************************
subroutine generate_coefficients
	use global_c
	implicit none
!
	integer :: i, j
!
	do i = 1, 8
		coeff(i, 1) = 9d0 / dsqrt(4d0*pi)
		do j = 2, 4
			coeff(i, j) = StevensFactors_Dy(i, j-1)/dsqrt((4d0*pi)/(2d0*((j-1d0)*2d0)+1d0))
		end do
	end do
	
end subroutine

subroutine print_stevens_factors(verbose)
	use global_c
	implicit none
!
	logical, intent(in) :: verbose
	integer :: i, j
!
	if ( verbose .eqv. .TRUE. ) then
		write(*,*) "> STEVENS FACTORS FOR Dy(III):"
		write(*,*)
		call write_rline(65)
		write(*,'(14X,A12,2X,A16,A16,A16)') '|JM> state', 'A2', 'A4', 'A6'
		call write_rline(65)
		do i = 1, 8
			write(*,'(14X,A,I2,A,2X,3(F16.6))') '  | +/- ', (8-i)*2+1, '/2 > ', (StevensFactors_Dy(i, j), j=1, 3)
		end do
		call write_rline(65)
		write(*,*)
		
		write(*,*) "> GENERATED COEFFICIENTS FOR SPHERICAL HARMONICS:"
		write(*,*)
		call write_rline(75)
		write(*,'(8X,A,4A15)') '|JM> state', 'c(0)', 'c(2)', 'c(4)', 'c(6)'
		call write_rline(75)
		do i = 1, 8
			write(*,'(4X,A,I2,A,4(F15.6))') '  | +/- ',(8-i)*2+1, '/2 > ', ( coeff(i, j), j=1, 4)
		end do	
		call write_rline(75)
		write(*,*)
	end if
	
end subroutine

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

subroutine Radial_4f_WF(r, RadWF)
	use global_c
	implicit none
!
	double precision, intent(in) :: r
	double precision, intent(out) :: RadWF
!
	call Radial_4f_WF_ANORCC(r, RadWF)
	
end subroutine

subroutine Radial_4f_WF_ANORCC(r, RadWF)
	use global_c
	implicit none
!
	double precision, intent(in) :: r
	double precision, intent(out) :: RadWF
	double precision, parameter :: NormFac = 3.5338398682016225d0
!	

	if ( r .gt. 0d0 ) then
	
		RadWF = NormFac * &
			( 0.00282264d0 * dexp(-r**2 * 191.936524d0) + &
			  0.01054334d0 * dexp(-r**2 * 91.6221274d0) + &
			  0.04227740d0 * dexp(-r**2 * 43.8611428d0) + &
			  0.11522458d0 * dexp(-r**2 * 21.5452745d0) + &
			  0.24128770d0 * dexp(-r**2 * 10.5469112d0) + &
			  0.33641748d0 * dexp(-r**2 * 5.01110350d0) + &
			  0.33089588d0 * dexp(-r**2 * 2.30220189d0) + &
			  0.23406372d0 * dexp(-r**2 * 1.00615809d0) + &
			  0.10588540d0 * dexp(-r**2 * 0.40632415d0) + &
			  0.01641108d0 * dexp(-r**2 * 0.16252966d0) - &
			  0.00020473d0 * dexp(-r**2 * 0.06501186d0) )
		
	else
	
		RadWF = 0d0
		
	end if
	
end subroutine

function Radial_Expectation_Value_4f(k) result(res)
	use global_c
	implicit none
!
	integer, intent(in) :: k
	double precision :: res
	double precision, dimension(6) :: RadIntegral = &
		(/ &
			0.755478019325158d0, & ! < r^1 >
			0.725518019723808d0, & ! < r^2 >
			0.879250181934042d0, & ! < r^3 >
			1.322310787871935d0, & ! < r^4 >
			2.401518960767103d0, & ! < r^5 >
			5.101902744743668d0  & ! < r^6 >
		/)
!

	select case ( k ) 
	
		case ( 0 )
			res = 1d0
		case ( 1:6 )
		
			res = RadIntegral(k)
		
		case default
			write(*,*) "ERROR! <r^k> value not implemented!"
			res = 0d0
			stop
	
	end select
	
end function

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
! * 4f Radial Wave Function
! *************************************************************************
subroutine Radial_4f_WF(r, RadWF)
	use global_c
	implicit none
!
	double precision, intent(in) :: r
	double precision, intent(out) :: RadWF
!
	RadWF = 0d0
	
	select case ( iRadWF )
		case ( 1 )
			call Radial_4f_WF_ANORCC(r, RadWF)
		case ( 2 )
			call Radial_4f_WF_FreemanWatson(r, RadWF)
		case ( 3 )
			call Radial_4f_WF_MB(r, RadWF)
		case default
			write(*,*) "ERROR! Invalid radial wave function selected!"
	end select
	
end subroutine

! *************************************************************************
! * ANO-RCC based function for Dy(0)
! *   M. Böhme, W. Plass, 
! *   J. Comput. Chem., 2018, 39, 2697--2712
! *************************************************************************
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

! *************************************************************************
! * (renormalized) Radial 4f wave function for Dy(III):
! *   A.J. Freeman, R.E. Watson, 
! *   Phys. Rev., American Physical Society, 1962, 127, 2058--2075
! *************************************************************************
subroutine Radial_4f_WF_FreemanWatson(r, RadWF)
        use global_c
        implicit none
!
        double precision, intent(in) :: r
        double precision, intent(out) :: RadWF
        double precision, parameter :: FreemanWatsonNormFac = 1.000000063653808d0        
!        

        if ( r .gt. 0d0 ) then
        
                ! r**3 because P4f**2 = RadWF**2 * r**2
                RadWF = FreemanWatsonNormFac * &
                        ( 2480.4013d0 * dexp(-r*13.463d0) + &
                          448.83699d0 * dexp(-r* 7.529d0) + &
                          55.967002d0 * dexp(-r* 5.019d0) + &
                          2.3524738d0 * dexp(-r* 2.762d0) ) * r**3
                
        else
        
                RadWF = 0d0
                
        end if
        
end subroutine

! *************************************************************************
! * ANO-RCC / DKH2 / CASSCF - based function for Dy(III)
! * M. Böhme 2019 (unpublished)
! *************************************************************************
subroutine Radial_4f_WF_MB(r, RadWF)
        use global_c
        implicit none
!
        double precision, intent(in) :: r
        double precision, intent(out) :: RadWF
        double precision, parameter :: NormFac = 1.0000000d0        
!        

        if ( r .gt. 0d0 ) then
        
                ! r**3 because P4f**2 = RadWF**2 * r**2
                RadWF = NormFac * &
                        ( 2083.97277986782d0 * dexp(-r*16.74038946349620d0) + &
                          1095.79449625534d0 * dexp(-r* 9.62332852824475d0) + &
                          123.435880413911d0 * dexp(-r* 5.49609640503113d0) + &
                          3.91211713806515d0 * dexp(-r* 2.95944608005325d0) ) * r**3
                
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
	double precision, dimension(6) :: ANORCCIntegrals = &
		(/ &
			0.7530082024615980d0, & ! < r^1 >
			0.7929514365528688d0, & ! < r^2 >
			1.0930798149695940d0, & ! < r^3 >
			1.8767954644411686d0, & ! < r^4 >
			3.8616607663314610d0, & ! < r^5 >
			9.2387984090389160d0  & ! < r^6 >
		/)
	double precision, dimension(6) :: FWIntegrals = &
		(/ &
			0.755478019325158d0, & ! < r^1 >
			0.725518019723808d0, & ! < r^2 >
			0.879250181934042d0, & ! < r^3 >
			1.322310787871935d0, & ! < r^4 >
			2.401518960767103d0, & ! < r^5 >
			5.101902744743668d0  & ! < r^6 >
		/)
	double precision, dimension(6) :: CASSCF_4f_MB_Integrals = &
		(/ &
			0.778113353308d0, & ! < r^1 >
			0.768525333170d0, & ! < r^2 >
			0.950328550536d0, & ! < r^3 >
			1.445494084350d0, & ! < r^4 >
			2.660968434839d0, & ! < r^5 >
			5.920293113883d0  & ! < r^6 >
		/)
!
	select case ( k ) 
	
		case ( 0 )
			res = 1d0
		case ( 1:6 )
		
			select case ( iRadWF )
				case ( 1 )
					res = ANORCCIntegrals(k)
				case ( 2 )
					res = FWIntegrals(k)
				case ( 3 )
					res = CASSCF_4f_MB_Integrals(k)
				case default
					write(*,*) "ERROR! Invalid radial wave function selected!"
			end select
		
		case default
			write(*,*) "ERROR! <r^k> value not implemented!"
			res = 0d0
			stop
	
	end select
	
end function

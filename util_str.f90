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
! * prints the program header 
! *************************************************************************
subroutine print_header
	implicit none
!
	write(*,*)
	write(*,*) "*****************************************************************************"
	write(*,*) "*                                                                           *"
	write(*,*) "*                                 #######                                   *"
	write(*,*) "*                   #####  #    # #     # #    # ###### #####               *"  
	write(*,*) "*                   #    # #    # #     # #    # #      #    #              *" 
	write(*,*) "*                   #    # ###### #     # #    # #####  #    #              *"
	write(*,*) "*                   #####  #    # #     # #    # #      #####               *"
	write(*,*) "*                   #   #  #    # #     #  #  #  #      #   #               *"
	write(*,*) "*                   #    # #    # #######   ##   ###### #    #       v1.00  *"
	write(*,*) "*                                                                           *"
	write(*,*) "*****************************************************************************"
	write(*,*)
	write(*,*) "                                   - * -                                     "
	write(*,*)
	write(*,*) "          written 2014-2019 by Michael Böhme - boehme.mic@gmail.com          "
#ifdef __GFORTRAN__
	write(*,*) "                     build date: ", __DATE__, " ", __TIME__
#endif
#ifndef RELEASE
	write(*,*)
	write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	write(*,*) "!                                ATTENTION!                                 !"
	write(*,*) "!      THIS IS AN UNTESTED DEVELOPMENT VERSION! USE AT YOUR OWN RISK!       !"
	write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif
	write(*,*)
	write(*,*) "Please cite:"
	write(*,*) "  M. Böhme, W. Plass, J. Comput. Chem. 2018 (39), 2697--2712."
	write(*,*)        

end subroutine

! *************************************************************************
! * gets the chemical symbol for a given charge
! *************************************************************************
subroutine get_Element_by_Charge(charge, ostr)
	implicit none
!
	integer, intent(in) :: charge
	character (len=*), intent(out) :: ostr
	integer :: pos1, pos2, i
	character (len=1000) :: Elements = &
		"H:He:Li:Be:B:C:N:O:F:Ne:Na:Mg:Al:Si:P:S:Cl:Ar:" // &
		"K:Ca:Sc:Ti:V:Cr:Mn:Fe:Co:Ni:Cu:Zn:Ga:Ge:As:Se:Br:Kr:" // &
		"Rb:Sr:Y:Zr:Nb:Mo:Tc:Ru:Rh:Pd:Ag:Cd:In:Sn:Sb:Te:I:Xe:" // &
		"Cs:Ba:La:Ce:Pr:Nd:Pm:Sm:Eu:Gd:Tb:Dy:Ho:Er:Tm:Yb:Lu:"  // &
		"Hf:Ta:W:Re:Os:Ir:Pt:Au:Hg:Tl:Pb:Bi:Po:At:Rn:Fr:Ra:"   // &
		"Ac:Th:Pa:U:Np:Pu:Am:Cm:Bk:Cf:Es:Fm:Md:No:Lr:"         // &
		"Rf:Db:Sg:Bh:Hs:Mt:Ds:Rg:Cn:Nh:Fl:Mc:Lv:Ts:Og:XX"
!
	i = 0
	pos1 = 1
	ostr = 'Xx'
	
	if ( charge .eq. 0 ) then
		return
	end if
	
	do i = 1, 119
		pos2 = index(Elements(pos1:), ":")
		if ( pos2 .ne. 0 ) then
			ostr = trim(adjustl(Elements(pos1:pos1+pos2-2)))
			if ( i .eq. charge ) then
				exit
			end if
			pos1 = pos1 + pos2 
		else
			exit
		end if
	end do
	
end subroutine

! *************************************************************************
! * converts all characters of a string to upper case
! *************************************************************************
function to_upper(input) result(output)
	implicit none
!
	character(len=*), intent(in) :: input
	character(len=len(input)) :: output
	integer :: i
!
	do i = 1, len(input)
		select case ( iachar(input(i:i)) )
			case ( 97:122 )
				output(i:i) = achar(iachar(input(i:i))-32)
			case default
				output(i:i) = input(i:i)
		end select
	end do
	
end function

! *************************************************************************
! * converts the first character of a string to upper case
! *************************************************************************
function to_upper_first(input) result(output)
	implicit none
!
	character(len=*), intent(in) :: input
	character(len=len(input)) :: output
	integer :: i
!
	do i = 1, len(input)
		select case ( iachar(input(i:i)) )
			case ( 97:122 )
				if ( i .eq. 1 ) then
					output(i:i) = achar(iachar(input(i:i))-32)
				else
					output(i:i) = input(i:i)
				end if
			case default
				output(i:i) = input(i:i)
		end select
	end do
	
end function

! *************************************************************************
! * gets the shell type
! *************************************************************************
function get_shell_type(id) result(chr)
	implicit none
!
	integer, intent(in) :: id
	character(len=1) :: chr
!
	select case(id)
		case ( 1 )
			chr = "s"
		case ( 2 )
			chr = "p"
		case ( 3 )
			chr = "d"
		case ( 4 )
			chr = "f"
		case ( 5 )
			chr = "g"
		case ( 6 )
			chr = "h"
		case ( 7 )
			chr = "i"
		case default
			chr = "?"
	end select
end function

! *************************************************************************
! * writes a bar (right justified)
! *************************************************************************
subroutine write_rline(lchars)
	implicit none
!
	integer, intent(in) :: lchars
	integer :: i
	integer, parameter :: linelength = 80
!
	do i = linelength - lchars - 1, 1, -1
		write(*,'(X)',advance='no') 
	end do
	
	do i = 1, lchars
		write(*,'(A)',advance='no') '-'
	end do
	
	write(*,'(/)',advance='no')

end subroutine

! *************************************************************************
! * writes a bar (left justified)
! *************************************************************************
subroutine write_line(lchars)
	implicit none
!
	integer, intent(in) :: lchars
	integer :: i
	integer, parameter :: linelength = 80
!
	do i = 1, lchars
		write(*,'(A)',advance='no') '-'
	end do
	
	write(*,'(/)',advance='no')

end subroutine

! *************************************************************************
! * writes a boolean value
! *************************************************************************
subroutine write_option_bool(bool)
	implicit none
!
	logical, intent(in) :: bool
!
	if ( bool .eqv. .TRUE. ) then
		write(*,'(A14)') "enabled"
	else
		write(*,'(A14)') "disabled"
	end if

end subroutine

! *************************************************************************
! * writes 'DONE!'
! *************************************************************************
subroutine write_done
	implicit none
!
	write(*,'(70X,A)') "... DONE!"
	write(*,*)
end subroutine


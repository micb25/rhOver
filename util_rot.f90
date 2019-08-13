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
! * rotation about Euler angles alpha, beta, gamma 
! *************************************************************************
subroutine rotate_euler_3d(alpha, beta, gamma, srcX, srcY, srcZ, destX, destY, destZ)
	implicit none
!
	double precision, intent(in) :: alpha, beta, gamma, srcX, srcY, srcZ
	double precision, intent(out) :: destX, destY, destZ
!
	destX = srcZ*dsin(alpha)*dsin(beta) + srcY*(-(dcos(beta)*dcos(gamma)*dsin(alpha)) - dcos(alpha)*dsin(gamma)) + &
		srcX*(dcos(alpha)*dcos(gamma) - dcos(beta)*dsin(alpha)*dsin(gamma))
	destY = -(srcZ*dcos(alpha)*dsin(beta)) + srcX*(dcos(gamma)*dsin(alpha) + dcos(alpha)*dcos(beta)*dsin(gamma)) + &
		srcY*(dcos(alpha)*dcos(beta)*dcos(gamma) - dsin(alpha)*dsin(gamma))
	destZ = srcZ*dcos(beta) + srcY*dcos(gamma)*dsin(beta) + srcX*dsin(beta)*dsin(gamma)

end subroutine

! *************************************************************************
! * inverse rotation about Euler angles alpha, beta, gamma 
! *************************************************************************
subroutine rotate_euler_3d_inv(alpha, beta, gamma, srcX, srcY, srcZ, destX, destY, destZ)
	implicit none
!
	double precision, intent(in) :: alpha, beta, gamma, srcX, srcY, srcZ
	double precision, intent(out) :: destX, destY, destZ
!
	destX = srcZ*dsin(beta)*dsin(gamma) + srcY*(dcos(gamma)*dsin(alpha) + dcos(alpha)*dcos(beta)*dsin(gamma)) + &
		srcX*(dcos(alpha)*dcos(gamma) - dcos(beta)*dsin(alpha)*dsin(gamma))
	destY = srcZ*dcos(gamma)*dsin(beta) + srcX*(-(dcos(beta)*dcos(gamma)*dsin(alpha)) - dcos(alpha)*dsin(gamma)) + &
		srcY*(dcos(alpha)*dcos(beta)*dcos(gamma) - dsin(alpha)*dsin(gamma))
	destZ = srcZ*dcos(beta) - srcY*dcos(alpha)*dsin(beta) + srcX*dsin(alpha)*dsin(beta)

end subroutine

! *************************************************************************
! * calculates the angle between two vectors
! *************************************************************************
function calc_vec_angle(vax, vay, vaz, vbx, vby, vbz) result(res)
	use global_c
	implicit none
!
	double precision, intent(in) :: vax, vay, vaz, vbx, vby, vbz
	double precision :: res
!
	res = 180d0/pi * dacos( (vax*vbx+vay*vby+vaz*vbz) / &
		( dsqrt(vax**2 + vay**2 + vaz**2) * dsqrt(vbx**2 + vby**2 + vbz**2)) ) 
	if ( res .gt. 90d0 ) then
		res= 180d0 - res
	end if 
	
end function

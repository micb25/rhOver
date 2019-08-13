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
! * spherical harmonics up to l = 15 with m = -l, ..., +l
! *************************************************************************
function SphericalHarmonicY(l, m, theta, phi) result(res)
	use global_c
	implicit none
!	
	integer, intent(in) :: l, m
	double precision, intent(in) :: theta, phi
	double complex :: res
	double precision :: costheta, sintheta
	double complex :: prefac
!	
	if ( l .ne. 0 ) then
		costheta = dcos(theta)
	end if
	
	if ( m .ne. 0 ) then
	
		sintheta = dsin(theta)
		prefac = cdexp( m * cmplx(0.0, 1.0) * phi ) 
		if ( m .gt. 0 ) then
			prefac = prefac * ( - 1 )**m
		end if
		
	end if

	select case ( l )
	
		case ( 0 )
			res = 0.5d0 / dsqrt( pi )
			return
			
		case ( 1 )
			
			select case ( m )
				case ( -1, 1)
					res = 0.5d0 * prefac * dsqrt( 3.0d0 / ( 2.0d0 * pi) ) * sintheta
					return
				case ( 0 )
					res = 0.5d0 * dsqrt( 3.0d0 / pi ) * costheta
					return
			end select
			
		case ( 2 )
			
			select case ( m )
				case ( -2, 2 )
					res = 0.25d0 * prefac * dsqrt( 15.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**2
					return
				case ( -1, 1 )
					res = 0.50d0 * prefac * dsqrt( 15.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta * costheta
					return
				case ( 0 )
					res = 0.25d0 * dsqrt( 5.0d0 / pi ) * ( 3.d0 * costheta**2 - 1.0d0 )
					return
			end select
			
		case ( 3 )
		
			select case ( m )
				case ( -3, 3 )
					res = 0.125d0 * prefac * dsqrt( 35.0d0 / pi ) * &
						sintheta**3
					return
				case ( -2, 2 )
					res = 0.250d0 * prefac * dsqrt( 105.d0 / ( 2.0d0 * pi ) ) * &
						sintheta**2 * costheta
					return
				case ( -1, 1 )
					res = 0.125d0 * prefac * dsqrt( 21.d0 / pi ) * &
						sintheta * ( 5.0d0 * costheta**2 - 1.0d0 )
					return
				case ( 0 )
					res = 0.25d0 * dsqrt( 7.0d0 / pi ) * ( 5.d0 * costheta**3 - 3.0d0 * costheta )
					return
			end select
			
		case ( 4 )
		
			select case ( m )
				case ( -4, 4 )
					res = 3.0d0 / 16.0d0 * prefac * dsqrt( 35.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**4
					return
				case ( -3, 3 )
					res = 3.0d0 / 8.0d0 * prefac * dsqrt( 35.0d0 / ( pi ) ) * &
						sintheta**3 * costheta
					return
				case ( -2, 2 )
					res = 3.0d0 / 8.0d0 * prefac * dsqrt( 5.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**2 * ( 7.0d0 * costheta**2 - 1 )
					return
				case ( -1, 1 )
					res = 3.0d0 / 8.0d0 * prefac * dsqrt( 5.0d0 / ( pi ) ) * &
						costheta * sintheta * ( 7.0d0  * costheta**2  - 3.0d0 )
					return
				case ( 0 )
					res = 3.0d0 * ( 35d0*costheta**4 - 30d0*costheta**2 + 3d0) / ( 16d0 * dsqrt( pi ) )
					return
			end select
			
		case ( 5 )
			
			select case ( m )
				case ( -5, 5 )
					res = 3.0d0 / 32.d0 * prefac * dsqrt( 77.0d0 / ( pi ) ) * &
						sintheta**5
					return
				case ( -4, 4 )
					res = 3.0d0 / 16.d0 * prefac * dsqrt( 385.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**4 * costheta
					return
				case ( -3, 3 )
					res = 1.0d0 / 32.d0 * prefac * dsqrt( 385.0d0 / ( pi ) ) * &
						sintheta**3 * ( 9.0d0 * costheta**2 - 1.0d0 )
					return
				case ( -2, 2 )
					res = 1.0d0 / 8.d0 * prefac * dsqrt( 1155.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**2 * costheta * ( 3.0d0*costheta**2 - 1.0d0 )
					return
				case ( -1, 1 )
					res = 1.0d0 / 16.d0 * prefac * dsqrt( 165.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta * ( 21.0d0 * costheta**4 - 14.0d0 * costheta**2 + 1.0d0 )
					return
				case ( 0 )
					res = 1.0d0 / 16.d0 * dsqrt( 11.0d0 / pi ) * &
						( 63.0d0*costheta**5 - 70.0d0*costheta**3 + 15.0d0*costheta )
					return
			end select
			
		case ( 6 )
		
			select case ( m )
				case ( -6, 6 )
					res = prefac / 64.0d0 * dsqrt( 3003.0d0 / pi ) * &
						sintheta**6
					return
				case ( -5, 5 )
					res = 3.0d0 / 32.0d0 * prefac * dsqrt( 1001.0d0 / pi ) * &
						sintheta**5 * costheta
					return
				case ( -4, 4 )
					res = 3.0d0 / 32.0d0 * prefac * dsqrt( 91.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**4 * ( 11.0d0 * costheta**2 - 1.d0 )
					return
				case ( -3, 3 )
					res = 1.0d0 / 32.0d0 * prefac * dsqrt( 1365.0d0 / pi ) * &
						sintheta**3 * costheta * ( 11.0d0 * costheta**2 - 3.d0 )
					return
				case ( -2, 2 )
					res = 1.0d0 / 64.0d0 * prefac * dsqrt( 1365.0d0 / pi ) * &
						sintheta**2 * ( 33.0d0 * costheta**4 - 18.d0 * costheta**2 + 1.d0 )
					return
				case ( -1, 1 )
					res = 1.0d0 / 16.0d0 * prefac * dsqrt( 273.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta * costheta * ( 33.0d0 * costheta**4 - 30.d0 * costheta**2 + 5.d0 )
					return
				case ( 0 )
					res = dsqrt(13.d0 / pi) / 32.d0 * ( &
						+ 231.d0 * costheta**6 - 315.0d0 * costheta**4 + 105.d0 * costheta**2 - 5.0d0 )
					return
			end select
			
		case ( 7 )
		
			select case ( m )
				case ( -7, 7 )
					res = 3.0d0 / 64.0d0 * prefac * dsqrt( 715d0 / ( 2d0*pi ) ) * &
						sintheta**7
					return
				case ( -6, 6 )
					res = 3.0d0 / 64.0d0 * prefac * dsqrt( 5005d0 / ( pi ) ) * &
						costheta * sintheta**6
					return
				case ( -5, 5 )
					res = 3.0d0 / 64.0d0 * prefac * dsqrt( 385d0 / ( 2d0*pi ) ) * &
						sintheta**5 * (13d0*costheta**2 - 1d0)			
					return
				case ( -4, 4 )
					res = 3.0d0 / 32.0d0 * prefac * dsqrt( 385d0 / ( 2d0*pi ) ) * &
						costheta * sintheta**4 * (13d0*costheta**2 - 3d0)
					return			
				case ( -3, 3 )
					res = 3.0d0 / 64.0d0 * prefac * dsqrt( 35.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta**3 * (143d0*costheta**4 - 66d0*costheta**2 + 3d0)
					return
				case ( -2, 2 )					
					res = 3.0d0 / 64.0d0 * prefac * dsqrt( 35.0d0 / ( pi ) ) * &
						costheta * sintheta**2 * (143d0*costheta**4 - 110d0*costheta**2 + 15d0)			
					return
				case ( -1, 1 )					
					res = 1.0d0 / 64.0d0 * prefac * dsqrt( 105.0d0 / ( 2.0d0 * pi ) ) * &
						sintheta * (429d0*costheta**6 - 495d0*costheta**4 + 135d0*costheta**2 - 5d0)
					return			
				case ( 0 )
					res = 1d0 / 32d0 * dsqrt(15.d0 / pi) * ( &
						429d0*costheta**7 - 693d0*costheta**5 + 315d0*costheta**3 - 35d0*costheta)
					return
			end select
			
		case ( 8 )
		
			select case ( m )
				case ( -8, 8 )
					res = 3d0 / 256d0 * prefac * dsqrt( 12155d0 / ( 2d0*pi ) ) * &
						sintheta**8 
					return			
				case ( -7, 7 )
					res = 3d0 / 64d0 * prefac * dsqrt( 12155d0 / ( 2d0*pi ) ) * &
						costheta * sintheta**7 
					return	
				case ( -6, 6 )
					res = 1d0 / 128d0 * prefac * dsqrt( 7293d0 / ( pi ) ) * &
						sintheta**6 * ( 15d0*costheta**2 - 1d0 )			
					return
				case ( -5, 5 )
					res = 3d0 / 64d0 * prefac * dsqrt( 17017d0 / ( 2d0*pi ) ) * &
						costheta * sintheta**5 * ( 5d0*costheta**2 - 1d0 )
					return
				case ( -4, 4 )
					res = 3d0 / 128d0 * prefac * dsqrt( 1309d0 / ( 2d0*pi ) ) * &
						sintheta**4 * ( 65d0*costheta**4 - 26d0*costheta**2 + 1d0 )			
					return
				case ( -3, 3 )
					res = 1d0 / 64d0 * prefac * dsqrt( 19635d0 / ( 2d0*pi ) ) * &
						costheta * sintheta**3 * ( 39d0*costheta**4 - 26d0*costheta**2 + 3d0 )
					return
				case ( -2, 2 )
					res = 3d0 / 128d0 * prefac * dsqrt( 595d0 / ( pi ) ) * &
						sintheta**2 * ( 143d0*costheta**6 - 143d0*costheta**4 + 33d0*costheta**2 - 1d0 )
					return
				case ( -1, 1 )
					res = 3d0 / 64d0 * prefac * dsqrt( 17d0 / ( 2d0*pi ) ) * &
						sintheta * costheta * ( 715d0*costheta**6 - 1001d0*costheta**4 + 385d0*costheta**2 - 35d0 )
					return
				case ( 0 )
					res = 1d0/256d0 * dsqrt(17d0/pi) * ( &
						6435d0*costheta**8 - 12012d0*costheta**6 + 6930d0*costheta**4 - 1260d0*costheta**2 + 35d0)
					return
			end select
			
		case ( 9 )
		
			select case ( m )
				case ( -9, 9 )
					res = 1d0 / 512d0 * prefac * dsqrt( 230945d0 / ( pi ) ) * &
						sintheta**9
					return
				case ( -8, 8 )
					res = 3d0 / 256d0 * prefac * dsqrt( 230945d0 / ( 2d0*pi ) ) * &
						sintheta**8 * costheta
					return
				case ( -7, 7 )
					res = 3d0 / 512d0 * prefac * dsqrt( 13585d0 / ( pi ) ) * &
						sintheta**7 * ( 17d0*costheta**2 - 1d0 )			
					return
				case ( -6, 6 )
					res = 1d0 / 128d0 * prefac * dsqrt( 40755d0 / ( pi ) ) * &
						sintheta**6 * costheta * ( 17d0*costheta**2 - 3d0 )			
					return
				case ( -5, 5 )
					res = 3d0 / 256d0 * prefac * dsqrt( 2717d0 / ( pi ) ) * &
						sintheta**5 * ( 85d0*costheta**4 - 30d0*costheta**2 + 1d0 )
					return
				case ( -4, 4 )
					res = 3d0 / 128d0 * prefac * dsqrt( 95095d0 / ( 2d0*pi ) ) * &
						sintheta**4 *costheta * ( 17d0*costheta**4 - 10d0*costheta**2 + 1d0 )			
					return
				case ( -3, 3 )
					res = 1d0 / 256d0 * prefac * dsqrt( 21945d0 / ( pi ) ) * &
						sintheta**3 * ( 221d0*costheta**6 - 195d0*costheta**4 + 39d0*costheta**2 - 1d0 )
					return
				case ( -2, 2 )
					res = 3d0 / 128d0 * prefac * dsqrt( 1045d0 / ( pi ) ) * &
						sintheta**2 * costheta * ( 221d0*costheta**6 - 273d0*costheta**4 + 91d0*costheta**2 - 7d0 )
					return
				case ( -1, 1 )
					res = 3d0 / 256d0 * prefac * dsqrt( 95d0 / ( 2d0*pi ) ) * &
						sintheta * ( 2431d0*costheta**8 - 4004d0*costheta**6 + 2002d0*costheta**4 - 308d0*costheta**2 + 7d0 )
					return
				case ( 0 )
					res = 1d0/256d0 * dsqrt(19d0/pi) * ( &
						12155d0*costheta**9 - 25740d0*costheta**7 + 18018d0*costheta**5 - 4620d0*costheta**3 + 315d0*costheta)		
					return
			end select
			
		case ( 10 )
		
			select case ( m )
				case ( -10, 10 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 969969d0 / ( pi ) ) * &
						sintheta**10
					return			
				case ( -9, 9 )
					res = 1d0 / 512d0 * prefac * dsqrt( 4849845d0 / ( pi ) ) * &
						sintheta**9 * costheta
					return
				case ( -8, 8 )
					res = 1d0 / 512d0 * prefac * dsqrt( 255255d0 / ( 2d0*pi ) ) * &
						sintheta**8 * ( 19d0*costheta**2 - 1d0 )			
					return
				case ( -7, 7 )
					res = 3d0 / 512d0 * prefac * dsqrt( 85085d0 / ( pi ) ) * &
						sintheta**7 * costheta * ( 19d0*costheta**2 - 3d0 )			
					return
				case ( -6, 6 )
					res = 3d0 / 1024d0 * prefac * dsqrt( 5005d0 / ( pi ) ) * &
						sintheta**6 * ( 323d0*costheta**4 - 102d0*costheta**2 + 3d0 )			
					return
				case ( -5, 5 )
					res = 3d0 / 256d0 * prefac * dsqrt( 1001d0 / ( pi ) ) * &
						sintheta**5 * costheta * ( 323d0*costheta**4 - 170d0*costheta**2 + 15d0 )
					return
				case ( -4, 4 )
					res = 3d0 / 256d0 * prefac * dsqrt( 5005d0 / ( 2d0*pi ) ) * &
						sintheta**4 * ( 323d0*costheta**6 - 255d0*costheta**4 + 45d0*costheta**2 - 1d0 )
					return
				case ( -3, 3 )
					res = 3d0 / 256d0 * prefac * dsqrt( 5005d0 / ( pi ) ) * &
						sintheta**3 * costheta * ( 323d0*costheta**6 - 357d0*costheta**4 + 105d0*costheta**2 - 7d0 )
					return
				case ( -2, 2 )
					res = 3d0 / 512d0 * prefac * dsqrt( 385d0 / ( 2d0*pi ) ) * &
						sintheta**2 * ( 4199d0*costheta**8 - 6188d0*costheta**6 + 2730d0*costheta**4 - 364d0*costheta**2 + 7d0 )
					return
				case ( -1, 1 )
					res = 1d0 / 256d0 * prefac * dsqrt( 1155d0 / ( 2d0*pi ) ) * &
						sintheta * costheta * ( 4199d0*costheta**8 - 7956d0*costheta**6 + 4914d0*costheta**4 - 1092d0*costheta**2 + 63d0 )
					return
				case ( 0 )
					res = 1d0/512d0 * dsqrt(21d0/pi) * ( &
						46189d0*costheta**10 - 109395d0*costheta**8 + 90090d0*costheta**6 - 30030d0*costheta**4 + 3465d0*costheta**2 - 63d0)
					return
			end select
			
		case ( 11 )
		
			select case ( m )	
				case ( -11, 11 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 2028117d0 / ( 2d0*pi ) ) * &
						sintheta**11
					return
				case ( -10, 10 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 22309287d0 / ( pi ) ) * &
						sintheta**10 * costheta
					return
				case ( -9, 9 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 1062347d0 / ( 2d0*pi ) ) * &
						sintheta**9 * ( 21d0*costheta**2 - 1d0 )
					return
				case ( -8, 8 )
					res = 1d0 / 512d0 * prefac * dsqrt( 15935205d0 / ( 2d0*pi ) ) * &
						sintheta**8 * costheta * ( 7d0*costheta**2 - 1d0 )
					return
				case ( -7, 7 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 838695d0 / ( 2d0*pi ) ) * &
						sintheta**7 * ( 133d0*costheta**4 - 38d0*costheta**2 + 1d0 )
					return
				case ( -6, 6 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 167739d0 / ( pi ) ) * &
						sintheta**6 * costheta * ( 399d0*costheta**4 - 190d0*costheta**2 + 15d0 )			
					return
				case ( -5, 5 )
					res = 3d0 / 1024d0 * prefac * dsqrt( 3289d0 / ( 2d0*pi ) ) * &
						sintheta**5 * ( 2261d0*costheta**6 - 1615d0*costheta**4 + 255d0*costheta**2 - 5d0 )
					return
				case ( -4, 4 )
					res = 3d0 / 256d0 * prefac * dsqrt( 23023d0 / ( 2d0*pi ) ) * &
						sintheta**4 * costheta * ( 323d0*costheta**6 - 323d0*costheta**4 + 85d0*costheta**2 - 5d0 )
					return
				case ( -3, 3 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 345345d0 / ( pi ) ) * &
						sintheta**3 * ( 969d0*costheta**8 - 1292d0*costheta**6 + 510d0*costheta**4 - 60d0*costheta**2 + 1d0 )
					return
				case ( -2, 2 )
					res = 1d0 / 512d0 * prefac * dsqrt( 49335d0 / ( 2d0*pi ) ) * &
						sintheta**2 * costheta * ( 2261d0*costheta**8 - 3876d0*costheta**6 + 2142d0*costheta**4 - 420d0*costheta**2 + 21d0 )
					return
				case ( -1, 1 )
					res = 1d0 / 1024d0 * prefac * dsqrt( 759d0 / ( pi ) ) * &
						sintheta * ( 29393d0*costheta**10 - 62985d0*costheta**8 + 46410d0*costheta**6 - 13650d0*costheta**4 + 1365d0*costheta**2 - 21d0 )
					return
				case ( 0 )
					res = 1d0/512d0 * dsqrt(23d0/pi) * ( &
						88179d0*costheta**11 - 230945d0*costheta**9 + 218790d0*costheta**7 - 90090d0*costheta**5 + 15015d0*costheta**3 - 693d0*costheta)
					return
			end select
			
		case ( 12 )
		
			select case ( m )
				case ( -12, 12 )
					res = 5d0 / 4096d0 * prefac * dsqrt( 676039d0 / ( pi ) ) * &
						sintheta**12
					return
				case ( -11, 11 )
					res = 5d0 / 1024d0 * prefac * dsqrt( 2028117d0 / ( 2d0*pi ) ) * &
						sintheta**11 * costheta
					return
				case ( -10, 10 )
					res = 5d0 / 2048d0 * prefac * dsqrt( 88179d0 / ( pi ) ) * &
						sintheta**10 * ( 23d0*costheta**2 - 1d0)
					return
				case ( -9, 9 )
					res = 5d0 / 1024d0 * prefac * dsqrt( 323323d0 / ( 2d0*pi ) ) * &
						sintheta**9 * costheta * ( 23d0*costheta**2 - 3d0)
					return		
				case ( -8, 8 )
					res = 5d0 / 2048d0 * prefac * dsqrt( 138567d0 / ( 2d0*pi ) ) * &
						sintheta**8 * ( 161d0*costheta**4 - 42d0*costheta**2 + 1d0)			
					return
				case ( -7, 7 )
					res = 5d0 / 1024d0 * prefac * dsqrt( 138567d0 / ( 2d0*pi ) ) * &
						sintheta**7 * costheta * ( 161d0*costheta**4 - 70d0*costheta**2 + 5d0)
					return
				case ( -6, 6 )
					res = 5d0 / 2048d0 * prefac * dsqrt( 2431d0 / ( pi ) ) * &
						sintheta**6 * ( 3059d0*costheta**6 - 1995d0*costheta**4 + 285d0*costheta**2 - 5d0)
					return
				case ( -5, 5 )
					res = 15d0 / 1024d0 * prefac * dsqrt( 17017d0 / ( 2d0*pi ) ) * &
						sintheta**5 * costheta * ( 437d0*costheta**6 - 399d0*costheta**4 + 95d0*costheta**2 - 5d0)
					return
				case ( -4, 4 )
					res = 15d0 / 4096d0 * prefac * dsqrt( 1001d0 / ( pi ) ) * &
						sintheta**4 * ( 7429d0*costheta**8 - 9044d0*costheta**6 + 3230d0*costheta**4 - 340d0*costheta**2 + 5d0)
					return
				case ( -3, 3 )
					res = 5d0 / 1024d0 * prefac * dsqrt( 1001d0 / ( pi ) ) * &
						sintheta**3 * costheta * ( 7429d0*costheta**8 - 11628d0*costheta**6 + 5814d0*costheta**4 - 1020d0*costheta**2 + 45d0)
					return
				case ( -2, 2 )
					res = 5d0 / 1024d0 * prefac * dsqrt( 3003d0 / ( 2d0*pi ) ) * &
						sintheta**2 * ( 7429d0*costheta**10 - 14535d0*costheta**8 + 9690d0*costheta**6 - 2550d0*costheta**4 + 225d0*costheta**2 - 3d0)
					return
				case ( -1, 1 )
					res = 5d0 / 1024d0 * prefac * dsqrt( 39d0 / ( pi ) ) * &
						sintheta * costheta * ( 52003d0*costheta**10 - 124355d0*costheta**8 + 106590d0*costheta**6 - 39270d0*costheta**4 + 5775d0*costheta**2 - 231d0)
					return
				case ( 0 )
					res = 5d0/2048d0 * dsqrt(1d0/pi) * ( &
						676039d0*costheta**12 - 1939938d0*costheta**10 + 2078505d0*costheta**8 - 1021020d0*costheta**6 + 225225d0*costheta**4 - 18018d0*costheta**2 + 231d0)
					return
			end select
			
		case ( 13 )
		
			select case ( m )
				case ( -13, 13 )
					res = 15d0 / 4096d0 * prefac * dsqrt( 156009d0 / ( 2d0*pi ) ) * &
						sintheta**13
					return
				case ( -12, 12 )
					res = 15d0 / 4096d0 * prefac * dsqrt( 2028117d0 / ( pi ) ) * &
						sintheta**12 * costheta
					return
				case ( -11, 11 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 2028117d0 / ( 2d0*pi ) ) * &
						sintheta**11 * ( 25d0*costheta**2 - 1d0)
					return
				case ( -10, 10 )
					res = 3d0 / 2048d0 * prefac * dsqrt( 2028117d0 / ( pi ) ) * &
						sintheta**10 * costheta * ( 25d0*costheta**2 - 3d0)
					return
				case ( -9, 9 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 88179d0 / ( pi ) ) * &
						sintheta**9 * ( 575d0*costheta**4 - 138d0*costheta**2 + 3d0)
					return
				case ( -8, 8 )
					res = 3d0 / 2048d0 * prefac * dsqrt( 4849845d0 / ( 2d0*pi ) ) * &
						sintheta**8 * costheta * ( 115d0*costheta**4 - 46d0*costheta**2 + 3d0)
					return
				case ( -7, 7 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 692835d0 / ( pi ) ) * &
						sintheta**7 * ( 805d0*costheta**6 - 483d0*costheta**4 + 63d0*costheta**2 - 1d0)
					return
				case ( -6, 6 )
					res = 3d0 / 2048d0 * prefac * dsqrt( 969969d0 / ( pi ) ) * &
						sintheta**6 * costheta * ( 575d0*costheta**6 - 483d0*costheta**4 + 105d0*costheta**2 - 5d0)
					return
				case ( -5, 5 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 51051d0 / ( 2d0*pi ) ) * &
						sintheta**5 * ( 10925d0*costheta**8 - 12236d0*costheta**6 + 3990d0*costheta**4 - 380d0*costheta**2 + 5d0)
					return
				case ( -4, 4 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 51051d0 / ( pi ) ) * &
						sintheta**4 * costheta * ( 10925d0*costheta**8 - 15732d0*costheta**6 + 7182d0*costheta**4 - 1140d0*costheta**2 + 45d0)
					return
				case ( -3, 3 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 15015d0 / ( 2d0*pi ) ) * &
						sintheta**3 * ( 37145d0*costheta**10 - 66861d0*costheta**8 + 40698d0*costheta**6 - 9690d0*costheta**4 + 765d0*costheta**2 - 9d0)
					return
				case ( -2, 2 )
					res = 3d0 / 1024d0 * prefac * dsqrt( 1365d0 / ( 2d0*pi ) ) * &
						sintheta**2 * costheta * ( 37145d0*costheta**10 - 81719d0*costheta**8 + 63954d0*costheta**6 - 21318d0*costheta**4 + 2805d0*costheta**2 - 99d0)			
					return
				case ( -1, 1 )
					res = 3d0 / 2048d0 * prefac * dsqrt( 273d0 / ( 2d0*pi ) ) * &
						sintheta * ( 185725d0*costheta**12 - 490314d0*costheta**10 + 479655d0*costheta**8 - 213180d0*costheta**6 + 42075d0*costheta**4 - 2970d0*costheta**2 + 33d0)
					return
				case ( 0 )
					res = 3d0/2048d0 * dsqrt(3d0/pi) * ( &
						1300075d0*costheta**13 - 4056234d0*costheta**11 + 4849845d0*costheta**9 - 2771340d0*costheta**7 + 765765d0*costheta**5 - 90090d0*costheta**3 + 3003d0*costheta)
					return
			end select
				
		case ( 14 )
		
			select case ( m )
				case ( -14, 14 )
					res = 15d0 / 8192d0 * prefac * dsqrt( 646323d0 / ( 2d0*pi ) ) * &
						sintheta**14
					return
				case ( -13, 13 )
					res = 15d0 / 4096d0 * prefac * dsqrt( 4524261d0 / ( 2d0*pi ) ) * &
						sintheta**13 * costheta
					return
				case ( -12, 12 )
					res = 5d0 / 8192d0 * prefac * dsqrt( 1508087d0 / ( pi ) ) * &
						sintheta**12 * ( 27d0*costheta**2 - 1d0)
					return
				case ( -11, 11 )
					res = 5d0 / 4096d0 * prefac * dsqrt( 58815393d0 / ( 2d0*pi ) ) * &
						sintheta**11 * costheta * ( 9d0*costheta**2 - 1d0)
					return
				case ( -10, 10 )
					res = 1d0 / 8192d0 * prefac * dsqrt( 58815393d0 / ( 2d0*pi ) ) * &
						sintheta**10 * ( 225d0*costheta**4 - 50d0*costheta**2 + 1d0)
					return
				case ( -9, 9 )
					res = 1d0 / 4096d0 * prefac * dsqrt( 98025655d0 / ( pi ) ) * &
						sintheta**9 * costheta * ( 135d0*costheta**4 - 50d0*costheta**2 + 3d0)
					return
				case ( -8, 8 )
					res = 1d0 / 4096d0 * prefac * dsqrt( 12785955d0 / ( 2d0*pi ) ) * &
						sintheta**8 * ( 1035d0*costheta**6 - 575d0*costheta**4 + 69d0*costheta**2 - 1d0)
					return
				case ( -7, 7 )
					res = 1d0 / 4096d0 * prefac * dsqrt( 20092215d0 / ( pi ) ) * &
						sintheta**7 * costheta * ( 1035d0*costheta**6 - 805d0*costheta**4 + 161d0*costheta**2 - 7d0)
					return
				case ( -6, 6 )
					res = 1d0 / 8192d0 * prefac * dsqrt( 46881835d0 / ( 2d0*pi ) ) * &
						sintheta**6 * ( 3105d0*costheta**8 - 3220d0*costheta**6 + 966d0*costheta**4 - 84d0*costheta**2 + 1d0)
					return
				case ( -5, 5 )
					res = 3d0 / 4096d0 * prefac * dsqrt( 9376367d0 / ( 2d0*pi ) ) * &
						sintheta**5 * costheta * ( 1725d0*costheta**8 - 2300d0*costheta**6 + 966d0*costheta**4 - 140d0*costheta**2 + 5d0)
					return
				case ( -4, 4 )
					res = 3d0 / 8192d0 * prefac * dsqrt( 2467465d0 / ( pi ) ) * &
						sintheta**4 * ( 6555d0*costheta**10 - 10925d0*costheta**8 + 6118d0*costheta**6 - 1330d0*costheta**4 + 95d0*costheta**2 - 1d0)
					return
				case ( -3, 3 )
					res = 1d0 / 4096d0 * prefac * dsqrt( 224315d0 / ( 2d0*pi ) ) * &
						sintheta**3 * costheta * ( 58995d0*costheta**10 - 120175d0*costheta**8 + 86526d0*costheta**6 - 26334d0*costheta**4 + 3135d0*costheta**2 - 99d0)
					return
				case ( -2, 2 )
					res = 1d0 / 8192d0 * prefac * dsqrt( 39585d0 / ( 2d0*pi ) ) * &
						sintheta**2 * ( 334305d0*costheta**12 - 817190d0*costheta**10 + 735471d0*costheta**8 - 298452d0*costheta**6 + 53295d0*costheta**4 - 3366d0*costheta**2 + 33d0)
					return
				case ( -1, 1 )
					res = 1d0 / 2048d0 * prefac * dsqrt( 3045d0 / ( 2d0*pi ) ) * &
						sintheta * costheta * ( 334305d0*costheta**12 - 965770d0*costheta**10 + 1062347d0*costheta**8 - 554268d0*costheta**6 + 138567d0*costheta**4 - 14586d0*costheta**2 + 429d0)
					return
				case ( 0 )
					res = 1d0/4096d0 * dsqrt(29d0/pi) * ( &
						5014575d0*costheta**14 - 16900975d0*costheta**12 + 22309287d0*costheta**10 - 14549535d0*costheta**8 + 4849845d0*costheta**6 - 765765d0*costheta**4 + 45045d0*costheta**2 - 429d0)
					return
			end select
			
		case ( 15 )
		
			select case ( m )
				case ( -15, 15 )
					res = 3d0 / 16384d0 * prefac * dsqrt( 33393355d0 / ( pi ) ) * &
						sintheta**15
					return
				case ( -14, 14 )
					res = 15d0 / 8192d0 * prefac * dsqrt( 20036013d0 / ( 2d0*pi ) ) * &
						sintheta**14 * costheta 
					return
				case ( -13, 13 )
					res = 15d0 / 16384d0 * prefac * dsqrt( 690897d0 / ( pi ) ) * &
						sintheta**13 * ( 29d0*costheta**2 - 1d0)
					return
				case ( -12, 12 )
					res = 15d0 / 8192d0 * prefac * dsqrt( 1612093d0 / ( pi ) ) * &
						sintheta**12 * costheta * ( 29d0*costheta**2 - 3d0)
					return
				case ( -11, 11 )
					res = 5d0 / 16384d0 * prefac * dsqrt( 4836279d0 / ( pi ) ) * &
						sintheta**11 * ( 261d0*costheta**4 - 54d0*costheta**2 + 1d0)
					return
				case ( -10, 10 )
					res = 1d0 / 8192d0 * prefac * dsqrt( 314358135d0 / ( 2d0*pi ) ) * &
						sintheta**10 * costheta * ( 261d0*costheta**4 - 90d0*costheta**2 + 5d0)
					return
				case ( -9, 9 )
					res = 1d0 / 16384d0 * prefac * dsqrt( 104786045d0 / ( pi ) ) * &
						sintheta**9 * ( 1305d0*costheta**6 - 675d0*costheta**4 + 75d0*costheta**2 - 1d0)
					return
				case ( -8, 8 )
					res = 1d0 / 4096d0 * prefac * dsqrt( 44908305d0 / ( 2d0*pi ) ) * &
						sintheta**8 * costheta * ( 1305d0*costheta**6 - 945d0*costheta**4 + 175d0*costheta**2 - 7d0)
					return
				case ( -7, 7 )
					res = 1d0 / 16384d0 * prefac * dsqrt( 1952535d0 / ( pi ) ) * &
						sintheta**7 * ( 30015d0*costheta**8 - 28980d0*costheta**6 + 8050d0*costheta**4 - 644d0*costheta**2 + 7d0)
					return
				case ( -6, 6 )
					res = 1d0 / 8192d0 * prefac * dsqrt( 21477885d0 / ( 2d0*pi ) ) * &
						sintheta**6 * costheta * ( 10005d0*costheta**8 - 12420d0*costheta**6 + 4830d0*costheta**4 - 644d0*costheta**2 + 21d0)
					return
				case ( -5, 5 )
					res = 3d0 / 16384d0 * prefac * dsqrt( 10023013d0 / ( pi ) ) * &
						sintheta**5 * ( 10005d0*costheta**10 - 15525d0*costheta**8 + 8050d0*costheta**6 - 1610d0*costheta**4 + 105d0*costheta**2 - 1d0)
					return
				case ( -4, 4 )
					res = 3d0 / 8192d0 * prefac * dsqrt( 4555915d0 / ( pi ) ) * &
						sintheta**4 * costheta * ( 10005d0*costheta**10 - 18975d0*costheta**8 + 12650d0*costheta**6 - 3542d0*costheta**4 + 385d0*costheta**2 - 11d0)
					return
				case ( -3, 3 )
					res = 1d0 / 16384d0 * prefac * dsqrt( 719355d0 / ( pi ) ) * &
						sintheta**3 * ( 190095d0*costheta**12 - 432630d0*costheta**10 + 360525d0*costheta**8 - 134596d0*costheta**6 + 21945d0*costheta**4 - 1254d0*costheta**2 + 11d0)
					return
				case ( -2, 2 )
					res = 1d0 / 8192d0 * prefac * dsqrt( 55335d0 / ( 2d0*pi ) ) * &
						sintheta**2 * costheta * ( 570285d0*costheta**12 - 1533870d0*costheta**10 + 1562275d0*costheta**8 - 749892d0*costheta**6 + 171171d0*costheta**4 - 16302d0*costheta**2 + 429d0)
					return
				case ( -1, 1 )
					res = 1d0 / 16384d0 * prefac * dsqrt( 465d0 / ( pi ) ) * &
						sintheta * ( 9694845d0*costheta**14 - 30421755d0*costheta**12 + 37182145d0*costheta**10 - 22309287d0*costheta**8 + 6789783d0*costheta**6 - 969969d0*costheta**4 + 51051d0*costheta**2 - 429d0)
					return
				case ( 0 )
					res = 1d0/4096d0 * dsqrt(31d0/pi) * ( &
						9694845d0*costheta**15 - 35102025d0*costheta**13 + 50702925d0*costheta**11 - 37182145d0*costheta**9 + 14549535d0*costheta**7 - 2909907d0*costheta**5 + 255255d0*costheta**3 - 6435d0*costheta**1)
					return
			end select
			
	end select

	if ( ( m .gt. l ) .OR. ( m .lt. -l ) ) then
		write(*,*) "ERROR! Invalid spherical harmonic!"
		write(*,*)
		stop
	end if
	
	write(*,*) "ERROR! Spherical harmonic not implemented!"
	write(*,*)
	stop
	
end function

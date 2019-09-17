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

module data_mo
	use global_c
	implicit none
	
	character(sMaxBuffer) :: FTitle
	
	type t_atom
		integer :: id
		integer :: charge
		logical :: skip
		logical :: ghost
		character(5) :: str
		double precision :: x
		double precision :: y
		double precision :: z
		double precision :: theta
		double precision :: phi
		double precision :: dist
		double precision :: mulliken
		double precision :: ccharge
	end type t_atom
	
	type(t_atom), dimension(iMaxAtoms) :: Atoms
	integer :: NA = 0, NEl = 0
	
	type t_cgto
		integer :: id
		integer :: atomid
		integer :: shelltype
		integer :: npgto
		integer :: npgto_beg
		integer :: n
		integer :: subtype
		double precision :: ox, oy, oz
		logical :: deleted
		double precision, allocatable, dimension(:,:) :: coeff
	end type t_cgto
	
	type(t_cgto), allocatable, dimension(:) :: CGTOs
	integer :: NCGTO
	
	type t_pgto
		integer :: id
		integer :: atomid
		integer :: cgtoid
		integer :: cgtonpgto
		integer :: shelltype
		integer :: subtype
		integer :: subnum
		integer :: n
		double precision :: coeff, coeff_d, coeff_alpha
		double precision :: ox, oy, oz
		logical :: deleted
	end type t_pgto
	
	type(t_pgto), allocatable, dimension(:) :: PGTOs
	
	type t_mo
		integer :: id
		integer :: motype
		character(sMaxBuffer) :: name
		double precision :: energy
		double precision :: occup
		double precision, allocatable, dimension(:) :: mocoeff
	end type t_mo
	
	type(t_mo), allocatable, dimension(:) :: MOs
	integer :: NMO = 0
	
	double precision, dimension(:,:), allocatable :: AtomsDistMat, AtomsDistMatSq
	double precision, dimension(:,:), allocatable :: S_PGTO, S_CGTO, D_CGTO, O_CGTO12, DSMat, O_CGTO, O_CGTOT
	
	double precision, dimension(:,:), allocatable :: DensPMat
		
	double precision :: CSThresC, CSThresD, CSThresI
	
	real :: TimerA, TimerB

end module

subroutine print_bs_info
	use data_mo
	implicit none
!
	write(*,'(5X,A)') "> BASIS SET INFORMATION:"
	write(*,*)
	write(*,'(A35,I14)') "number of atoms: ", NA
	write(*,'(A35,I14)') "primitive GTOs: ", NumPGTO
	write(*,'(A35,I14)') "contracted GTOs: ", NumCGTO
	write(*,*)
	write(*,*)
	
end subroutine

subroutine check_ghost_atoms
	use data_mo
	implicit none
!
	integer :: i, j, nGhosts
!
	nGhosts = 0
  
	do j = 1, NA
		Atoms(j)%ghost = .TRUE.
		
		do i = 1, NumCGTO
			if ( CGTOs(i)%atomid .eq. j ) then
				Atoms(j)%ghost = .FALSE.
				exit
			end if
		end do
		
		if ( Atoms(j)%ghost .eqv. .TRUE. ) then
			nGhosts = nGhosts + 1
			Atoms(j)%charge = 0d0
		end if
	end do
	
	if ( nGhosts .gt. 0 ) then
		write(*,'(A35,I14)') "ghost atom(s) found: ", nGhosts
	end if
	
	if ( ODelete4f .eqv. .TRUE. ) then
		write(*,*) "> REMOVING 4f SHELL CONTRIBUTIONS (AS REQUESTED) ..."
		write(*,*)
		do j = 1, NumCGTO
			if ( ( Atoms(CGTOs(j)%atomid)%dist .lt. 0.01d0 ) ) then
			
				Atoms(CGTOs(j)%atomid)%ghost = .TRUE.
				
				! delete 4f shell
				if ( ( ODelete4f .eqv. .TRUE. ) .and. ( CGTOs(j)%shelltype .eq. 4 ) .and. ( CGTOs(j)%n .eq. 4 ) ) then
					write(*,'(10X,A,I4,A)') " ... CGTO #", j, " has been removed!"
					
					do i = 1, CGTOs(j)%npgto
						  CGTOs(j)%coeff(3,i) = 0d0
						  CGTOs(j)%coeff(2,i) = 0d0
						  CGTOs(j)%coeff(1,i) = 0d0
						  CGTOs(j)%deleted = .TRUE.
					end do
					
					do i = CGTOs(j)%npgto_beg, CGTOs(j)%npgto_beg + CGTOs(j)%npgto - 1
						PGTOs(i)%coeff_alpha = 0d0
						PGTOs(i)%coeff_d = 0d0
						PGTOs(i)%coeff = 0d0
						PGTOs(i)%deleted = .TRUE.
					end do
				end if
				
			end if
		end do
	else if ( ODeleteAll .eqv. .TRUE. ) then
		write(*,*) "> REMOVING ALL SHELLS FROM DY(III) CENTER (AS REQUESTED) ..."
		write(*,*)
		do j = 1, NumCGTO
			if ( ( Atoms(CGTOs(j)%atomid)%dist .lt. 0.01d0 ) ) then
			
				Atoms(CGTOs(j)%atomid)%ghost = .TRUE.
				
				! delete shell
				write(*,'(10X,A,I4,A)') " ... CGTO #", j, " has been removed!"
					
				do i = 1, CGTOs(j)%npgto
					  CGTOs(j)%deleted = .TRUE.
				end do
				
				do i = CGTOs(j)%npgto_beg, CGTOs(j)%npgto_beg + CGTOs(j)%npgto - 1
					PGTOs(i)%deleted = .TRUE.
				end do
				
			end if
		end do
! 	else
! 		do j = 1, NumCGTO
! 			if ( ( Atoms(CGTOs(j)%atomid)%dist .lt. 0.01d0 ) ) then
! 			
! 				Atoms(CGTOs(j)%atomid)%ghost = .TRUE.
! 			end if
! 		end do
	end if
	
	
	if ( ( ODelete4f .eqv. .TRUE. ) .or. ( ODeleteAll .eqv. .TRUE. ) ) then
		write(*,*)
		write(*,*)
	end if
	
end subroutine

subroutine read_molden_file(filename)
	use data_mo
	use global_c
	implicit none
!
	character (len=sMaxBuffer), intent(in) :: filename
	character (len=sMaxBuffer) :: to_upper, to_upper_first
	integer :: iost, iAtom, iBuffer, iPGTO, i, j, oshell, onshell, line = 0, mode = 0
	double precision :: dBuffer, tdensity
	character(len=sMaxBuffer) :: sBuffer, sBuffer2, sBuffer3, sShellType
	logical :: SkipR = .FALSE.
!
	tdensity = 0.d0
	NCGTO = 0
	
	write(*,'(/5X,3A/)',advance='no') "> READING MOLDEN FILE '", trim(filename), "' ..."
	
	open(unit=uMolF, file=filename, action='read', status='old', iostat=iost)
	if ( iost .ne. 0 ) then
		write(*,*) "ERROR! Can't open molden file!"
		stop
	end if
	
	! read file header
	read(uMolF, '(A)') sBuffer
	line = line + 1
	sBuffer = to_upper(sBuffer)
	if ( sBuffer .ne. "[MOLDEN FORMAT]" ) then
		write(*,*)
		write(*,*) "ERROR! This is not a molden file!"
		close(uMolF)
		stop
	end if
	do
		if ( SkipR .eqv. .FALSE. ) then
			line = line + 1
			read(uMolF, '(A)',iostat=iost) sBuffer
		else
			SkipR = .FALSE.
		end if
		
		if ( iost .ne. 0 ) then
			! end of file
			exit
		else
		
			sBuffer = to_upper(sBuffer)
			
			if ( len_trim(sBuffer) .ne. 0 ) then

				if ( sBuffer .eq. "[TITLE]" ) then
					line = line + 1
					read(uMolF, '(A)', iostat=iost) FTitle
					continue
				else if ( sBuffer(1:7) .eq. "[ATOMS]" ) then
				
					! TODO: check for AU keyword				
					mode = 1
					
					do
						line = line + 1
						read(uMolF, '(A)',iostat=iost) sBuffer
						
						if ( sBuffer(1:1) .ne. "[" ) then
						
							NA = NA + 1
							read(sBuffer, *) Atoms(NA)%str, Atoms(NA)%id, Atoms(NA)%charge, Atoms(NA)%x, Atoms(NA)%y, Atoms(NA)%z
							
							Atoms(NA)%str = to_upper_first(Atoms(NA)%str)
							Atoms(NA)%skip = .FALSE.
							Atoms(NA)%ghost = .FALSE.
							
						else
							SkipR = .TRUE.
							exit
						end if
					
					end do
					
					continue
				
				else if ( sBuffer .eq. "[STO]" ) then
				
					write(*,*)
					write(*,*) "ERROR! STO basis sets are not supported!"
					close(uMolF)
					stop
					
				else if ( sBuffer .eq. "[GTO]" ) then
				
					if ( mode .ne. 1 ) then
						write(*,*)
						write(*,*) "ERROR! Molden file is corrupt!"
						close(uMolF)
						stop
					end if
					mode = 2

					allocate(CGTOs(iMaxCGTOs),stat=iost)
					
					do
					
						line = line + 1
						read(uMolF, '(A)',iostat=iost) sBuffer
						if ( len_trim(sBuffer) .ne. 0 ) then
						
							if ( sBuffer(1:1) .ne. "[" ) then
							
								read(sBuffer, *) iAtom, iBuffer
								oshell = 0
								onshell = 0
								
								do
								
									line = line + 1
									read(uMolF, '(A)',iostat=iost) sBuffer
								
									if ( len_trim(sBuffer) .ne. 0 ) then
									
										NCGTO = NCGTO + 1
										CGTOs(NCGTO)%atomid = iAtom
										CGTOs(NCGTO)%ox = Atoms(iAtom)%x
										CGTOs(NCGTO)%oy = Atoms(iAtom)%y
										CGTOs(NCGTO)%oz = Atoms(iAtom)%z
										CGTOs(NCGTO)%id = NCGTO
										
										read(sBuffer, *) sShellType, iPGTO, dBuffer
										
										
										sShellType = to_upper(sShellType)
										select case ( sShellType ) 
											case ( "S" )
												NPGTO = NPGTO + 1
												CGTOs(NCGTO)%shelltype = 1
											case ( "P" )
												NPGTO = NPGTO + 3
												CGTOs(NCGTO)%shelltype = 2
											case ( "D" )
												NPGTO = NPGTO + 6
												CGTOs(NCGTO)%shelltype = 3
											case ( "F" )
												NPGTO = NPGTO + 10
												CGTOs(NCGTO)%shelltype = 4
											case ( "G" ) 
												NPGTO = NPGTO + 15
												CGTOs(NCGTO)%shelltype = 5
											case default
												write(*,*)
												write(*,*) "ERROR! Only shell type S,P,D,F,G are supported!"
												close(uMolF)
												stop
										end select
										
										if ( CGTOs(NCGTO)%shelltype .eq. oshell ) then
										
											onshell = onshell + 1
										
										else
										
											oshell = CGTOs(NCGTO)%shelltype
											onshell = CGTOs(NCGTO)%shelltype 
										
										end if
										
										oshell = CGTOs(NCGTO)%shelltype
										
										CGTOs(NCGTO)%n = onshell
										CGTOs(NCGTO)%npgto = iPGTO
										allocate(CGTOs(NCGTO)%coeff(3,iPGTO))
										
										do i = 1, iPGTO
										
											line = line + 1
											read(uMolF, '(A)',iostat=iost) sBuffer
											read(sBuffer,*) CGTOs(NCGTO)%coeff(1, i), CGTOs(NCGTO)%coeff(2, i)
											select case ( CGTOs(NCGTO)%shelltype )
												case ( 1 )
													CGTOs(NCGTO)%coeff(3, i) = CGTOs(NCGTO)%coeff(2, i)
												case ( 2 )
													CGTOs(NCGTO)%coeff(3, i) = CGTOs(NCGTO)%coeff(2, i)
												case ( 3 ) 
													CGTOs(NCGTO)%coeff(3, i) = CGTOs(NCGTO)%coeff(2, i)
												case ( 4 ) 
													CGTOs(NCGTO)%coeff(3, i) = CGTOs(NCGTO)%coeff(2, i)
												case ( 5 ) 
													CGTOs(NCGTO)%coeff(3, i) = CGTOs(NCGTO)%coeff(2, i) 
												case default
													write(*,*)
													write(*,*) "ERROR! Unsupported shell type found!"
													close(uMolF)
													stop
											end select
										
										end do
										
									else
									
										exit
										
									end if
								end do
							else
								exit
							end if
						end if
						
					end do
					
					SkipR = .TRUE.
										
				else if ( sBuffer .eq. "[MO]" ) then
				
					if ( mode .ne. 2 ) then
						write(*,*)
						write(*,*) "ERROR! Molden file seems to be invalid!"
						close(uMolF)
						stop
					end if
					
					mode = 3
					allocate(MOs(iMaxMOs))
					
					do
					
						line = line + 1
						read(uMolF, '(A)',iostat=iost) sBuffer
						if ( iost < 0 ) then
						
							! end of file
							exit
						
						end if
						if ( len_trim(sBuffer) .ne. 0 ) then
						
							NMO = NMO + 1
							sBuffer = adjustl(to_upper(sBuffer))
							read(sBuffer, *) sBuffer2, sBuffer3
						
							if ( sBuffer2 .ne. "SYM=" ) then
								write(*,*)
								write(*,*) "ERROR! Can't read MO parameter (1)!", sBuffer2
								close(uMolF)
								stop
							
							end if
							
							MOs(NMO)%id = NMO
							MOs(NMO)%name = sBuffer3
							
							line = line + 1
							read(uMolF, '(A)') sBuffer
							sBuffer = adjustl(to_upper(sBuffer))
							read(sBuffer, *) sBuffer2, sBuffer3
							
							if ( sBuffer2 .ne. "ENE=" ) then
								write(*,*)
								write(*,*) "ERROR! Can't read MO parameters (2)!"
								close(uMolF)
								stop
							
							end if
							
							read(sBuffer3, *) MOs(NMO)%energy
							
							line = line + 1
							read(uMolF, '(A)') sBuffer
							sBuffer = adjustl(to_upper(sBuffer))
							read(sBuffer, *) sBuffer2, sBuffer3
							
							if ( sBuffer2 .ne. "SPIN=" ) then
								write(*,*)
								write(*,*) "ERROR! Can't read MO parameters (3)!"
								close(uMolF)
								stop
							
							end if
							
							if ( sBuffer3 .eq. "ALPHA" ) then
								MOs(NMO)%motype = 1
							else
								MOs(NMO)%motype = 2
								if ( OUHF .eqv. .FALSE. ) then
									OUHF = .TRUE.
									write(*,*)
									write(*,'(7X,A)') "An UHF/UKS wave function was found!"
									write(*,*)
								end if
							end if
							
							line = line + 1
							read(uMolF, '(A)') sBuffer
							sBuffer = adjustl(to_upper(sBuffer))
							read(sBuffer, *) sBuffer2, sBuffer3
							
							if ( sBuffer2 .ne. "OCCUP=" ) then
								write(*,*)
								write(*,*) "ERROR! Can't read MO parameters (4)!"
								close(uMolF)
								stop
							
							end if
							
							read(sBuffer3, *) MOs(NMO)%occup
							allocate(MOs(NMO)%mocoeff(NPGTO))
							
							do i = 1, NPGTO
							
								line = line + 1
								read(uMolF, '(A)') sBuffer
								read(sBuffer, *) iBuffer, dBuffer
								
								if ( iBuffer .ne. i ) then
									write(*,*)
									write(*,*) "ERROR! Can't read MO parameters (5)!"
									close(uMolF)
									stop
								
								end if
								
								MOs(NMO)%mocoeff(i) = dBuffer
							
							end do
! 							
						else
							! empty line
							exit
						end if
					end do
				else
				
					write(*,*) "WARNING! Unknown keyword in Molden file found!"
					write(*,*) sBuffer
					close(uMolF)
					stop
					
				end if
			
			end if
		
		end if
	
	end do
	
	close(uMolF)
		
	allocate(AtomsDistMat(NA,NA))
	allocate(AtomsDistMatSq(NA,NA))
	
	AtomsDistMat = 0d0
	AtomsDistMatSq = 0d0
	
	do i = 1, NA
		do j = i + 1, NA
			AtomsDistMatSq(i, j) = (Atoms(i)%x - Atoms(j)%x)**2 + (Atoms(i)%y - Atoms(j)%y)**2 + (Atoms(i)%z - Atoms(j)%z)**2
			if ( AtomsDistMatSq(i, j) .lt. 1d-15 ) then
				AtomsDistMatSq(i, j) = 0d0
			end if
			
			AtomsDistMatSq(j, i) = AtomsDistMatSq(i, j)
			
			AtomsDistMat(i, j) = dsqrt( AtomsDistMatSq(i, j) )
			AtomsDistMat(j, i) = AtomsDistMat(i, j)
		end do
	end do
	
	TotalNuclearCharge = 0.0d0
	
	do i = 1, NA
		if ( Atoms(i)%ghost .eqv. .FALSE. ) then
			TotalNuclearCharge = TotalNuclearCharge + Atoms(i)%charge
		end if
	
	end do
	
	if ( CAtomID .gt. 0 ) then
		EDyX = Atoms(CAtomID)%x
		EDyY = Atoms(CAtomID)%y
		EDyZ = Atoms(CAtomID)%z
	end if
	
	do i = 1, NA
		Atoms(i)%dist = dsqrt( (EDyX - Atoms(i)%x )**2 + (EDyY - Atoms(i)%y )**2 + (EDyZ - Atoms(i)%z )**2 )
		Atoms(i)%theta = dacos( (Atoms(i)%z-EDyZ) / Atoms(i)%dist )
		Atoms(i)%phi = datan2( Atoms(i)%y-EDyY, Atoms(i)%x-EDyX )
	end do
	
	call write_done

end subroutine

function get_pgto_norm_factor(alpha, nx, ny, nz) result(res)
	use data_mo
	implicit none
!
	double precision, intent(in) :: alpha
	integer, intent(in) :: nx, ny, nz
	integer :: n_fac
	double precision :: res
!
	
	res = 	(2d0*alpha/Pi)**(3d0/4d0) * &
			dsqrt( (8d0*alpha)**(nx+ny+nz) * n_fac(nx) * n_fac(ny) * n_fac(nz) / &
			( n_fac(2*nx) * n_fac(2*ny) * n_fac(2*nz) ) )

end function

subroutine reorder_cgtos
	use data_mo
	implicit none
!
	integer :: i, j, k, l, m, n, nprim, ncont, nx, ny, nz, pgto_get_l
	double precision :: get_pgto_norm_factor
!
	nprim = 0
	ncont = 0
	
	! count contracted and primitive GTOs, cartesian only
	do i = 1, NCGTO
	
		select case ( CGTOs(i)%shelltype )
		
			case ( 1 )
				nprim = nprim + 1 * CGTOs(i)%npgto
				ncont = ncont + 1
			case ( 2 )
				nprim = nprim + 3 * CGTOs(i)%npgto
				ncont = ncont + 3
			case ( 3 )
				nprim = nprim + 6 * CGTOs(i)%npgto
				ncont = ncont + 6
			case ( 4 )
				nprim = nprim + 10 * CGTOs(i)%npgto
				ncont = ncont + 10
			case ( 5 )
				nprim = nprim + 15 * CGTOs(i)%npgto
				ncont = ncont + 15
			case default
				write(*,*) "ERROR! Shell type not implemented!"
				stop
		
		end select
	
	end do
	
	NumPGTO = nprim
	NumCGTO = ncont
	
	if ( NumCGTO > iMaxCGTOs ) then
		write(*,*) "ERROR! Too many CGTOs! Try to increase iMaxCGTOs!"
		stop
	end if
	
	allocate(PGTOs(NumPGTO))
	
	j = 0
	n = 1
	
	do i = 1, NumCGTO
	
		select case ( CGTOs(i)%shelltype )
		
			case ( 1 )
				k = 1
			case ( 2 )
				k = 3
			case ( 3 )
				k = 6
			case ( 4 )
				k = 10
			case ( 5 )
				k = 15
			case default
				k = 0
		end select
		
		do l = 1, k
			
			do m = 1, CGTOs(i)%npgto
			
				j = j + 1

				PGTOs(j)%id = j
				PGTOs(j)%atomid = CGTOs(i)%atomid
				PGTOs(j)%cgtoid = n
				PGTOs(j)%shelltype = CGTOs(i)%shelltype
				
				PGTOs(j)%ox = Atoms(CGTOs(i)%atomid)%x 
				PGTOs(j)%oy = Atoms(CGTOs(i)%atomid)%y 
				PGTOs(j)%oz = Atoms(CGTOs(i)%atomid)%z 
			
				PGTOs(j)%coeff_alpha = CGTOs(i)%coeff(1, m)
				PGTOs(j)%coeff_d = CGTOs(i)%coeff(2, m)
				PGTOs(j)%coeff = CGTOs(i)%coeff(2, m) 
				PGTOs(j)%cgtonpgto = CGTOs(i)%npgto
				PGTOs(j)%subtype = l
				PGTOs(j)%subnum = m
				PGTOs(j)%n = CGTOs(i)%n
				PGTOs(j)%deleted = .FALSE.
				
				nx = pgto_get_l(PGTOs(j)%shelltype, PGTOs(j)%subtype, 1)
				ny = pgto_get_l(PGTOs(j)%shelltype, PGTOs(j)%subtype, 2)
				nz = pgto_get_l(PGTOs(j)%shelltype, PGTOs(j)%subtype, 3)
				
				PGTOs(j)%coeff_d = PGTOs(j)%coeff_d * get_pgto_norm_factor(PGTOs(j)%coeff_alpha, nx, ny, nz)
				
				! MOLDEN files from TURBOMOLE need a correction
				select case ( PGTOs(j)%shelltype )
					case ( 3 ) 						
						PGTOs(j)%coeff_d = PGTOs(j)%coeff_d * dsqrt(3d0)
					case ( 4 )
						PGTOs(j)%coeff_d = PGTOs(j)%coeff_d * dsqrt(15d0)
					case ( 5 )
						PGTOs(j)%coeff_d = PGTOs(j)%coeff_d * dsqrt(105d0)
				end select				
				
			end do
			
			n = n + 1
			
		end do
		
	end do
	
	do i = 1, NCGTO	
		deallocate(CGTOs(i)%coeff)	
	end do
	
	deallocate(CGTOs)
	allocate(CGTOs(NumCGTO))
	
	do i = 1, NumCGTO
		do j = 1, NumPGTO
			if ( PGTOs(j)%cgtoid .eq. i ) then
			
				CGTOs(i)%id = PGTOs(j)%cgtoid
				CGTOs(i)%atomid = PGTOs(j)%atomid
				CGTOs(i)%shelltype = PGTOs(j)%shelltype
				CGTOs(i)%npgto = PGTOs(j)%cgtonpgto
				CGTOs(i)%ox = PGTOs(j)%ox
				CGTOs(i)%oy = PGTOs(j)%oy
				CGTOs(i)%oz = PGTOs(j)%oz
				CGTOs(i)%n = PGTOs(j)%n
				CGTOs(i)%subtype = PGTOs(j)%subtype
				if ( PGTOs(j)%subnum .eq. 1 ) then
					CGTOs(i)%npgto_beg = j
				end if
				CGTOs(i)%deleted = .FALSE.
				
				allocate(CGTOs(i)%coeff(3, CGTOs(i)%npgto))
				exit
			
			end if
		end do
	end do
	
	do i = 1, NumPGTO
		CGTOs(PGTOs(i)%cgtoid)%coeff(1, PGTOs(i)%subnum) = PGTOs(i)%coeff_alpha
		CGTOs(PGTOs(i)%cgtoid)%coeff(2, PGTOs(i)%subnum) = PGTOs(i)%coeff
	end do
	
end subroutine

subroutine print_pgto_info
	use data_mo
	implicit none
!
	double precision :: res
	integer :: i
	character (len=1) :: stype
	character (len=8) :: subtype
!
	res = 0d0
	i = 0
		
	write(*,*) "PRIMITIVE GTO BASIS SET INFO:"
	write(*,*)
	
	write(*,'(8X,A)') "---------------------------------------------------------------------"
	write(*,'(8X,A5,4A8,2A16)') "#", "CGTO", "ATOM", "SHELL", "TYPE", "ALPHA", "D"
	write(*,'(8X,A)') "---------------------------------------------------------------------"
	
	do while ( i .lt. NumPGTO )
	
		i = i + 1
		
		stype = "?"
		subtype = " "
		
		select case ( PGTOs(i)%shelltype )
			case ( 1 )
				stype = "s"
			case ( 2 )
				stype = "p"
				select case ( PGTOs(i)%subtype ) 
					case ( 1 )
						subtype = "x"
					case ( 2 )
						subtype = "y"
					case ( 3 )
						subtype = "z"
					case default
						subtype = "?"
				end select
			case ( 3 )
				stype = "d"
				select case ( PGTOs(i)%subtype ) 
					case ( 1 )
						subtype = "xx"
					case ( 2 )
						subtype = "yy"
					case ( 3 )
						subtype = "zz"
					case ( 4 )
						subtype = "xy"
					case ( 5 )
						subtype = "xz"
					case ( 6 )
						subtype = "yz"
					case default
						subtype = "??"
				end select
			case ( 4 )
				stype = "f"
				select case ( PGTOs(i)%subtype ) 
					case ( 1 )
						subtype = "xxx"
					case ( 2 )
						subtype = "yyy"
					case ( 3 )
						subtype = "zzz"
					case ( 4 )
						subtype = "xyy"
					case ( 5 )
						subtype = "xxy"
					case ( 6 )
						subtype = "xxz"
					case ( 7 )
						subtype = "xzz"
					case ( 8 )
						subtype = "yzz"
					case ( 9 )
						subtype = "yyz"
					case ( 10 )
						subtype = "xyz"
					case default
						subtype = "???"
				end select
			case default
				stype = "?"
		
		end select
		
		stype = adjustr(stype)
		subtype = adjustr(subtype)
		
		write(*,'(8X,I5,2I8,6X,I1,A1,A8,2F16.6)') i, PGTOs(i)%cgtoid, PGTOs(i)%atomid, PGTOs(i)%n, stype, subtype, PGTOs(i)%coeff_alpha, PGTOs(i)%coeff
		
	end do

	write(*,'(8X,A)') "---------------------------------------------------------------------"
	write(*,*)
	
end subroutine

subroutine print_cgto_info
	use data_mo
	implicit none
!
	double precision :: res
	integer :: i
	character (len=1) :: stype
	character (len=8) :: subtype
!
	res = 0d0
	i = 0
		
	write(*,*) "CONTRACTED GTO BASIS SET INFO:"
	write(*,*)
	
	write(*,'(48X,A)') "-----------------------------"
	write(*,'(48X,A5,4A8,2A16)') "#", "ATOM", "SHELL", "TYPE"
	write(*,'(48X,A)') "-----------------------------"
	
	do while ( i .lt. NumCGTO )
	
		i = i + 1
		
		stype = "?"
		subtype = " "
		
		select case ( CGTOs(i)%shelltype )
			case ( 1 )
				stype = "s"
			case ( 2 )
				stype = "p"
				select case ( CGTOs(i)%subtype ) 
					case ( 1 )
						subtype = "x"
					case ( 2 )
						subtype = "y"
					case ( 3 )
						subtype = "z"
					case default
						subtype = "?"
				end select
			case ( 3 )
				stype = "d"
				select case ( CGTOs(i)%subtype ) 
					case ( 1 )
						subtype = "xx"
					case ( 2 )
						subtype = "yy"
					case ( 3 )
						subtype = "zz"
					case ( 4 )
						subtype = "xy"
					case ( 5 )
						subtype = "xz"
					case ( 6 )
						subtype = "yz"
					case default
						subtype = "??"
				end select
			case ( 4 )
				stype = "f"
				select case ( CGTOs(i)%subtype ) 
					case ( 1 )
						subtype = "xxx"
					case ( 2 )
						subtype = "yyy"
					case ( 3 )
						subtype = "zzz"
					case ( 4 )
						subtype = "xyy"
					case ( 5 )
						subtype = "xxy"
					case ( 6 )
						subtype = "xxz"
					case ( 7 )
						subtype = "xzz"
					case ( 8 )
						subtype = "yzz"
					case ( 9 )
						subtype = "yyz"
					case ( 10 )
						subtype = "xyz"
					case default
						subtype = "???"
				end select
			case default
				stype = "?"
		
		end select
		
		stype = adjustr(stype)
		subtype = adjustr(subtype)
		
		write(*,'(48X,I5,I8,6X,I1,A1,A8)') i, CGTOs(i)%atomid, CGTOs(i)%n, stype, subtype
		
	end do

	write(*,'(48X,A)') "-----------------------------"
	write(*,*)
	
end subroutine

function pgto_chi(mu, x, y, z) result(res)
	use data_mo
	implicit none
!
	integer, intent(in) :: mu
	double precision, intent(in) :: x, y, z
	double precision :: res
	integer :: i
	double precision :: x0, y0, z0, drho, r2
!
	res = 0d0
	drho = 0d0
	i = 0
	
	do while ( i .lt. NumPGTO )
	
		i = i + 1
		
		if ( PGTOs(i)%cgtoid .eq. mu ) then
		
			x0 = x - PGTOs(i)%ox 
			y0 = y - PGTOs(i)%oy 
			z0 = z - PGTOs(i)%oz 

			r2 = x0**2 + y0**2 + z0**2
		
			if ( PGTOs(i)%shelltype .gt. 1 ) then
		
				select case ( PGTOs(i)%shelltype ) 
				
					case ( 2 )
						select case ( PGTOs(i)%subtype )
							case ( 1 )
								drho = x0
							case ( 2 )
								drho = y0
							case ( 3 )
								drho = z0
						end select
					case ( 3 )
						select case ( PGTOs(i)%subtype )
						
							case ( 1 )
								drho = x0*x0
							case ( 2 )
								drho = y0*y0
							case ( 3 ) 
								drho = z0*z0
							case ( 4 )
								drho = x0*y0
							case ( 5 )
								drho = x0*z0
							case ( 6 )
								drho = y0*z0
						end select
					case ( 4 )
						select case ( PGTOs(i)%subtype )
							case ( 1 )
								drho = x0*x0*x0
							case ( 2 )
								drho = y0*y0*y0
							case ( 3 ) 
								drho = z0*z0*z0
							case ( 4 )
								drho = x0*y0*y0
							case ( 5 )
								drho = x0*x0*y0
							case ( 6 )
								drho = x0*x0*z0
							case ( 7 )
								drho = x0*z0*z0
							case ( 8 )
								drho = y0*z0*z0
							case ( 9 )
								drho = y0*y0*z0
							case ( 10 )
								drho = x0*y0*z0
						end select
					case ( 5 )
						select case ( PGTOs(i)%subtype )
							case ( 1 )
								drho = x0*x0*x0*x0 !xxxx
							case ( 2 )
								drho = y0*y0*y0*y0 !yyyy
							case ( 3 ) 
								drho = z0*z0*z0*z0 !zzzz
							case ( 4 )
								drho = x0*x0*x0*y0 !xxxy
							case ( 5 )
								drho = x0*x0*x0*z0 !xxxz
							case ( 6 )
								drho = y0*y0*y0*x0 !yyyx
							case ( 7 )
								drho = y0*y0*y0*z0 !yyyz
							case ( 8 )
								drho = z0*z0*z0*x0 !zzzx
							case ( 9 )
								drho = z0*z0*z0*y0 !zzzy
							case ( 10 )
								drho = x0*x0*y0*y0 !xxyy
							case ( 11 )
								drho = x0*x0*z0*z0 !xxzz
							case ( 12 )
								drho = y0*y0*z0*z0 !yyzz
							case ( 13 )
								drho = x0*x0*y0*z0 !xxyz
							case ( 14 )
								drho = y0*y0*x0*z0 !yyxz
							case ( 15 )
								drho = z0*z0*x0*y0 !zzxy
						end select
					case default
						write(*,*) "WARNING! Shell type not implemented!"
						
				end select
			
				res = res + drho * PGTOs(i)%coeff_d * dexp( - PGTOs(i)%coeff_alpha * r2 )
			
			else 
			
				res = res + PGTOs(i)%coeff_d * dexp( - PGTOs(i)%coeff_alpha * r2 )
			
			end if
			
		end if
	
	end do
	
	
end function

subroutine calc_density_point_dmat4(x, y, z, rho)
	use data_mo
	implicit none
!
	double precision, intent(in) :: x, y, z
	double precision, intent(out) :: rho
	integer :: i, j, l, aid
	double precision :: drho, x0, y0, z0, r2, x02, y02, z02, exppart
	double precision, parameter :: ThresDens = 1d-7
	double precision, parameter :: ThresDMat = 1d-7
	double precision, dimension(iMaxCGTOs) :: cgtoval
!
	rho = 0d0
	aid = 0
	
	do i = 1, NumCGTO
	
		if ( CGTOs(i)%atomid .ne. aid ) then
			x0 = x - CGTOs(i)%ox 
			y0 = y - CGTOs(i)%oy 
			z0 = z - CGTOs(i)%oz 
			
			x02 = x0**2
			y02 = y0**2
			z02 = z0**2
			
			r2 = x02 + y02 + z02
			aid = CGTOs(i)%atomid
		end if
		
		cgtoval(i) = 0d0
		drho = 0d0
		
		do l = CGTOs(i)%npgto_beg, CGTOs(i)%npgto_beg + CGTOs(i)%npgto - 1
			exppart = PGTOs(l)%coeff_d * dexp( - PGTOs(l)%coeff_alpha * r2 ) 
			if ( PGTOs(l)%shelltype .gt. 1 ) then
				select case ( PGTOs(l)%shelltype ) 
					case ( 2 ) ! p
						select case ( PGTOs(l)%subtype )
							case ( 1 )
								drho = x0 
							case ( 2 )
								drho = y0 
							case ( 3 )
								drho = z0 
						end select
					case ( 3 ) ! d 
						select case ( PGTOs(l)%subtype )
							case ( 1 )
								drho = x0*x0 !xx
							case ( 2 )
								drho = y0*y0 !yy
							case ( 3 )
								drho = z0*z0 !zz
							case ( 4 )
								drho = x0*y0 !xy
							case ( 5 )
								drho = x0*z0 !xz
							case ( 6 )
								drho = y0*z0 !yz
						end select
					case ( 4 )
						select case ( PGTOs(l)%subtype )
							case ( 1 )
								drho = x0*x0*x0 !xxx
							case ( 2 )
								drho = y0*y0*y0 !yyy
							case ( 3 ) 
								drho = z0*z0*z0 !zzz
							case ( 4 )
								drho = x0*y0*y0 !xyy
							case ( 5 )
								drho = x0*x0*y0 !xxy
							case ( 6 )
								drho = x0*x0*z0 !xxz
							case ( 7 )
								drho = x0*z0*z0 !xzz
							case ( 8 )
								drho = y0*z0*z0 !yzz
							case ( 9 )
								drho = y0*y0*z0 !yyz
							case ( 10 )
								drho = x0*y0*z0 !xyz
						end select
					case ( 5 )
						select case ( PGTOs(i)%subtype )
							case ( 1 )
								drho = x0*x0*x0*x0 !xxxx
							case ( 2 )
								drho = y0*y0*y0*y0 !yyyy
							case ( 3 ) 
								drho = z0*z0*z0*z0 !zzzz
							case ( 4 )
								drho = x0*x0*x0*y0 !xxxy
							case ( 5 )
								drho = x0*x0*x0*z0 !/ 105d0**2  !xxxz
							case ( 6 )
								drho = y0*y0*y0*x0 !yyyx
							case ( 7 )
								drho = y0*y0*y0*z0 !yyyz
							case ( 8 )
								drho = z0*z0*z0*x0 !/ 105d0**2  !zzzx
							case ( 9 )
								drho = z0*z0*z0*y0 !zzzy
							case ( 10 )
								drho = x0*x0*y0*y0 !xxyy
							case ( 11 )
								drho = x0*x0*z0*z0 !/ 105d0**2 !xxzz
							case ( 12 )
								drho = y0*y0*z0*z0 !yyzz
							case ( 13 )
								drho = x0*x0*y0*z0 !xxyz
							case ( 14 )
								drho = y0*y0*x0*z0 !yyxz
							case ( 15 )
								drho = z0*z0*x0*y0 !zzxy
						end select
					case default
						write(*,*) "Error! shell type not implemented!"
						stop
				end select
				cgtoval(i) = cgtoval(i) + drho * exppart
			else
				cgtoval(i) = cgtoval(i) + exppart
			end if
		end do
  
	end do
	
	do i = 1, NumCGTO
		do j = 1, NumCGTO
			rho = rho + 2d0 * cgtoval(i) * cgtoval(j) * D_CGTO(i, j) 
		end do
	end do 
	
end subroutine	
	
function n_fac(n) result(res)
	implicit none
	integer, intent(in) :: n
	integer :: res, i
	
	if ( n .ge. 0 ) then
	
		res = 1
	
		do i = 1, n
			res = res * i
		end do
	
	else
		res = 0
		write(*,*) "internal Error!", n
		stop 
	end if
	
end function

function n_fac2(n) result(res)
	implicit none
	integer, intent(in) :: n
	integer(8) :: res

!					 -1  0  1  2  3  4   5   6    7    8    9       10    11    12     13     14      15       16       17,       18,       19,        20,         21,         22
	integer(8), dimension(31) :: N_ = (/ &
		1, 1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, &
		10395,46080,135135,645120,2027025,10321920,34459425,185794560,654729075,3715891200,13749310575,81749606400, &
		316234143225,1961990553600,7905853580625,51011754393600,213458046676875,1428329123020800,6190283353629375 /)
	
	if ( ( n .ge. -1 ) .AND. ( n .le. 29 ) ) then
		res = N_( n + 2 )
	else
		write(*,*) "internal Error!1", n
		stop 
	end if
end function

function l_over_k(l, k) result(n)
	implicit none
	integer, intent(in) :: l, k
	integer :: n, n_fac
	
	if ( ( l .ge. 0 ) .AND. ( k .ge. 0 ) ) then
	
		n = n_fac( l ) / ( n_fac(k) * n_fac( l - k ) )
	
	else
		n = 0
		write(*,*) "internal Error!"
		stop 
	end if

end function

function pgto_get_l(shelltype, subtype, xyz) result(n)

	implicit none
	integer, intent(in) :: shelltype, subtype, xyz
	integer :: n
	
!						   xx, yy, zz, xy, xz, yz
	integer, dimension(3, 6) :: N3_ = reshape( (/ &
						(/  2,  0,  0,  1,  1,  0 /), &	! x
						(/  0,  2,  0,  1,  0,  1 /), &	! y
						(/  0,  0,  2,  0,  1,  1 /)  &	! z
					/), shape(N3_), order=(/2, 1/) )
	
! 						 xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
	integer, dimension(3,10) :: N4_ = reshape( (/  &
						(/ 3,   0,   0,   1,   2,   2,   1,   0,   0,   1 /), &		! x
						(/ 0,   3,   0,   2,   1,   0,   0,   1,   2,   1 /), &		! y
						(/ 0,   0,   3,   0,   0,   1,   2,   2,   1,   1 /)  &	! z
					/), shape(N4_), order=(/2, 1/) )
!   						xxxx, yyyy, zzzz, xxxy, xxxz, yyyx, yyyz, zzzx, zzzy, xxyy, xxzz, yyzz, xxyz, yyxz, zzxy
	integer, dimension(3,15) :: N5_ = reshape( (/  &
						(/ 4,    0,    0,    3,    3,    1,    0,    1,    0,    2,    2,    0,    2,    1,    1 /), &	! x
						(/ 0,    4,    0,    1,    0,    3,    3,    0,    1,    2,    0,    2,    1,    2,    1 /), &	! y
						(/ 0,    0,    4,    0,    1,    0,    1,    3,    3,    0,    2,    2,    1,    1,    2 /)  &
					/), shape(N5_), order=(/2, 1/) )
	
	
	if ( ( xyz .ge. 1 ) .AND. ( xyz .le. 3 ) ) then
	
		select case ( shelltype )
			case ( 1 )
				n = 0
				return
			case ( 2 ) 
				if ( subtype .eq. xyz ) then
					n = 1
				else
					n = 0
				end if
				return
			case ( 3 )
				if ( subtype .le. 6 ) then
					n = N3_(xyz, subtype)
					return
				end if
			case ( 4 )
				if ( subtype .le. 10 ) then
					n = N4_(xyz, subtype)
					return
				end if
			case ( 5 )
				if ( subtype .le. 15 ) then
					n = N5_(xyz, subtype)
					return
				end if
		end select
	end if

	write(*,*) "internal Error!"
	stop 
	
end function

function f_binom_factors(j, l, m, a, b) result(res)
	implicit none
!
	integer, intent(in) :: j, l, m
	double precision, intent(in) :: a, b
	double precision :: res
	integer :: l_over_k
	integer :: k1
!
	res = 0d0

	do k1 = max(0, j-m), min(j, l)
	
		res = res + l_over_k(l, k1) * l_over_k(m, j-k1) * a**(l-k1) * b**(m+k1-j)
	
	end do
	
end function

function G_binom_factors(J, gl1, gl2, gax, gbx, cx, px, gamma) result(res)
	implicit none
!
	integer, intent(in) :: J, gl1, gl2
	double precision, intent(in) :: gax, gbx, cx, px, gamma
	integer :: i, r, u, n_fac, l1, l2
	double precision :: res, eps, f_binom_factors, ax, bx
!
	res = 0d0
	
	if ( gl1 .ge. gl2 ) then
		ax = gax
		bx = gbx
		l1 = gl1
		l2 = gl2
	else
		ax = gbx
		bx = gax
		l1 = gl2
		l2 = gl1
	end if
	
	! implemented up to f-orbitals to speed up the calculations
	select case ( J )
		case ( 0 )
			select case ( l1+l2 )
				case ( 0 ) 
					res = 1d0
					return
				case ( 1 )
					res = px-ax
					return
				case ( 2 )
					if ( l1 .gt. l2 ) then
						res = (px-ax)**2 + 1d0 / (2d0*gamma)
						return
					else
						res = (px-ax)*(px-bx) + 1d0 / (2d0*gamma)
						return
					end if
				case ( 3 )
					if ( l1 .eq. 3 ) then
						res = (px-ax)**3 + 3d0*(px-ax) / (2d0*gamma)
						return
					else
						res = (px-bx)*(px-ax)**2 + ( 2d0*(px-ax) + (px-bx) ) / (2d0*gamma)
						return
					end if
				case ( 4 )
					res = f_binom_factors(0, l1, l2, px-ax, px-bx) + f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/(4d0*gamma**2)
					return
				case ( 5 )
					res = f_binom_factors(0, l1, l2, px-ax, px-bx) + f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/4d0 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma**2
					return
				case ( 6 )
					res = f_binom_factors(0, l1, l2, px-ax, px-bx) + f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/4d0 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma**2 + 15d0/(8d0*gamma**3)
					return
				case ( 7 )
					res = f_binom_factors(0, l1, l2, px-ax, px-bx) + f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/4d0 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma**2 + 15d0/8d0 * f_binom_factors(6, l1, l2, px-ax, px-bx) / gamma**3 
					return
			end select
		case ( 1 )
			select case ( l1+l2 )
				case ( 1 )
					res = -(px-cx) 
					return
				case ( 2 )
					if ( l1 .gt. l2 ) then
						res = (cx-px) * 2d0 * (px-ax) - 1d0 / (2d0*gamma)
						return
					else
						res = (cx-px) * ( (px-ax) + (px-bx) ) - 1d0 / (2d0*gamma)
						return
					end if
				case ( 3 )
					if ( l1 .eq. 3 ) then
						res = (cx-px) * 3d0*(px-ax)**2 - 3d0*(px-ax)/(2d0*gamma) - 3d0/2d0*(px-cx)/gamma
						return
					else
						res = (cx-px) * ((px-ax)**2 + 2d0*(px-ax)*(px-bx)) - (2d0*(px-ax)+(px-bx))/(2d0*gamma) - 3d0/2d0*(px-cx)/gamma
						return
					end if
				case ( 4 )
					res = (cx-px) * f_binom_factors(1, l1, l2, px-ax, px-bx) - f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/2d0 * (cx-px) * f_binom_factors(3, l1, l2, px-ax, px-bx) / gamma - 3d0/2d0 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma**2
					return
				case ( 5 )
					res = (cx-px) * f_binom_factors(1, l1, l2, px-ax, px-bx) - f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/2d0 * (cx-px) * f_binom_factors(3, l1, l2, px-ax, px-bx) / gamma - 3d0/2d0 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma**2 - &
					      15d0/4d0 * (px-cx) * f_binom_factors(5, l1, l2, px-ax, px-bx) / gamma**2
					return
				case ( 6 )
					res = (cx-px) * f_binom_factors(1, l1, l2, px-ax, px-bx) - f_binom_factors(2, l1, l2, px-ax, px-bx) / (2d0*gamma) + &
					      3d0/2d0 * (cx-px) * f_binom_factors(3, l1, l2, px-ax, px-bx) / gamma - 3d0/2d0 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma**2 - &
					      15d0/4d0 * (px-cx) * f_binom_factors(5, l1, l2, px-ax, px-bx) / gamma**2 - 45d0/8d0 * f_binom_factors(6, l1, l2, px-ax, px-bx) / gamma**3
					return
			end select
		case ( 2 )
			select case ( l1+l2 )
				case ( 2 )
					res = +(px-cx)**2 
					return
				case ( 3 ) 
					if ( l1 .eq. 3 ) then
						res = (px-cx)**2 * 3d0 * (px-ax) + 3d0/2d0 * (px-cx) / gamma
						return
					else
						res = (px-cx)**2 * (2d0 * (px-ax) + (px-bx) ) + 3d0/2d0 * (px-cx) / gamma
						return
					end if
				case ( 4 ) 
					res = (px-cx)**2 * f_binom_factors(2, l1, l2, px-ax, px-bx) + 3d0/2d0 * (px-cx) * f_binom_factors(3, l1, l2, px-ax, px-bx) / gamma + &
						3d0*(px-cx)**2 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma + 3d0/(4d0 * gamma**2)
					return
				case ( 5 ) 
					res = (px-cx)**2 * f_binom_factors(2, l1, l2, px-ax, px-bx) + 3d0*(px-cx)*f_binom_factors(3, l1, l2, px-ax, px-bx)/(2d0*gamma) + &
						3d0 * f_binom_factors(4, l1, l2, px-ax, px-bx)/(4d0*gamma**2) + 3d0*(px-cx)**2* f_binom_factors(4, l1, l2, px-ax, px-bx)/gamma + &
						15d0 * (px-cx) * f_binom_factors(5, l1, l2, px-ax, px-bx)/(2d0*gamma**2)
					return
				case ( 6 ) 
					res = (px-cx)**2 * f_binom_factors(2, l1, l2, px-ax, px-bx) + 3d0*(px-cx)*f_binom_factors(3, l1, l2, px-ax, px-bx)/(2d0*gamma) + &
						3d0 * f_binom_factors(4, l1, l2, px-ax, px-bx)/(4d0*gamma**2) + 3d0*(px-cx)**2* f_binom_factors(4, l1, l2, px-ax, px-bx)/gamma + &
						15d0 * (px-cx) * f_binom_factors(5, l1, l2, px-ax, px-bx)/(2d0*gamma**2) + 45d0*f_binom_factors(6, l1, l2, px-ax, px-bx)/(8d0*gamma**3) + &
						45d0*(px-cx)**2*f_binom_factors(6, l1, l2, px-ax, px-bx)/(4d0*gamma**2)
					return
			end select
		case ( 3 )
			select case ( l1+l2 )
				case ( 3 )
					res = -(px-cx)**3 
					return
				case ( 4 )
					if ( l1 .eq. 4 ) then
						res = - (px-cx)**3 * 4d0 * (px-ax) - 3d0*(px-cx)**2 / gamma
						return
					else if ( l1 .eq. 3 ) then
						res = - (px-cx)**3 * ( 3d0*(px-ax) + (px-bx) ) - 3d0*(px-cx)**2 / gamma
						return
					else
						res = - (px-cx)**3 * ( 2d0*(px-ax) + 2d0*(px-bx) ) - 3d0*(px-cx)**2 / gamma
						return
					end if
				case ( 5 )
					res = - (px-cx)**3 * f_binom_factors(3, l1, l2, px-ax, px-bx) - 3d0*(px-cx)**2 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma - &
						15d0/4d0*(px-cx) * f_binom_factors(5, l1, l2, px-ax, px-bx) / gamma**2 - 5d0*(px-cx)**3 * f_binom_factors(5, l1, l2, px-ax, px-bx) / gamma
					return
				case ( 6 )
					res = - (px-cx)**3 * f_binom_factors(3, l1, l2, px-ax, px-bx) - 3d0*(px-cx)**2 * f_binom_factors(4, l1, l2, px-ax, px-bx) / gamma - &
						15d0/4d0*(px-cx) * f_binom_factors(5, l1, l2, px-ax, px-bx) / gamma**2 - 5d0*(px-cx)**3 * f_binom_factors(5, l1, l2, px-ax, px-bx) / gamma - &
						15d0*f_binom_factors(6, l1, l2, px-ax, px-bx)/(8d0*gamma**3) - 45d0*(px-cx)**2*f_binom_factors(6, l1, l2, px-ax, px-bx)/(2d0*gamma**2)
					return
			end select
		case ( 4 )
			select case ( l1+l2 )
				case ( 4 )
					res = +(px-cx)**4 
					return
				case ( 5 )
					res = +(px-cx)**4 * f_binom_factors(4, l1, l2, px-ax, px-bx) + 5d0*(px-cx)**3*f_binom_factors(5, l1, l2, px-ax, px-bx)/ gamma
					return
				case ( 6 )
					res = +(px-cx)**4 * f_binom_factors(4, l1, l2, px-ax, px-bx) + 5d0*(px-cx)**3*f_binom_factors(5, l1, l2, px-ax, px-bx)/ gamma + 45d0*(px-cx)**2*f_binom_factors(6, l1, l2, px-ax, px-bx)/(4d0*gamma**2) + &
					      15d0*(px-cx)**4*f_binom_factors(6, l1, l2, px-ax, px-bx)/(2d0*gamma)
					return
			end select
		case ( 5 )
			select case ( l1+l2 )
				case ( 5 )
					res = -(px-cx)**5 
					return
				case ( 6 )
					res = +(px-cx)**5 * f_binom_factors(5, l1, l2, px-ax, px-bx) - 15d0*(px-cx)**4*f_binom_factors(6, l1, l2, px-ax, px-bx)/(2d0*gamma)
					return
			end select
		case ( 6 )
			select case ( l1+l2 )
				case ( 6 )
					res = +(px-cx)**6
					return
			end select
	end select
	
	! general recursion formula, but very slow
	eps = 1 / (4d0*gamma)
	do i = 0, l1+l2
		do r = 0, i/2
			do u = 0, (i-2*r)/2
			
				if ( i-2*r-u .eq. J ) then
				
				      res = res + (-1d0)**i * f_binom_factors(i, l1, l2, px-ax, px-bx) * &
					  (-1d0)**u * n_fac(i) * (px-cx)**(i-2*r-2*u) * eps**(r+u) &
					  / ( n_fac(r) * n_fac(u) * n_fac(i-2*r-2*u) )
				
				end if
							
			end do
		end do
	end do
	
end function

! *************************************************************************
! * Binomial factors for integral evaluation, hard-coded up to f-shells
! *  see (3.5) in H. Taketa, S. Huzinaga, K. Ohata, J. Phys. Soc. Jpn. 1966, (21) 
! *************************************************************************
function H_binom_factors(L, l1, l2, a, b, gamma) result(res)
	implicit none
!
	integer, intent(in) :: L, l1, l2
	double precision, intent(in) :: a, b, gamma
	integer :: i, r, n_fac, vl1, vl2
	double precision :: res, va, vb, f_binom_factors
!
	! check for special cases
	if ( l1+l2 .eq. 0 ) then
		res = 1d0
		return
	else if ( L .eq. l1+l2 ) then
		res = 1d0/(4d0*gamma)**L
		return
	end if
	
	res = 0d0

	if ( l1 .ge. l2 ) then
		va = a
		vb = b
		vl1 = l1
		vl2 = l2
	else 
		va = b
		vb = a
		vl1 = l2
		vl2 = l1
	end if
	
	select case ( vl1+vl2 )
		case ( 1 )
			res = va
		case ( 2 )
			if ( vl1 .eq. 2 ) then
				if ( L .eq. 0 ) then
					res = va**2 + 1d0/(2d0*gamma)
				else
					res = va/(2d0*gamma)
				end if
			else
				if ( L .eq. 0 ) then
					res = va*vb + 1d0/(2d0*gamma)
				else
					res = (va+vb)/(4d0*gamma)
				end if
			end if
		case ( 3 )
			if ( vl1 .eq. 3 ) then
				select case ( L )
					case ( 0 )
						res = va**3 + (3d0*va)/(2d0*gamma)
					case ( 1 )
						res = (3d0*va**2)/(4d0*gamma) + 3d0/(8d0*gamma**2)
					case ( 2 )
						res = (3d0*va)/(16d0*gamma**2)
				end select
			else
				select case ( L )
					case ( 0 )
						res = va**2*vb + (2d0*va+vb)/(2d0*gamma)
					case ( 1 ) 
						res = (va**2 + 2d0*va*vb)/(4d0*gamma) + 3d0/(8d0*gamma**2)
					case ( 2 )
						res = (2d0*va+vb)/(16d0*gamma**2)
				end select
			end if
		case ( 4 )
			if ( vl1 .eq. 4 ) then
				select case ( L )
					case ( 0 )
						res = va**4 + 3d0/(4d0*gamma**2) + (3d0*va**2)/gamma
					case ( 1 )
						res = va**3/gamma + (3d0*va)/(2d0*gamma**2)
					case ( 2 )
						res = 3d0/(16d0*gamma**3) + (3d0*va**2)/(8d0*gamma**2)
					case ( 3 )
						res = va/(16d0*gamma**3)
				end select
			else if ( vl1 .eq. 3 ) then
				select case ( L )
					case ( 0 )
						res = va**3*vb + (3d0*va**2 + 3d0*va*vb)/(2d0*gamma) + 3d0/(4d0*gamma**2)
					case ( 1 )
						res = (3d0*(3d0*va + vb))/(8d0*gamma**2) + (va**3 + 3d0*va**2*vb)/(4d0*gamma)
					case ( 2 )
						res = (3d0*va**2 + 3d0*va*vb)/(16d0*gamma**2) + 3d0/(16d0*gamma**3)
					case ( 3 )
						res = (3d0*va + vb)/(64d0*gamma**3)
				end select
			else
				select case ( L )
					case ( 0 )
						res = 3d0/(4d0*gamma**2) + va**2*vb**2 + (va**2 + 4d0*va*vb + vb**2)/(2d0*gamma)
					case ( 1 )
						res = (3d0*(2d0*va + 2d0*vb))/(8d0*gamma**2) + (2d0*va**2*vb + 2d0*va*vb**2)/(4d0*gamma)
					case ( 2 )
						res = 3d0/(16d0*gamma**3) + (va**2 + 4d0*va*vb + vb**2)/(16d0*gamma**2)
					case ( 3 )
						res = (2d0*va + 2d0*vb)/(64d0*gamma**3)
				end select
			end if
		case ( 5 )
			if ( vl1 .eq. 5 ) then
				select case ( L )
					case ( 0 )
						res = (15d0*va)/(4d0*gamma**2) + (5d0*va**3)/gamma + va**5
					case ( 1 )
						res = 15d0/(16d0*gamma**3) + (15d0*va**2)/(4d0*gamma**2) + (5d0*va**4)/(4d0*gamma)
					case ( 2 )
						res = (15d0*va)/(16d0*gamma**3) + (5d0*va**3)/(8d0*gamma**2)
					case ( 3 )
						res = 5d0/(64d0*gamma**4) + (5d0*va**2)/(32d0*gamma**3)
					case ( 4 )
						res = (5d0*va)/(256d0*gamma**4)
				end select
			else if ( vl1 .eq. 4 ) then
				select case ( L )
					case ( 0 )
						res = va**4*vb + (3d0*(4d0*va + vb))/(4d0*gamma**2) + (4d0*va**3 + 6d0*va**2*vb)/(2d0*gamma)
					case ( 1 )
						res = 15d0/(16d0*gamma**3) + (3d0*(6d0*va**2 + 4d0*va*vb))/(8d0*gamma**2) + (va**4 + 4d0*va**3*vb)/(4d0*gamma)
					case ( 2 )
						res = (3d0*(4d0*va + vb))/(16d0*gamma**3) + (4d0*va**3 + 6d0*va**2*vb)/(16d0*gamma**2)
					case ( 3 )
						res = 5d0/(64d0*gamma**4) + (6d0*va**2 + 4d0*va*vb)/(64d0*gamma**3)
					case ( 4 )
						res = (4d0*va + vb)/(256d0*gamma**4)
				end select
			else
				select case ( L )
					case ( 0 )
						res = va**3*vb**2 + (3d0*(3d0*va + 2d0*vb))/(4d0*gamma**2) + (va**3 + 6d0*va**2*vb + 3d0*va*vb**2)/(2d0*gamma)
					case ( 1 )
						res = 15d0/(16d0*gamma**3) + (3d0*(3d0*va**2 + 6d0*va*vb + vb**2))/(8d0*gamma**2) + (2d0*va**3*vb + 3d0*va**2*vb**2)/(4d0*gamma)
					case ( 2 )
						res = (3d0*(3d0*va + 2d0*vb))/(16d0*gamma**3) + (va**3 + 6d0*va**2*vb + 3d0*va*vb**2)/(16d0*gamma**2)
					case ( 3 )
						res = 5d0/(64d0*gamma**4) + (3d0*va**2 + 6d0*va*vb + vb**2)/(64d0*gamma**3)
					case ( 4 )
						res = (3d0*va + 2d0*vb)/(256d0*gamma**4)
				end select
			end if
		case ( 6 )
			if ( vl1 .eq. 6 ) then
				select case ( L )
					case ( 0 )
						res = 15d0/(8d0*gamma**3) + (45d0*va**2)/(4d0*gamma**2) + (15d0*va**4)/(2d0*gamma) + va**6
					case ( 1 )
						res = (45d0*va)/(8d0*gamma**3) + (15d0*va**3)/(2d0*gamma**2) + (3d0*va**5)/(2d0*gamma)
					case ( 2 )
						res = 45d0/(64d0*gamma**4) + (45d0*va**2)/(16d0*gamma**3) + (15d0*va**4)/(16d0*gamma**2)
					case ( 3 )
						res = (15d0*va)/(32d0*gamma**4) + (5d0*va**3)/(16d0*gamma**3)
					case ( 4 )
						res = 15d0/(512d0*gamma**5) + (15d0*va**2)/(256d0*gamma**4)
					case ( 5 )
						res = (3d0*va)/(512d0*gamma**5)
				end select
			else if ( vl1 .eq. 5 ) then
				select case ( L )
					case ( 0 )
						res = 15d0/(8d0*gamma**3) + va**5*vb + (3d0*(10d0*va**2 + 5d0*va*vb))/(4d0*gamma**2) + (5d0*va**4 + 10d0*va**3*vb)/(2d0*gamma)
					case ( 1 )
						res = (15d0*(5d0*va + vb))/(16d0*gamma**3) + (3d0*(10d0*va**3 + 10d0*va**2*vb))/(8d0*gamma**2) + (va**5 + 5d0*va**4*vb)/(4d0*gamma)
					case ( 2 )
						res = 45d0/(64d0*gamma**4) + (3d0*(10d0*va**2 + 5d0*va*vb))/(16d0*gamma**3) + (5d0*va**4 + 10d0*va**3*vb)/(16d0*gamma**2)
					case ( 3 )
						res = (5d0*(5d0*va + vb))/(64d0*gamma**4) + (10d0*va**3 + 10d0*va**2*vb)/(64d0*gamma**3)
					case ( 4 )
						res = 15d0/(512d0*gamma**5) + (10d0*va**2 + 5d0*va*vb)/(256d0*gamma**4)
					case ( 5 )
						res = (5d0*va + vb)/(1024d0*gamma**5)
				end select
			else if ( vl1 .eq. 4 ) then
				select case ( L )
					case ( 0 )
						res = 15d0/(8d0*gamma**3) + va**4*vb**2 + (3d0*(6d0*va**2 + 8d0*va*vb + vb**2))/(4d0*gamma**2) + (va**4 + 8d0*va**3*vb + 6d0*va**2*vb**2)/(2d0*gamma)
					case ( 1 )
						res = (15d0*(4d0*va + 2d0*vb))/(16d0*gamma**3) + (3d0*(4d0*va**3 + 12d0*va**2*vb + 4d0*va*vb**2))/(8d0*gamma**2) + (2d0*va**4d0*vb + 4d0*va**3*vb**2)/(4d0*gamma)
					case ( 2 )
						res = 45d0/(64d0*gamma**4) + (3d0*(6d0*va**2 + 8d0*va*vb + vb**2))/(16d0*gamma**3) + (va**4 + 8d0*va**3*vb + 6d0*va**2*vb**2)/(16d0*gamma**2)
					case ( 3 )
						res = (5d0*(4d0*va + 2d0*vb))/(64d0*gamma**4) + (4d0*va**3 + 12d0*va**2*vb + 4d0*va*vb**2)/(64d0*gamma**3)
					case ( 4 )
						res = 15d0/(512d0*gamma**5) + (6d0*va**2 + 8d0*va*vb + vb**2)/(256d0*gamma**4)
					case ( 5 )
						res = (4d0*va + 2d0*vb)/(1024d0*gamma**5)
				end select
			else
				select case ( L )
					case ( 0 )
						res = 15d0/(8d0*gamma**3) + va**3*vb**3 + (3d0*(3d0*va**2 + 9d0*va*vb + 3*vb**2))/(4d0*gamma**2) + (3d0*va**3*vb + 9d0*va**2*vb**2 + 3d0*va*vb**3)/(2d0*gamma)
					case ( 1 )
						res = (15d0*(3d0*va + 3d0*vb))/(16d0*gamma**3) + (3d0*(va**3 + 9d0*va**2*vb + 9d0*va*vb**2 + vb**3))/(8d0*gamma**2) + (3d0*va**3*vb**2 + 3d0*va**2*vb**3)/(4d0*gamma)
					case ( 2 )
						res = 45d0/(64d0*gamma**4) + (3d0*(3d0*va**2 + 9d0*va*vb + 3d0*vb**2))/(16d0*gamma**3) + (3d0*va**3*vb + 9d0*va**2*vb**2 + 3d0*va*vb**3)/(16d0*gamma**2)
					case ( 3 )
						res = (5d0*(3d0*va + 3d0*vb))/(64d0*gamma**4) + (va**3 + 9d0*va**2*vb + 9d0*va*vb**2 + vb**3)/(64d0*gamma**3)
					case ( 4 )
						res = 15d0/(512d0*gamma**5) + (3d0*va**2 + 9d0*va*vb + 3d0*vb**2)/(256d0*gamma**4)
					case ( 5 )
						res = (3d0*va + 3d0*vb)/(1024d0*gamma**5)
				end select
			end if
		case default
			! (4*gamma)**(i-r) instead of (4*gamma)**(i-2*r)
			do i = 0, vl1+vl2
				do r = 0, i/2
					if ( i - 2*r .eq. L ) then
						res = res + n_fac(i) * f_binom_factors(i, vl1, vl2, va, vb) / &
							( dble(n_fac(r)) * dble(n_fac(i-2*r)) * (4*gamma)**(i-r) )
					end if
				end do
			end do
			
	end select

end function

subroutine build_ortho_matrix
	use data_mo
	implicit none
!
	double precision, dimension(:,:), allocatable :: S_CGTOEigV, S_TMPA, S_TMPB
	double precision, dimension(:), allocatable :: S_CGTOEigVal
	double precision, dimension(:), allocatable :: S_CGTOWork
	integer :: i, j, iost, lapack_status
!
	write(*,'(5X,A)') "> DIAGONALIZING OVERLAP MATRIX..."
	write(*,*)
	
	if ( ( iost .eq. 0 ) .and. ( allocated(S_CGTOEigV) .eqv. .FALSE. ) ) then
		allocate(S_CGTOEigV(NumCGTO, NumCGTO), stat=iost)
	end if
	if ( ( iost .eq. 0 ) .and. ( allocated(S_TMPA) .eqv. .FALSE. ) ) then
		allocate(S_TMPA(NumCGTO, NumCGTO), stat=iost)
	end if
	if ( ( iost .eq. 0 ) .and. ( allocated(S_TMPB) .eqv. .FALSE. ) ) then
		allocate(S_TMPB(NumCGTO, NumCGTO), stat=iost)
	end if
	
	if ( iost .ne. 0 ) then
		write(*,*) "Error! Can't allocate memory for matrices!"
		write(*,*)
		stop
	end if
	
	S_CGTOEigV = S_CGTO
	
	allocate(S_CGTOEigVal(NumCGTO))
	allocate(S_CGTOWork(max(1,3*NumCGTO-1)))
	
	call DSYEV('V', 'L', NumCGTO, S_CGTOEigV, NumCGTO, S_CGTOEigVal, S_CGTOWork, size(S_CGTOWork), lapack_status)
	if ( lapack_status .ne. 0 ) then
		write(*,*) "Error! LAPACK returned error code: ", lapack_status
		stop
	end if
	
	S_TMPA = 0d0
	S_TMPB = 0d0
	
	do i = 1, NumCGTO
		S_TMPA(i,i) = 1d0 / dsqrt(S_CGTOEigVal(i))
		do j = 1, i
			S_TMPB(i, j) = S_CGTOEigV(j, i)
			S_TMPB(j, i) = S_CGTOEigV(i, j)
		end do
	end do
	
	O_CGTO = matmul(matmul(S_CGTOEigV, S_TMPA), S_TMPB)
	O_CGTOT = transpose(O_CGTO)
	
	do i = 1, NumCGTO
		S_TMPA(i,i) = dsqrt(S_CGTOEigVal(i))
	end do
	
	O_CGTO12 = matmul(matmul(S_CGTOEigV, S_TMPA), S_TMPB)
	
	deallocate(S_TMPB)
	deallocate(S_TMPA)
	deallocate(S_CGTOEigV)
	deallocate(S_CGTOEigVal)
	deallocate(S_CGTOWork)

end subroutine

subroutine print_eri_statistics
	use global_c
	use data_mo
	implicit none
!
	integer :: i, j, k, l, ij, kl
	integer (8) :: N_PGTO_ERIs, N_CGTO_ERIs, N_CGTO_SymERIs
!	
	N_PGTO_ERIs = NumPGTO**4
	N_CGTO_ERIs = NumCGTO**4
	
	write(*,*)
	write(*,*) "PRIMITIVE TWO-ELECTRON REPULSION INTEGRAL STATISTICS: "
	write(*,*)
	write(*,'(A35,I16)') "primitive integrals: ", N_PGTO_ERIs 
	write(*,'(A35,I16)') "contracted integrals: ", N_CGTO_ERIs
	write(*,*)
	
	N_CGTO_SymERIs = 0
	
	!$OMP PARALLEL DO REDUCTION(+:N_CGTO_SymERIs) PRIVATE(i, j, ij, k, l, kl) SCHEDULE(dynamic)
	do i = 1, NumCGTO
		do j = 1, i
			ij = (i-1)*i/2 + j-1
			do k = 1, NumCGTO
				do l = 1, k
					kl = (k-1)*k/2 + l-1
					if ( ij .ge. kl ) then
						N_CGTO_SymERIs = N_CGTO_SymERIs + 1
					end if
				end do
			end do
		end do
	end do
	!$OMP END PARALLEL DO
	
	write(*,'(A35,I16)') "unique contracted integrals: ", N_CGTO_SymERIs
	write(*,*)
	
end subroutine

subroutine calc_mulliken_charges(verbose)
	use data_mo
	use global_c
	implicit none
!
	logical, intent(in) :: verbose
	double precision :: TotalDens, SumMC, fCharge
	integer :: i, j
	character(len=5) :: ghost
!
	if ( verbose .eqv. .TRUE. ) then
		if ( NCustomCharges .eq. 0 ) then
			write(*,'(X,A/)') "> ATOMS AND CALCULATED MULLIKEN CHARGES:"
		else
			write(*,'(X,A/)') "> ATOMS AND CORRESPONDING CUSTOM CHARGES:"
		end if
		call write_rline(75)
		write(*,'(6X,A4,A5,4A16)') "#", "type", "x / angstrom", "y / angstrom", "z / angstrom", "charge"
		call write_rline(75)
	end if
	
	SumMC = 0d0

	do i = 1, NA
		TotalDens = 0d0
		do j = 1, NumCGTO
			if ( CGTOs(j)%atomid .eq. i ) then
				TotalDens = TotalDens + DSMat(j, j)
			end if
		end do
		
		if ( OUHF .eqv. .FALSE. ) then
			Atoms(i)%mulliken = Atoms(i)%charge - TotalDens * 2d0
		else
			Atoms(i)%mulliken = Atoms(i)%charge - TotalDens * 1d0
		end if
		
		if ( NCustomCharges .eq. 0 ) then
			fCharge = Atoms(i)%mulliken
		else
			fCharge = Atoms(i)%ccharge
		end if
		
		SumMC = SumMC + fCharge 
		
		if ( verbose .eqv. .TRUE. ) then
			if ( Atoms(i)%ghost .eqv. .TRUE. ) then
				ghost = "GHOST"
			else
				ghost = ""
			end if
			write(*,'(X,A5,I4,X,A3,X,3F16.8,F16.6)') ghost, Atoms(i)%id, Atoms(i)%str, Atoms(i)%x * au2pm, Atoms(i)%y * au2pm, Atoms(i)%z * au2pm, fCharge
		end if
	end do
	
	if ( verbose .eqv. .TRUE. ) then
		call write_rline(75)
		write(*,'(48X,A,F15.6)') "sum of charges: ", SumMC
		call write_rline(75)
		write(*,*)
		write(*,*)
	end if

end subroutine

subroutine check_for_ecp
	use data_mo
	use global_c
	implicit none
!
	integer :: i
!
	if ( OExpert .eqv. .FALSE. ) then
		do i = 1, NA
			if ( dabs(Atoms(i)%mulliken) .gt. 8.0d0 ) then
				write(*,*) "*****************************************************************************"
				write(*,*) "*                                  ERROR!                                   *"
				write(*,*) "*****************************************************************************"
				write(*,*) 
				write(*,*) "It seems that you have applied an effective core potential (ECP)!"
				write(*,*) "The results you will get might be wrong since ECPs are not supported."
				write(*,*) 
				write(*,*) "Either you utilize full electron basis sets to your calculations or"
				write(*,*) "use the 'EXPERT' keyword in your input file to ignore this error."
				write(*,*)
				write(*,*)
				stop
			end if
		end do
	end if

end subroutine

subroutine calc_norm_wf(verbose)
	use data_mo
	use global_c
	implicit none
!
	logical, intent(in) :: verbose
	double precision :: TotalDens, TotalNorm, SNorm
	integer :: i, j, k
!	
	TotalDens = 0d0
	TotalNorm = 0d0
	
	if ( verbose .eqv. .TRUE. ) then
		write(*,*)
		write(*,*) "CALCULATING NORM OF MOs..."		
		write(*,'(39X,A)') "--------------------------------------"
		write(*,'(39X,A8,A20,A10)') "MO #", "NORM", "OCCUP."
		write(*,'(39X,A)') "--------------------------------------"	
	end if
	
	do k = 1, NMO
	
		SNorm = 0d0
		
		do i = 1, NumCGTO
			do j = 1, NumCGTO
				SNorm = SNorm + S_CGTO(i, j) * MOs(k)%mocoeff(j) * MOs(k)%mocoeff(i)
			end do
		end do
	
		TotalNorm = TotalNorm + SNorm
		TotalDens = TotalDens + SNorm * MOs(k)%occup
		
		if ( verbose .eqv. .TRUE. ) then
			write(*,'(39X,I8,F20.13,F10.4)') k, SNorm, MOs(k)%occup
		end if
	
	end do
	
	TotDensMO = TotalDens
	TotNormWF = TotalNorm / dble(NMO)
	
	if ( verbose .eqv. .TRUE. ) then
		write(*,'(39X,A)') "--------------------------------------"	
		write(*,*)
		write(*,'(A47,F20.13)') "total norm of WF:", TotNormWF 
		write(*,'(A47,F20.13)') "total density from MOs:", TotDensMO 
		write(*,*)
	end if

end subroutine

subroutine allocate_matrices
	use data_mo
	use global_c
	implicit none
!
	integer :: iost
	integer :: mem
	
	iost = 0
	
	write(*,'(5X,A)') "> ALLOCATING MEMORY ..."
	write(*,*)
	
	! overlap matrix
	if ( ( iost .eq. 0 ) .and. ( allocated(S_PGTO) .eqv. .FALSE. ) ) then
		allocate(S_PGTO(NumPGTO, NumPGTO), stat=iost)
	end if	
	if ( ( iost .eq. 0 ) .and. ( allocated(S_CGTO) .eqv. .FALSE. ) ) then
		allocate(S_CGTO(NumCGTO, NumCGTO), stat=iost)
	end if
	
	! orthogonalization matrix
	if ( ( iost .eq. 0 ) .and. ( allocated(O_CGTO) .eqv. .FALSE. ) ) then
		allocate(O_CGTO(NumCGTO, NumCGTO), stat=iost)
	end if
	
	! S^(1/2) matrix
	if ( ( iost .eq. 0 ) .and. ( allocated(O_CGTO12) .eqv. .FALSE. ) ) then
		allocate(O_CGTO12(NumCGTO, NumCGTO), stat=iost)
	end if
	
	! transposed orthogonalization matrix
	if ( ( iost .eq. 0 ) .and. ( allocated(O_CGTOT) .eqv. .FALSE. ) ) then
		allocate(O_CGTOT(NumCGTO, NumCGTO), stat=iost)
	end if
	
	! D*S matrix
	if ( ( iost .eq. 0 ) .and. ( allocated(DSMat) .eqv. .FALSE. ) ) then
		allocate(DSMat(NumCGTO, NumCGTO), stat=iost)
	end if
	
	! density matrix
	if ( ( iost .eq. 0 ) .and. ( allocated(D_CGTO) .eqv. .FALSE. ) ) then
		allocate(D_CGTO(NumCGTO, NumCGTO), stat=iost)
	end if
	
	if ( iost .ne. 0 ) then
		write(*,*) "Error! Can't allocate memory for matrices!"
		write(*,*)
		stop
	end if
	
	mem = 	sizeof(S_PGTO) + sizeof(S_CGTO) + &
		sizeof(O_CGTO) + sizeof(O_CGTOT) + &
		sizeof(O_CGTO12) + &
		sizeof(DSMat) + sizeof(D_CGTO) 
		
	write(*,'(A35,F14.2,A)') "allocated memory: ", mem/1024d0/1024d0, " MB"
	write(*,*)
	write(*,*)
	
end subroutine

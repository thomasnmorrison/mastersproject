!!!! zeta.f90 !!!!
!!!! Compute the zeta parameter
!!!! 

! to do: Clean out redundant variables/subroutines
! to do: output equation of state
! to do: add dzeta_lat_part to outputs																	- done
! to do: declare variables for partial epsilon 													- done
! to do: add calculation of partial epsilon															- done						
! to do: fix partial zeta integration to account for partial epsilon		- done
! to do: make output for partial zeta, partial epsilon, etc.						- done
! to do: declare variables for linear/non-linear parts of zeta split		- done
! to do: calculate linear/non-linear parts of dzeta splits							- done
! to do: calculate linear/non-linear parts of zeta splits								- done
! to do: output linear/non-linear parts of zeta splits									- done

! to do: make subroutine that decomposes the fields into various components ie longitudinal and transverse
! to do: make subroutine to calculate zeta split phi - chi
! to do: make subroutine to calculate zeta split longitudinal - transverse
! to do: make subroutine to calculate zeta split slow roll logintudinal - transverse

module zeta_mod

!	use, intrinsic :: iso_c_binding
!  include 'fftw3.f03'
!	use fftw3
!	use params
#include "macros.h"
	use hamiltonian
	use sample_sites_mod
	use moments_mod

	real(dl) :: zeta = 0.0_dl
	real(dl), dimension(2) :: dzeta = (/0.0_dl,0.0_dl/)
	real(dl), dimension(IRANGE) :: zeta_lat
	real(dl), dimension(2,IRANGE) :: dzeta_lat !dzeta/dtau
	!real(dl), dimension(2,IRANGE) :: dzeta_lat_temp	! 1st component is numerator of dzeta, 2nd component is denominator of dzeta
	real(dl), dimension(nfld,IRANGE) :: epsilon_part	! components rho*x_A*epsilon_A
	real(dl), dimension(2*nfld,IRANGE) :: zeta_part			! zeta component calculated from individual fields (and linear/non-linear)
	real(dl), dimension(2,2*nfld,IRANGE) :: dzeta_part	! dzeta/dtau component calculated from individual fields
	
	! additional pointers needed for calculating dzeta
	real(C_DOUBLE), pointer :: ggrad1(:,:,:)
!	real(C_DOUBLE), pointer :: ggrad2(:,:,:)
!	real(C_DOUBLE), pointer :: ggrad3(:,:,:)
!	real(C_DOUBLE), pointer :: ggrad4(:,:,:)
	real(C_DOUBLE), pointer :: zf1(:,:,:)
	real(C_DOUBLE), pointer :: zf2(:,:,:)
!	real(C_DOUBLE), pointer :: zf3(:,:,:)!new
!	real(C_DOUBLE), pointer :: zf4(:,:,:)!new
!	real(C_DOUBLE), pointer :: zfp1(:,:,:)
!	real(C_DOUBLE), pointer :: zfp2(:,:,:)
	real(C_DOUBLE), pointer :: zlap1(:,:,:)
	real(C_DOUBLE), pointer :: zlap2(:,:,:)
	complex(C_DOUBLE_COMPLEX), pointer :: Fk3(:,:,:)

contains
	!!!! Step zeta integration !!!!
!	subroutine zeta_step(zeta, dzeta, dt)
!		real(dl), dimension(2) :: dzeta
!		real(dl) :: dt
!		
!		zeta = zeta + 0.5_dl*sum(dzeta)*dt	!!
!	end subroutine zeta_step


	subroutine zeta_init()
		zeta = 0.0_dl
		dzeta(:) = 0.0_dl
		zeta_lat(IRANGE) = 0.0_dl
		dzeta_lat(:,IRANGE) = 0.0_dl
	end subroutine zeta_init

	!!!! Step zeta integration using trapazoid rule
	subroutine zeta_step(dt, nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt
		integer, dimension(1:3), intent(in) :: nsize
		integer :: i
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		!call get_dzeta_test(nsize, dk, Fk1, Fk2, Fk3, planf, planb)			! call for dzeta calc no smoothing
		!call get_dzeta_smooth(nsize, dk, Fk1, Fk2, Fk3, planf, planb)		! call for dzeta calc from smoothed T
		call get_dzeta_part(nsize, dk, Fk1, Fk2, Fk3, planf, planb)				! call for dzeta calc with components
		zeta = zeta + 0.5_dl*sum(dzeta)*dt
		sample_zeta(:) = sample_zeta(:) + 0.5_dl*(sample_dzeta(1,:)+sample_dzeta(2,:))*dt
		!sample_zeta_f(:,:) = sample_zeta_f(:,:) + 0.5_dl*(sample_dzeta_f(1,:,:)+sample_dzeta_f(2,:,:))*dt
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + 0.5_dl*(dzeta_lat(1,IRANGE)+dzeta_lat(2,IRANGE))*dt
		do i=1,nfld
			!zeta_part(i,IRANGE) = zeta_part(i,IRANGE) + 0.5_dl*(dzeta_part(1,i,IRANGE) &
			!											+ epsilon_part(i,IRANGE)*dzeta_part(2,i,IRANGE)/(sum(epsilon_part(:,IRANGE) ,dim=1)))*dt
			zeta_part(2*i-1,IRANGE) = zeta_part(2*i-1,IRANGE) + 0.5_dl*(dzeta_part(1,2*i-1,IRANGE) &
														+ epsilon_part(i,IRANGE)*dzeta_part(2,2*i-1,IRANGE)/(sum(epsilon_part(:,IRANGE) ,dim=1)))*dt
			zeta_part(2*i,IRANGE) = zeta_part(2*i,IRANGE) + 0.5_dl*(dzeta_part(1,2*i,IRANGE) &
														+ epsilon_part(i,IRANGE)*dzeta_part(2,2*i,IRANGE)/(sum(epsilon_part(:,IRANGE) ,dim=1)))*dt
		enddo
	end subroutine zeta_step

!	subroutine set_fld_temp()
!		if (nfld == 2) then
!			zf1 = fld(1,IRANGE)
!			zf2 = fld(2,IRANGE)
!			zfp1 = fldp(1,IRANGE)
!			zfp2 = fldp(2,IRANGE)
!		elseif (nfld ==1) then
!			zf1 = fld(1,IRANGE)
!			zf2 = 0.0_dl*fld(1,IRANGE)
			!zf2 = fld(1,IRANGE)!THIS IS WRONG ONLY FOR TESTING
!			zfp1 = fldp(1,IRANGE)
!			zfp2 = 0.0_dl*fldp(1,IRANGE)
			!zfp2 = fldp(1,IRANGE)!THIS IS WRONG ONLY FOR TESTING
!		endif
!	end subroutine set_fld_temp

	! subroutint to set temp fields in get_dzeta_test
!#ifdef GHOST
!	subroutine set_fld_temp_test(ind)
!		integer :: ind
		
!		if (ind <= nfld) then
!			zf1 = fld(ind,IRANGE)
!			zfp1 = fldp(ind,IRANGE)
!		elseif (ind > nfld) then
!			zf1 = ghst(ind-nfld,IRANGE)
!			zfp1 = ghstp(ind-nfld,IRANGE)
!	end subroutine set_fld_temp_test
!#endif

	!!!! Updates new value of dzeta into dzeta(2), stores old value of dzeta into dzeta(1)
!	subroutine get_dzeta(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
!		integer, dimension(1:3), intent(in) :: nsize
!		real*8 :: dk	
!    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
!		type(C_PTR) :: planf, planb

!		real(dl) :: hub
!		integer :: n1, n2, n3
!		integer :: i,j

!		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
!		!hub = -ysclp/(6.0_dl*yscl**2)/dble(n1)/dble(n2)/dble(n3) !!!!Double check that this should not be dividing number of lattice sites
!		hub = -ysclp/(6.0_dl*yscl**2)

!		dzeta(1) = dzeta(2)	!store old value
!		dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
!		sample_dzeta(1,:) = sample_dzeta(2,:)
		!sample_dzeta_f(1,:,:) = sample_dzeta_f(2,:,:)

!		call set_fld_temp()
!		zlap1 = zfp1
!		zlap2 = zf1
!		call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
		!if (nfld==2) then
!			zlap1 = zfp2
!			zlap2 = zf2
!			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad2,Fk1,Fk2,Fk3, planf, planb)
		!endif
!		zlap1 = zf1
!		call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)
		!if (nfld==2) then
!			zlap2 = zf2
!			call laplacian_spectral(n1,n2,n3,zlap2,Fk1,dk,planf,planb)
		!endif
!		zf3 = zf1
!		call graddotgrad(nsize,dk,zf1,zf3,ggrad3,Fk1,Fk2,Fk3, planf, planb)
		!if (nfld==2) then
!			zf4 = zf2
!			call graddotgrad(nsize,dk,zf2,zf4,ggrad4,Fk1,Fk2,Fk3, planf, planb)
		!endif		
		!!!! Use graddotgrad instead of gradient_squared_3d_spectral, to be on the same side should use a field copy for one input
		!call gradient_squared_3d_spectral([n1,n2,n3],zf1,Fk1,Fk2,ggrad3,dk,planf,planb)
		!call gradient_squared_3d_spectral([n1,n2,n3],zf2,Fk1,Fk2,ggrad4,dk,planf,planb)

		!zf1 = (yscl*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2) + hub*yscl**4*(ggrad3 + ggrad4)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))!!!Included extra term
		!if (nfld==2) then
!			zf1 = (yscl**2*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))!!!Without extra term and squaring yscl to account for conformal time
			! First order piece only
			!zf1 = (yscl**2*(sum(fldp(1,IRANGE))/nvol*zlap1)) / (3.0_dl*sum(fldp(1,IRANGE)**2)/nvol + yscl**4*sum(ggrad3(IRANGE))/nvol)
		!endif
		!if (nfld==1) then
		!	zf1 = (yscl**2*(ggrad1 + zfp1*zlap1)) / (3.0_dl*zfp1**2 + yscl**4*(ggrad3))!!!Without extra term and squaring yscl to account for conformal time
		!endif
		!zf1 = (yscl*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2)

		
!		do i=1,n_sample
			!sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,1),sample_site(i,1))/dble(n1)/dble(n2)/dble(n3)!!!!!!!!!!!!SHOULD NOT DIVIDE BY LATTICE SITES!!!!!!!!!!
			!sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,1),sample_site(i,1))!looks like I'm choosing the wrong site
!			sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,2),sample_site(i,3))
!		enddo
!		dzeta(2) = sum(zf1(IRANGE))/dble(n1)/dble(n2)/dble(n3)
!		dzeta_lat(2,IRANGE) = zf1(IRANGE)

	!testing
!		do i=1,nfld
!			if (i==1) then
!				zf2 = (yscl**2*(ggrad1 + zfp1*zlap1)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))
!			elseif (i==2) then
!				zf2 = (yscl**2*(ggrad2 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))
!			endif
!			do j=1,n_sample
				!sample_dzeta_f(2,i,j) = zf2(sample_site(j,1),sample_site(j,2),sample_site(j,3))
!			enddo	
!		enddo
		!print*, ggrad1(6,6,6)
		!print*, zfp1(6,6,6)
		!print*, zlap1(6,6,6)
		!print*, ggrad3(6,6,6)
		!print*, zf1(6,6,6)
!	end subroutine get_dzeta

! subroutine to get dzeta for an arbitrary number of fields
! to do: Allow for smoothing of numerator and denominator
! to do: 
	subroutine get_dzeta_test(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		dzeta(1) = dzeta(2)	!store old value
		dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
		sample_dzeta(1,:) = sample_dzeta(2,:)

		zf1(IRANGE) = 0._dl
		zf2(IRANGE) = 0._dl
		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		

			zf1 = zf1 + ggrad1 + fldp(i,IRANGE)*zlap1
		enddo
		! Get denominator of dzeta
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			
			zf2 = zf2 + 3.0_dl*fldp(i,IRANGE)**2 + yscl**4*ggrad1
		enddo			
		! Calculate dzeta
		dzeta_lat(2,IRANGE) = yscl**2 * zf1 / zf2
		dzeta(2) = sum(dzeta_lat(2, IRANGE))/dble(n1)/dble(n2)/dble(n3)

		!do i=1,n_sample
		!	sample_dzeta(2,i) = dzeta_lat(2,sample_site(i,1),sample_site(i,2),sample_site(i,3))
		!enddo

	end subroutine get_dzeta_test

! Subroutine to calculate dzeta and the components of dzeta for individual species
! Does not include smoothing
! to do: store old value of dzeta_lat_comp
! to do: calulate new value of dzeta_lat_comp
! to do: double check dzeta partial and epsilon partial calculations 
	subroutine get_dzeta_part(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		dzeta(1) = dzeta(2)	!store old value
		dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
		!dzeta_lat_comp(1,:,IRANGE) = dzeta_lat_comp(2,:,IRANGE)
		! Store old value of x_A*epsilon_A*dzeta_A/epsilon
		do i=1,nfld
			dzeta_part(1,2*i-1,IRANGE) = epsilon_part(i,IRANGE)*dzeta_part(2,2*i-1,IRANGE)/(sum(epsilon_part(:,IRANGE) ,dim=1))
			dzeta_part(1,2*i,IRANGE) = epsilon_part(i,IRANGE)*dzeta_part(2,2*i,IRANGE)/(sum(epsilon_part(:,IRANGE) ,dim=1))
		enddo
		sample_dzeta(1,:) = sample_dzeta(2,:)

		zf1(IRANGE) = 0._dl
		zf2(IRANGE) = 0._dl
		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		
			
			zf1 = zf1 + ggrad1 + fldp(i,IRANGE)*zlap1
			epsilon_part(i,IRANGE) = 1.5_dl*(fldp(i,IRANGE)**2)/yscl**6
			dzeta_part(2,2*i-1,IRANGE) =  fldp(i,IRANGE)*zlap1
			dzeta_part(2,2*i,IRANGE) = ggrad1
		enddo
		! Get denominator of dzeta
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			
			zf2 = zf2 + 3.0_dl*fldp(i,IRANGE)**2 + yscl**4*ggrad1
			epsilon_part(i,IRANGE) = epsilon_part(i,IRANGE) + 0.5_dl*(ggrad1)/yscl**2
			dzeta_part(2,2*i-1,IRANGE) = yscl**2 *dzeta_part(2,2*i-1,IRANGE) / (3.0_dl*fldp(i,IRANGE)**2 + yscl**4*ggrad1)
			dzeta_part(2,2*i,IRANGE) = yscl**2 *dzeta_part(2,2*i,IRANGE) / (3.0_dl*fldp(i,IRANGE)**2 + yscl**4*ggrad1)
		enddo			
		! Calculate dzeta
		dzeta_lat(2,IRANGE) = yscl**2 * zf1 / zf2
		dzeta(2) = sum(dzeta_lat(2, IRANGE))/dble(n1)/dble(n2)/dble(n3)

		!do i=1,n_sample
		!	sample_dzeta(2,i) = dzeta_lat(2,sample_site(i,1),sample_site(i,2),sample_site(i,3))
		!enddo
	end subroutine get_dzeta_part

! Subroutine to calculate dzeta using a smoothed stress energy tensor
! to do: use lat_smooth to smooth stress energy tensor
! to do: update variables for smoothing scale
	subroutine get_dzeta_smooth(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j
	
		real(dl), parameter :: eps_smooth = 1._dl

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		dzeta(1) = dzeta(2)	!store old value
		dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
		sample_dzeta(1,:) = sample_dzeta(2,:)

		zf1(IRANGE) = 0._dl
		zf2(IRANGE) = 0._dl
		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		

			zf1 = zf1 + ggrad1 + fldp(i,IRANGE)*zlap1
		enddo
		! Get denominator of dzeta
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			
			zf2 = zf2 + 3.0_dl*fldp(i,IRANGE)**2 + yscl**4*ggrad1
		enddo
		call lat_smooth(zf1, Fk1, eps_smooth, planf, planb)
		call lat_smooth(zf2, Fk1, eps_smooth, planf, planb)	
		! Calculate dzeta
		dzeta_lat(2,IRANGE) = yscl**2 * zf1 / zf2
		dzeta(2) = sum(dzeta_lat(2, IRANGE))/dble(n1)/dble(n2)/dble(n3)

		!do i=1,n_sample
		!	sample_dzeta(2,i) = dzeta_lat(2,sample_site(i,1),sample_site(i,2),sample_site(i,3))
		!enddo

	end subroutine get_dzeta_smooth

! Subroutine that takes in a field as calculated on the lattice and smoothing on some physical scale
! to do: Check Fourier convention on kernal
! to do: Check sig_s for 3 dimensions
	subroutine lat_smooth(f1, Fk1, eps, planf, planb)
		real(C_DOUBLE), pointer :: f1(:,:,:)								! Field to be smoothed
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)		! Pointer for FFT of field to be smoothed
		real(dl) :: eps																			! number times comoving hubble scale on which to smooth
		real(dl) :: sig_s																		! Smoothing scale (comoving size)
		real(dl) :: ker																			! Gaussian kernal for smoothing
		integer :: i,j,k,ii,jj,kk														! lattice labels and wave numbers
		real(dl) :: rad2																		! radius squared in Fourier space
		type(C_PTR) :: planf, planb													! FFTW forward and backward plans

		! compute sig_s
		sig_s = -eps/(ysclp / 6. / yscl)

		! compute forward FFT
		call fftw_execute_dft_r2c(planf, f1, Fk1)

		! Loop over lattice, multiply FFT by kernal
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx
					! compute kernal
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					ker = exp(-0.5*rad2*sig_s**2)
					Fk1(LATIND) = Fk1(LATIND)*ker
				enddo
			enddo
		enddo
		! compute backward FFT and normalize
		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1 = f1/nvol

	end subroutine lat_smooth


! calculates grad(f1).grad(f2) and outputs in gg (f1 and f2 are also overwritten)
! check if I can do this by updating a partial sum
	subroutine graddotgrad(nsize, dk, f1, f2, gg, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk
		real(C_DOUBLE), pointer :: f1(:,:,:), f2(:,:,:), gg(:,:,:)!, g1(:,:,:), g2(:,:,:), g3(:,:,:)
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:), Fk3(:,:,:)!, Fk(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: nn1, nn2, nn3
    integer :: i,j,k,ii,jj,kk
		!real(C_DOUBLE) :: kcur!new
    complex(C_DOUBLE_COMPLEX), parameter :: i_imag = (0.,1.)!new
		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1

		call fftw_execute_dft_r2c(planf, f1, Fk1)
		call fftw_execute_dft_r2c(planf, f2, Fk2)

		!do k=1,n3; if (k>nn3) then; kk=nn3+1-k;

		!!!!!!!! putting in i_imag instead of factoring it out, to test if the c2r fft handles things properly if i_imag is factored out
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble((i-1)*dk)*Fk1(i,j,k)
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble((i-1)*dk)*Fk2(i,j,k)
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		!gg = -f1*f2
		gg = f1*f2 !using i_imag

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(jj*dk)*Fk1(i,j,k)
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(jj*dk)*Fk2(i,j,k)
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)
		!gg = gg-f1*f2
		gg = gg+f1*f2 !using i_imag
		
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(kk*dk)*Fk1(i,j,k)
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(kk*dk)*Fk2(i,j,k)
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)
		!gg = gg-f1*f2
		gg = gg+f1*f2 !using i_imag
		gg = gg /dble(n1)**2 /dble(n2)**2 /dble(n3)**2

!		!!!!Alrternative calculation with fewer lines, more variables!!!!
!		call fftw_execute_dft_r2c(planf, f1, Fk)
!		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
!			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
!				do i=1,nn1
					! this is now a vector so I will need separate DFTs
!					Fk1(i,j,k) = dble((i-1)*dk)*Fk(i,j,k)
!					Fk2(i,j,k) = dble(jj*dk)*Fk(i,j,k)
!					Fk3(i,j,k) = dble(kk*dk)*Fk(i,j,k)					
!				enddo
!			enddo
!		enddo
!		call fftw_execute_dft_c2r(planb, Fk1, g1)
!		call fftw_execute_dft_c2r(planb, Fk2, g2)
!		call fftw_execute_dft_c2r(planb, Fk3, g3)
!		
!		call fftw_execute_dft_r2c(planf, f2, Fk)
!		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
!			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
!				do i=1,nn1
!					! this is now a vector so I will need separate DFTs
!					Fk1(i,j,k) = dble((i-1)*dk)*Fk(i,j,k)
!					Fk2(i,j,k) = dble(jj*dk)*Fk(i,j,k)
!					Fk2(i,j,k) = dble(kk*dk)*Fk(i,j,k)					
!				enddo
!			enddo
!		enddo
!		call fftw_execute_dft_c2r(planb, Fk1, f1)
!		call fftw_execute_dft_c2r(planb, Fk2, f2)
!		call fftw_execute_dft_c2r(planb, Fk3, gg)

!		f1(:,:,:) = -f1(:,:,:)*g1(:,:,:) - f2(:,:,:)*g2(:,:,:) - gg(:,:,:)*g3(:,:,:)
!		f1 = f1 /dble(n1) /dble(n2) /dble(n3)

	end subroutine graddotgrad

! function to calculate averaged characteristic function
!	function char_func(f1)
!		real(dl) :: char_func
!		real(dl) :: f1(:,:,:)
!		real(dl) :: char_lat(:,:,:)
!		real(dl) :: f1_mean
!		integer, dimension(1:3) :: n_size = (/nx,ny,nz/)
		
!		f1_mean = get_mean_3d(f1, n_size)
!		char_lat(IRANGE) = dexp(-log(yscl)*(f1(IRANGE)-f1_mean))
!		char_func = sum(char_lat(IRANGE))/dble(nx)/dble(ny)/dble(nz)
!!!!! complete this function
!	end function char_func

!!! Add characteristic function to output
	subroutine write_zeta(time)
		real(dl) :: time
		real(dl) :: zeta_moment1, zeta_moment2, zeta_moment3, zeta_moment4!, zeta_char
		integer, dimension(1:3) :: n_size = (/nx,ny,nz/)

		!n_size(1) = nx; n_size(2) = ny; n_size(3) = nz
		
		zeta_moment1 = get_mean_3d(zeta_lat, n_size)
		zeta_moment2 = get_moment_3d(zeta_lat, zeta_moment1, 2, n_size)
		zeta_moment3 = get_moment_3d(zeta_lat, zeta_moment1, 3, n_size)
		zeta_moment4 = get_moment_3d(zeta_lat, zeta_moment1, 4, n_size)
		!zeta_char = char_func(zeta_lat)
		write(96,'(30(ES22.15,2x))') time, dzeta(1), dzeta(2), zeta, zeta_moment1, zeta_moment2, zeta_moment3, zeta_moment4!, zeta_char
		!write(96,*)
	end subroutine write_zeta

! Subroutine to output moments of a field on the lattice
! to do: 
	subroutine write_moments(f1, n_file)
		real(C_DOUBLE), pointer :: f1(:,:,:)								! Field to calculate moments
		integer :: n_file																		! File number of output
		real(dl) :: moment1, moment2, moment3, moment4
		integer, dimension(1:3) :: n_size = (/nx,ny,nz/)
		
		moment1 = get_mean_3d(f1, n_size)
		moment2 = get_moment_3d(f1, moment1, 2, n_size)
		moment3 = get_moment_3d(f1, moment1, 3, n_size)
		moment4 = get_moment_3d(f1, moment1, 4, n_size)
		write(n_file,'(30(ES22.15,2x))') moment1, moment2, moment3, moment4
	end subroutine write_moments

! Subroutine to output variables dealing with partial zetas on a slice of the lattice
! to do: call this subroutine from evolve_spectral.f90
	subroutine write_zeta_partial(time, kslice, n_file)
		real(dl), intent(in) :: time	! conformal time of output
		integer :: kslice							! the z index of the slice to be output
		integer :: n_file							! file number of output
		integer :: i,j								! lattice index (x-y slice)

		do i=1,nx; do j=1,ny
#ifdef ONEFLD
			write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, epsilon_part(1,i,j,kslice), dzeta_part(2,1,i,j,kslice),&
			dzeta_part(2,2,i,j,kslice), zeta_part(1,i,j,kslice), zeta_part(2,i,j,kslice)
#elif TWOFLD2
			write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, epsilon_part(1,i,j,kslice), dzeta_part(2,1,i,j,kslice),&
			dzeta_part(2,2,i,j,kslice), zeta_part(1,i,j,kslice), zeta_part(2,i,j,kslice), epsilon_part(2,i,j,kslice), dzeta_part(2,3,i,j,kslice),&
			dzeta_part(2,4,i,j,kslice), zeta_part(3,i,j,kslice), zeta_part(4,i,j,kslice)!, zeta_part(1,i,j,kslice)+zeta_part(2,i,j,kslice)
#endif
		enddo; enddo
		write(n_file,*)

	end subroutine write_zeta_partial

end module zeta_mod







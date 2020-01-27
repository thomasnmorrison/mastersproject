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

! to do: write option in smoothing functions for normalization by are, or peak value


module zeta_mod

!	use, intrinsic :: iso_c_binding
!  include 'fftw3.f03'
!	use fftw3
!	use params
#include "macros.h"
	use hamiltonian
	use sample_sites_mod
	use moments_mod

	integer, parameter :: nzpart = 4												! number of zeta parts
	integer, parameter :: n_scale = 3												! number of filtering scales
	integer, parameter :: n_moment = 4											! number of moments to be calculated

	real(dl), dimension(n_scale) :: k_hor_smooth = (/0.25_dl,1.0_dl,4.0_dl/)
	integer, dimension(n_scale) :: ord = (/3,3,3/)

	real(dl) :: zeta = 0.0_dl
	!real(dl), dimension(2) :: dzeta = (/0.0_dl,0.0_dl/)
	real(dl) :: dzeta = 0.0_dl															!	dzeta/dtau lattice average
	real(dl), dimension(IRANGE) :: zeta_lat									! zeta on the lattice
	real(dl), dimension(IRANGE) :: zeta_pert
	real(dl), dimension(IRANGE) :: zeta_pert_rho
	real(dl), dimension(n_scale,IRANGE) :: zeta_smooth							! zeta smoothed on the lattice
	!real(dl), dimension(2,IRANGE) :: dzeta_lat 						!	dzeta/dtau
	real(dl), dimension(IRANGE) :: dzeta_lat 								!	dzeta/dtau for first order integrator
	real(dl), dimension(IRANGE) :: dzeta_smooth							! dzeta/dtau smoothed on the lattice
	real(dl), dimension(IRANGE) :: zeta_smooth_pert					! dzeta/dtau smoothed on the lattice
	!real(dl), dimension(2,IRANGE) :: dzeta_lat_temp				! 1st component is numerator of dzeta, 2nd component is denominator of dzeta
	real(dl), dimension(2,IRANGE) :: epsilon_part						! components rho*x_A*epsilon_A, slpit as kinetic, gradient
	!real(dl), dimension(2*nfld,IRANGE) :: zeta_part				! zeta component calculated from individual fields (and linear/non-linear)
	!real(dl), dimension(2,2*nfld,IRANGE) :: dzeta_part			! dzeta/dtau component calculated from individual fields
	real(dl), dimension(nzpart,IRANGE) :: zeta_part					! zeta component calculated from K+V, G, K, G+V
	real(dl), dimension(nzpart,IRANGE) :: dzeta_part				! dzeta/dtau component calculated from K+V, G, K, G+V

	real(dl), dimension(IRANGE) :: delta_alpha							! perturbation in alpha
	real(dl), dimension(IRANGE) :: delta_alpha0_prime				! initial perturbation in alpha in the uniform rho gauge

	real(dl) :: k2_kick_low = 0.0_dl													! lower k bound on zeta kick
	real(dl) :: k2_kick_low_pert = 0.0_dl											! lower k bound on zeta kick
	
	! additional pointers needed for calculating dzeta
	real(C_DOUBLE), pointer :: ggrad1(:,:,:)
	real(C_DOUBLE), pointer :: zf1(:,:,:)
	real(C_DOUBLE), pointer :: zf2(:,:,:)
	real(C_DOUBLE), pointer :: zlap1(:,:,:)
	real(C_DOUBLE), pointer :: zlap2(:,:,:)
	complex(C_DOUBLE_COMPLEX), pointer :: Fk3(:,:,:)

contains

	! Subroutine to initialize zeta and its currents and derivatives to zero
	subroutine zeta_init_zero()
		zeta = 0.0_dl
		dzeta = 0.0_dl	!dzeta(:) = 0.0_dl
		zeta_lat(IRANGE) = 0.0_dl
		dzeta_lat(IRANGE) = 0.0_dl	!dzeta_lat(:,IRANGE) = 0.0_dl
		zeta_part(:,IRANGE) = 0.0_dl
		dzeta_part(:,IRANGE) = 0.0_dl
		zeta_smooth(:,IRANGE) = 0.0_dl
		dzeta_smooth(IRANGE) = 0.0_dl
	end subroutine zeta_init_zero

	! Subroutine to initialize zeta to the linear perturbation result ofver the lattice
	! Initializes only fluctuations
	! to do: test initialization of delta_alpha and delta_alpha0_prime
	subroutine zeta_init_pert(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:), Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		real(dl) :: rho_p_0

		zeta = 0.0_dl
		dzeta = 0.0_dl
		!zeta_lat(IRANGE) = get_hubble()*nvol*yscl**3/sum(fldp(1,IRANGE))*fld(1,IRANGE)
		zeta_lat(IRANGE) = 0.0_dl
		dzeta_lat(IRANGE) = 0.0_dl
		zeta_part(:,IRANGE) = 0.0_dl
		dzeta_part(:,IRANGE) = 0.0_dl
		zeta_pert_rho(IRANGE) = 0.0_dl
		zeta_smooth(:,IRANGE) = 0.0_dl
		zeta_smooth_pert(IRANGE) = 0.0_dl
		dzeta_smooth(IRANGE) = 0.0_dl
		delta_alpha(IRANGE) = 0.0_dl	! new

		! Get rho, rho+p
		rho_p_0 = 0.0_dl
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			rho_p_0 = rho_p_0 + sum(fldp(i,IRANGE)**2/yscl**6 + ggrad1/(3._dl*yscl**2))/nvol
			zeta_lat(IRANGE) = zeta_lat(IRANGE) + 0.5_dl*fldp(i,IRANGE)**2/yscl**6 + 0.5_dl*ggrad1/yscl**2
		enddo
#ifdef ONEFLD
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + potential_test(fld(1,IRANGE))
#elif TWOFLD2
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + potential_test(fld(1,IRANGE),fld(2,IRANGE))
#endif
		! Calculate initial zeta
		zeta_lat = (zeta_lat - sum(zeta_lat)/nvol)/(3._dl*rho_p_0)
		delta_alpha0_prime = zeta_lat!/get_hubble()	! new
	end subroutine zeta_init_pert

	! Subroutine to initialize zeta based on its gradients
	! to do: declare variables
	! to do: write a gradient subroutine
	! to do: calculate rho
	! to do: calculate rho+P
	! to do: loop over lattice, calculate FFT of zeta
	! to do: check sign conension being used for FFT
	! to do: normalization constant
!	subroutine zeta_init(f1, Fk1, planf, planb)
!		real(C_DOUBLE), pointer :: f1(:,:,:)
!		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:)
!		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
!		type(C_PTR) :: planf, planb																		! fftw plans
!	
!		real(dl), parameter :: norm = 						! normalization for zeta
!		complex(C_DOUBLE_COMPLEX), parameter :: i_imag = (0._dl, 1._dl)
!
!		! calculate rho, rho+P
!		do i = 1,nfld
!
!		enddo
!		
!		! kk .neq. 0
!		! get grad 
!		do k=2,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
!			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
!				do i=1,nnx; ii = i-1
!					Fk1(LATIND) = -i_imag*Fk2(LATIND)/dble(kk)		! check sign
!				enddo
!			enddo
!		enddo
!		! kk = 0; jj .neq. 0
!		k = 1; kk = 0
!		do j=2,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
!			do i=1,nnx; ii = i-1
!				Fk1(LATIND) = -i_imag*Fk2(LATIND)/dble(jj)		! check sign
!			enddo
!		enddo
!		! kk = 0; jj = 0; ii .neq. 0
!		k = 1; kk = 0; j = 1; jj = 0
!		do i=1,nnx; ii = i-1
!			Fk1(LATIND) = -i_imag*Fk2(LATIND)/dble(ii)		! check sign
!		enddo
!		! kk = 0 ; jj = 0; ii = 0
!		k = 1; kk = 0; j = 1; jj = 0; i = 1; ii = 0		
!		Fk1(LATIND) = (0._dl, 0._dl)								! set mean zeta=0
!		
!		call fftw_execute_dft_c2r(planb, Fk1, f1)
!		zeta_lat = f1/norm
!
!	end subroutine zeta_init

	! Subroutinte that gives a "stochastic kick" of zeta going from the fluctuations to the condensate
	! Initializes zeta from the linear perturbation equation zeta = H/dphi*delta phi within a 
	! specified k-band and adds it to the condensate (zeta_smooth)
	subroutine zeta_kick_pert_h(dk, k2_low, k_high_fac, Fk1, planf, planb)
		real(dl) :: dk
		real(dl) :: k2_low, k2_high, k_high_fac																	! lower k bound and upper k bound as horizon fraction
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)
		type(C_PTR) :: planf, planb																		! fftw plans
		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
		real(dl) :: rad2																							! k-space radius squared

		zlap1(IRANGE) = fld(1,IRANGE)													! this part is assuming one field for linear pert theory
		call fftw_execute_dft_r2c(planf, zlap1, Fk1)
		k2_high = (k_high_fac*(yscl*get_hubble()))**2
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					if (rad2 <= k2_low .or. rad2 > k2_high) then
						Fk1(LATIND) = (0._dl, 0._dl)
					endif
				enddo
			enddo
		enddo
		Fk1 = Fk1*get_hubble()*nvol*yscl**3/sum(fldp(1,IRANGE))
		k2_low = k2_high														! store upper k bound for next iteration's lower kbound
		call fftw_execute_dft_c2r(planb, Fk1, zlap1)
		zeta_smooth_pert = zeta_smooth_pert + zlap1/nvol
	end subroutine zeta_kick_pert_h

	! Subroutine to move zeta modes from within a specified k-band into the condensate
	! calculates zeta perturbatively at each step
	subroutine zeta_kick_pert(dk, k2_low, k_high_fac, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dk
		real(dl) :: k2_low, k2_high, k_high_fac												! lower k bound and upper k bound as horizon fraction
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:), Fk3(:,:,:)
		type(C_PTR) :: planf, planb																		! fftw plans

		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
		real(dl) :: rad2																							! k-space radius squared
		real(dl) :: rho_p_0

		! Get rho, rho+p
		rho_p_0 = 0.0_dl
		laplace(IRANGE) = 0.0_dl
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad([nx,ny,nz],dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			rho_p_0 = rho_p_0 + sum(fldp(i,IRANGE)**2/yscl**6 + ggrad1/(3._dl*yscl**2))/nvol							! 3(rho+P) term
			laplace(IRANGE) = laplace(IRANGE) + 0.5_dl*fldp(i,IRANGE)**2/yscl**6 + 0.5_dl*ggrad1/yscl**2	! rho term
			dzeta_smooth(IRANGE) = fldp(i,IRANGE)**2/yscl**6 + ggrad1/yscl**2/3._dl	! testing
		enddo
#ifdef ONEFLD
		laplace(IRANGE) = laplace(IRANGE) + potential_test(fld(1,IRANGE))																! rho term
#elif TWOFLD2
		laplace(IRANGE) = laplace(IRANGE) + potential_test(fld(1,IRANGE),fld(2,IRANGE))									! rho term
#endif
		laplace = (laplace - sum(laplace)/nvol)/ (3._dl*rho_p_0)
		zeta_pert = laplace
		call delta_alpha_step(dk, Fk1, Fk2, Fk3, planf, planb)

		zlap1(IRANGE) = zeta_pert(IRANGE) + delta_alpha(IRANGE)
		k2_high = (k_high_fac*(yscl*get_hubble()))**2

		call filter_sharp_k(dk, k2_low, k2_high, zlap1, Fk1, planf, planb)

		k2_low = k2_high														! store upper k bound for next iteration's lower kbound
		zeta_smooth_pert = zeta_smooth_pert + zlap1!/nvol
	end subroutine zeta_kick_pert

	! Subroutine to move zeta modes from within a specified k-band into the condensate
	! calculates zeta perturbatively at each step by using delta rho and rho dot
	subroutine zeta_kick_pert_rho(dk, k2_low, k_high_fac, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dk
		real(dl) :: k2_low, k2_high, k_high_fac												! lower k bound and upper k bound as horizon fraction
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:), Fk3(:,:,:)
		type(C_PTR) :: planf, planb																		! fftw plans

		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
		real(dl) :: rad2																							! k-space radius squared

		zeta_pert_rho(IRANGE) = 0.0_dl
		laplace(IRANGE) = 0.0_dl
		!print*,"seg check1"
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad([nx,ny,nz],dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			laplace(IRANGE) = laplace(IRANGE) + 3._dl*fldp(i,IRANGE)**2/yscl**6 + ggrad1/yscl**2											! 3(rho+P)
			zeta_pert_rho(IRANGE) = zeta_pert_rho(IRANGE) + 0.5_dl*fldp(i,IRANGE)**2/yscl**6 + 0.5_dl*ggrad1/yscl**2	! 
			dzeta_smooth(IRANGE) = fldp(i,IRANGE)**2/yscl**6 + ggrad1/yscl**2/3._dl	! testing
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad([nx,ny,nz],dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			!print*,"seg check2"
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(nx,ny,nz,zlap1,Fk1,dk,planf,planb)
			laplace(IRANGE) = laplace(IRANGE) - ggrad1(IRANGE)/yscl**5 - zlap1(IRANGE)/yscl**5
			!laplace(IRANGE) = zeta_pert_rho(IRANGE) - ggrad1(IRANGE)/yscl**5 - zlap1(IRANGE)/yscl**5
			!print*,"seg check3"
		enddo
#ifdef ONEFLD
		zeta_pert_rho(IRANGE) = zeta_pert_rho(IRANGE) + potential_test(fld(1,IRANGE))
#elif TWOFLD2
		zeta_pert_rho(IRANGE) = zeta_pert_rho(IRANGE) + potential_test(fld(1,IRANGE),fld(2,IRANGE))
#endif
		zeta_pert_rho = (zeta_pert_rho - sum(zeta_pert_rho)/nvol)/ (laplace)
		!zeta_pert_rho = laplace		
		!print*,"seg check4"
		zlap1(IRANGE) = zeta_pert_rho(IRANGE)
		k2_high = (k_high_fac*(yscl*get_hubble()))**2
		call filter_sharp_k(dk, k2_low, k2_high, zlap1, Fk1, planf, planb)
		k2_low = k2_high														! store upper k bound for next iteration's lower kbound
		zeta_smooth_pert = zeta_smooth_pert + zlap1
		!print*,"seg check5"
	end subroutine zeta_kick_pert_rho

	! Subroutine to move zeta modes from within a specified k-band into the condensate
	! Assumes zeta_lat is undergoing a full time evolution from proper initial conditions
	subroutine zeta_kick(dk, k2_low, k_high_fac, Fk1, planf, planb)
		real(dl) :: dk
		real(dl) :: k2_low, k2_high, k_high_fac												! lower k bound and upper k bound as horizon fraction
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)
		type(C_PTR) :: planf, planb																		! fftw plans
		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
		real(dl) :: rad2																							! k-space radius squared

		zlap1(IRANGE) = zeta_lat(IRANGE)! - sum(zeta_lat(IRANGE))/nvol
		!call fftw_execute_dft_r2c(planf, zlap1, Fk1)
		k2_high = (k_high_fac*(yscl*get_hubble()))**2

		call filter_sharp_k(dk, k2_low, k2_high, zlap1, Fk1, planf, planb)

		!do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
		!	do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
		!		do i=1,nnx; ii = i-1
		!			rad2 = (ii**2 + jj**2 + kk**2)*dk**2
		!			if (rad2 <= k2_low .or. rad2 > k2_high) then
		!				Fk1(LATIND) = (0._dl, 0._dl)
		!			endif
		!		enddo
		!	enddo
		!enddo

		k2_low = k2_high														! store upper k bound for next iteration's lower kbound
		!call fftw_execute_dft_c2r(planb, Fk1, zlap1)
		zeta_smooth_pert = zeta_smooth_pert + zlap1!/nvol
	end subroutine zeta_kick

	! Set smoothed zeta by filtering the fine-grained zeta with a sharp-k filter
	subroutine zeta_filter_sharp_k(dk, k_high_fac, Fk1, ind, planf, planb)
		real(dl) :: dk
		real(dl) :: k2_high, k_high_fac																! upper k bound as horizon fraction
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)
		integer, intent(in) :: ind																		! scale index of zeta_smooth
		type(C_PTR) :: planf, planb																		! fftw plans
		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
		real(dl) :: rad2																							! k-space radius squared

		zlap1(IRANGE) = zeta_lat(IRANGE)
		k2_high = (k_high_fac*(yscl*get_hubble()))**2

		call filter_sharp_k(dk, 0.0_dl, k2_high, zlap1, Fk1, planf, planb)

		zeta_smooth(ind, IRANGE) = zlap1(IRANGE)
	end subroutine zeta_filter_sharp_k

	! Set smoothed zeta by filtering the fine-grained zeta with a sharp-k filter
	subroutine zeta_filter_sgauss_k(dk, k_high_fac, Fk1, ord, ind, planf, planb)
		real(dl) :: dk
		real(dl) :: k2_high, k_high_fac																! k standard deviation as horizon fraction
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)
		integer, intent(in) :: ord, ind																! order of sgaussian filter, scale index of zeta_smooth
		type(C_PTR) :: planf, planb																		! fftw plans
		integer :: i,j,k,ii,jj,kk																			! lattice labels and wave numbers
		real(dl) :: rad2																							! k-space radius squared

		zlap1(IRANGE) = zeta_lat(IRANGE)
		k2_high = (k_high_fac*(yscl*get_hubble()))**2

		call filter_sgauss_k(dk, 0.0_dl, k2_high, zlap1, Fk1, ord, planf, planb)

		zeta_smooth(ind, IRANGE) = zlap1(IRANGE)
	end subroutine zeta_filter_sgauss_k

	! Subroutine to step zeta total without referencing zeta parts
	subroutine zeta_step(dt, nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt
		integer, dimension(1:3), intent(in) :: nsize
		integer :: i
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		call get_dzeta(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		zeta = zeta + dzeta*dt
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + dzeta_lat(IRANGE)*dt

	end subroutine zeta_step

	! Subroutine to step zeta intgration using Euler method with K, G, K+V, G+V split (also steps zeta total)
	! n.b. this subroutine uses a only one value of dzeta
	! n.b. this subroutine integrates zeta_part without weighting
	subroutine zeta_step_kgv(dt, nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt
		integer, dimension(1:3), intent(in) :: nsize
		integer :: i
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		call get_dzeta_part_kgv(nsize, dk, Fk1, Fk2, Fk3, planf, planb)		! call for dzeta calc with K, G, V components

		zeta = zeta + dzeta*dt
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + dzeta_lat(IRANGE)*dt
		zeta_part(1,IRANGE) = zeta_part(1,IRANGE) + dzeta_part(1,IRANGE)*dt	! K unweighted
		zeta_part(2,IRANGE) = zeta_part(2,IRANGE) + dzeta_part(2,IRANGE)*dt	! G unweighted
		zeta_part(3,IRANGE) = zeta_part(3,IRANGE) + dzeta_part(3,IRANGE)*dt	! K+V unweighted
		zeta_part(4,IRANGE) = zeta_part(4,IRANGE) + dzeta_part(4,IRANGE)*dt	! G+V unweighted
	end subroutine zeta_step_kgv

	! Subroutine to step zeta integration with K+V, G split both unsmoothed and smoothed
	! to do: allow for different smoothings, ie T_{i0}, fields, dzeta
	subroutine zeta_step_kgv_smooth(dt, nsize, dk, r_smooth, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt				! time step
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk					! fundemental wavenumber
		real(dl) :: r_smooth	! Smoothing scale in hubbles
		!integer :: ind				! index to select smoothing type
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		!if (ind==1) then
			call get_dzeta_part_kgv_smooth(nsize, dk, r_smooth, Fk1, Fk2, Fk3, planf, planb)	
		!endif

		zeta = zeta + dzeta*dt																							! step volume average zeta
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + dzeta_lat(IRANGE)*dt					! step zeta on the lattice
		zeta_smooth(1,IRANGE) = zeta_smooth(1,IRANGE) + dzeta_smooth(IRANGE)*dt	! step smoothed zeta
		zeta_part(:,IRANGE) = zeta_part(:,IRANGE) + dzeta_part(:,IRANGE)*dt	! step partial zeta
	end subroutine zeta_step_kgv_smooth

	! Subroutine to step zeta with fileds split and interaction potential equally shared (also steps zeta total)
	subroutine zeta_step_fld_vint(dt, nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt
		integer, dimension(1:3), intent(in) :: nsize
		integer :: i
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		call get_dzeta_part_fld_vint(nsize, dk, Fk1, Fk2, Fk3, planf, planb)		! call for dzeta calc with K, G, V components

		zeta = zeta + dzeta*dt
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + dzeta_lat(IRANGE)*dt
		do i=1,nfld
			zeta_part(i,IRANGE) = zeta_part(i,IRANGE) + dzeta_part(i,IRANGE)*dt
		enddo
	end subroutine zeta_step_fld_vint

	! Subroutine to calculate delta_alpha given delta_alpha0_prime is stored
	subroutine delta_alpha_step(dk, Fk1, Fk2, Fk3, planf, planb)
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: i,j
		real(dl) :: rho_p_0

		! Get rho, rho+p
		rho_p_0 = 0.0_dl
		laplace(IRANGE) = 0.0_dl
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad([nx,ny,nz],dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			laplace(IRANGE) = laplace(IRANGE) + 0.5_dl*fldp(i,IRANGE)**2/yscl**6 + 0.5_dl*ggrad1/yscl**2
			rho_p_0 = rho_p_0 + sum(fldp(i,IRANGE)**2/yscl**6 + ggrad1/yscl**2/3._dl)/nvol
		enddo
#ifdef ONEFLD
		laplace(IRANGE) = laplace(IRANGE) + potential_test(fld(1,IRANGE))
#elif TWOFLD2
		laplace(IRANGE) = laplace(IRANGE) + potential_test(fld(1,IRANGE),fld(2,IRANGE))
#endif
		delta_alpha = delta_alpha0_prime -(laplace - sum(laplace)/nvol)/(3._dl*rho_p_0)
	end subroutine delta_alpha_step

! subroutine to get dzeta for an arbitrary number of fields
! n.b. I haven't tested this subroutine
	subroutine get_dzeta(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		dzeta_lat(IRANGE) = 0._dl
		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		

			dzeta_lat = dzeta_lat + ggrad1/yscl**4 + fldp(i,IRANGE)*zlap1/yscl**4
		enddo
		! Get denominator of dzeta
		laplace = 0.0_dl
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			laplace(IRANGE) = laplace(IRANGE) + fldp(i,IRANGE)**2/yscl**6 + ggrad1/(3._dl*yscl**2)
			!epsilon_part(i,IRANGE) = fldp(i,IRANGE)**2/yscl**6 + ggrad1/(3._dl*yscl**2)
		enddo			
		! Calculate dzeta
		!dzeta_lat(IRANGE) = dzeta_lat(IRANGE)/(3._dl*sum(epsilon_part(:,IRANGE), dim=1))
		!dzeta_lat(IRANGE) = dzeta_lat(IRANGE) - 3._dl*get_hubble()*laplace(IRANGE)*yscl!/(3._dl*laplace(IRANGE))	! denominator removed for testing
		dzeta_lat(IRANGE) = dzeta_lat(IRANGE)/(3._dl*laplace(IRANGE))
		dzeta_lat(IRANGE) = dzeta_lat(IRANGE) - sum(dzeta_lat(IRANGE))/nvol
		dzeta = sum(dzeta_lat(IRANGE))/dble(n1)/dble(n2)/dble(n3)

	end subroutine get_dzeta

! Subroutine to calculate dzeta and the components of dzeta for kinetic, gradient and potential terms
! n.b. watch the factors of yscl and to what factors they are applied
	subroutine get_dzeta_part_kgv(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		epsilon_part = 0.0_dl	!clear previous calculation
		dzeta_part = 0.0_dl		!clear previous calculation
		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		
			
			dzeta_part(1,IRANGE) = dzeta_part(1,IRANGE) + fldp(i,IRANGE)*(zlap1/yscl**4 - modeldv_test(i,fld(1,IRANGE),fld(2,IRANGE))/yscl**2)! K
			dzeta_part(2,IRANGE) = dzeta_part(2,IRANGE) + ggrad1/yscl**4																																			! G
			dzeta_part(3,IRANGE) = dzeta_part(3,IRANGE) + fldp(i,IRANGE)*zlap1/yscl**4																												! K+V
			dzeta_part(4,IRANGE) = dzeta_part(2,IRANGE) + ggrad1/yscl**4 + fldp(i,IRANGE)*modeldv_test(i,fld(1,IRANGE),fld(2,IRANGE))/yscl**2	! G+V
		enddo
		! Get denominator of dzeta
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			epsilon_part(1,IRANGE) = epsilon_part(1,IRANGE) + fldp(i,IRANGE)**2/yscl**6				! K
			epsilon_part(2,IRANGE) = epsilon_part(2,IRANGE) + ggrad1/(3._dl*yscl**2)					! G
		enddo			
		! Calculate dzeta
		dzeta_lat(IRANGE) = (dzeta_part(2,IRANGE)+dzeta_part(3,IRANGE)) / (3._dl*sum(epsilon_part(:,IRANGE) ,dim=1))
		dzeta = sum(dzeta_lat(IRANGE))/dble(n1)/dble(n2)/dble(n3)
		dzeta_part(1,IRANGE) = dzeta_part(1,IRANGE)/(3._dl*epsilon_part(1,IRANGE))	! K
		dzeta_part(2,IRANGE) = dzeta_part(2,IRANGE)/(3._dl*epsilon_part(2,IRANGE))	! G
		dzeta_part(3,IRANGE) = dzeta_part(3,IRANGE)/(3._dl*epsilon_part(1,IRANGE))	! K+V
		dzeta_part(4,IRANGE) = dzeta_part(4,IRANGE)/(3._dl*epsilon_part(2,IRANGE))	! G+V
		
	end subroutine get_dzeta_part_kgv

	! Subroutine to calculate dzeta and the components of dzeta for kinetic+potential, and gradient terms
	! for lattice average, per lattice site, and smoothed drho and rho+P on hubble scales
	! to do: add an input to select type of smoothing
	subroutine get_dzeta_part_kgv_smooth(nsize, dk, r_smooth, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk
		real(dl) :: r_smooth
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		epsilon_part = 0.0_dl	!clear previous calculation
		dzeta_part = 0.0_dl		!clear previous calculation
		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)			
			dzeta_part(1,IRANGE) = dzeta_part(1,IRANGE) + fldp(i,IRANGE)*zlap1/yscl**4	! K+V
			dzeta_part(2,IRANGE) = dzeta_part(2,IRANGE) + ggrad1/yscl**4								! G
		enddo
		! Get denominator of dzeta
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			epsilon_part(1,IRANGE) = epsilon_part(1,IRANGE) + fldp(i,IRANGE)**2/yscl**6	! K+V
			epsilon_part(2,IRANGE) = epsilon_part(2,IRANGE) + ggrad1/(3._dl*yscl**2)		! G
		enddo			
		! Calculate unsmoothed dzeta
		dzeta_lat(IRANGE) = (dzeta_part(1,IRANGE)+dzeta_part(2,IRANGE)) / (3._dl*sum(epsilon_part(:,IRANGE) ,dim=1))	! total
		dzeta = sum(dzeta_lat(IRANGE))/dble(n1)/dble(n2)/dble(n3)

		! Calculate smoothed dzeta
		zlap1 = dzeta_part(1,IRANGE)
		zlap2 = dzeta_part(2,IRANGE)
#ifdef GAUSS_SMOOTH
		call gauss_smooth_hub(zlap1, Fk1, Fk2, r_smooth, planf, planb)					! smoothing numerator, Gaussian
		call gauss_smooth_hub(zlap2, Fk1, Fk2, r_smooth, planf, planb)					! smoothing numerator, Gaussian
#endif
#ifdef K_TOPHAT_SMOOTH
		call sharp_k_smooth_hub(zlap1, Fk1, Fk2, r_smooth, planf, planb)					! smoothing numerator, sharp k
		call sharp_k_smooth_hub(zlap2, Fk1, Fk2, r_smooth, planf, planb)					! smoothing numerator, sharp k
#endif
		dzeta_part(3,IRANGE) = zlap1(IRANGE)															! smoothed K+V numerator
		dzeta_part(4,IRANGE) = zlap2(IRANGE)															! smoothed G numerator
		dzeta_smooth(IRANGE) = zlap1(IRANGE) + zlap2(IRANGE)							! smoothed total numerator
		zlap1 = epsilon_part(1,IRANGE)
		zlap2 = epsilon_part(2,IRANGE)
#ifdef GAUSS_SMOOTH
		call gauss_smooth_hub(zlap1, Fk1, Fk2, r_smooth, planf, planb)					! smoothing denominator, Gaussian
		call gauss_smooth_hub(zlap2, Fk1, Fk2, r_smooth, planf, planb)					! smoothing denominator, Gaussian
#endif
#ifdef K_TOPHAT_SMOOTH
		call sharp_k_smooth_hub(zlap1, Fk1, Fk2, r_smooth, planf, planb)					! smoothing denominator, sharp k
		call sharp_k_smooth_hub(zlap2, Fk1, Fk2, r_smooth, planf, planb)					! smoothing denominator, sharp k
#endif
		dzeta_part(3,IRANGE) = dzeta_part(3,IRANGE)/(3._dl*zlap1(IRANGE))	! smoothed K+V testing
		dzeta_part(4,IRANGE) = dzeta_part(4,IRANGE)/(3._dl*zlap2(IRANGE))	! smoothed G testing
		dzeta_smooth(IRANGE) = dzeta_smooth(IRANGE)/(3._dl*(zlap1(IRANGE)+zlap2(IRANGE)))	! smoothed total

		! Calculate dzeta parts
		dzeta_part(1,IRANGE) = dzeta_part(1,IRANGE)/(3._dl*epsilon_part(1,IRANGE))		! K+V
		dzeta_part(2,IRANGE) = dzeta_part(2,IRANGE)/(3._dl*epsilon_part(2,IRANGE))		! G

	end subroutine get_dzeta_part_kgv_smooth

! Subroutine to get dzeta from rho split into phi and chi parts with the interaction potential split half and half
! Also computes rho_A + P_A and dzeta_tot
	subroutine get_dzeta_part_fld_vint(nsize, dk, Fk1, Fk2, Fk3,planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

		! Get numerator of dzeta
		do i=1,nfld
			zlap1 = fldp(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
			zlap1 = fld(i,IRANGE)
			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		
			
			dzeta_part(i,IRANGE) = fldp(i,IRANGE)*zlap1/yscl**4 + ggrad1/yscl**4 
			dzeta_part(i,IRANGE) = dzeta_part(i,IRANGE) &
+ 0.5_dl*(fldp(1,IRANGE)*(-1._dl)**i*dV_int(fld(1,IRANGE),fld(2,IRANGE),1) + fldp(2,IRANGE)*(-1._dl)**(i+1)*dV_int(fld(1,IRANGE),fld(2,IRANGE),2))/yscl**2
		enddo
		
		! Get epsilon_part
		do i=1,nfld
			zlap1 = fld(i,IRANGE)
			zlap2 = fld(i,IRANGE)
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)

			epsilon_part(i,IRANGE) = fldp(i,IRANGE)**2/yscl**6 + ggrad1/(3._dl*yscl**2)
		enddo			

		! Calculate dzeta
		dzeta_lat(IRANGE) = sum(dzeta_part(:,IRANGE), dim=1)/(3._dl*sum(epsilon_part(:,IRANGE), dim=1))
		dzeta = sum(dzeta_lat(IRANGE))/dble(n1)/dble(n2)/dble(n3)
		do i=1, nfld
			dzeta_part(i,IRANGE) = dzeta_part(i,IRANGE)/(3._dl*epsilon_part(i,IRANGE))
		enddo

	end subroutine get_dzeta_part_fld_vint

! Subroutine to calculate dzeta using a smoothed stress energy tensor
! to do: use lat_smooth to smooth stress energy tensor
! to do: update variables for smoothing scale
!	subroutine get_dzeta_smooth(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
!		integer, dimension(1:3), intent(in) :: nsize
!		real*8 :: dk	
!    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
!		type(C_PTR) :: planf, planb

!		integer :: n1, n2, n3
!		integer :: i,j
	
!		real(dl), parameter :: eps_smooth = 1._dl

!		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)

!		!dzeta(1) = dzeta(2)	!store old value
!		!dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
!		!sample_dzeta(1,:) = sample_dzeta(2,:)

!		zf1(IRANGE) = 0._dl
!		zf2(IRANGE) = 0._dl
!		! Get numerator of dzeta
!		do i=1,nfld
!			zlap1 = fldp(i,IRANGE)
!			zlap2 = fld(i,IRANGE)
!			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
!			zlap1 = fld(i,IRANGE)
!			call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)		
!
!			zf1 = zf1 + ggrad1 + fldp(i,IRANGE)*zlap1
!		enddo
!		! Get denominator of dzeta
!		do i=1,nfld
!			zlap1 = fld(i,IRANGE)
!			zlap2 = fld(i,IRANGE)
!			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
!			
!			zf2 = zf2 + 3.0_dl*fldp(i,IRANGE)**2 + yscl**4*ggrad1
!		enddo
!		call gauss_smooth_hub(zf1, Fk1, eps_smooth, planf, planb)
!		call gauss_smooth_hub(zf2, Fk1, eps_smooth, planf, planb)	
!		! Calculate dzeta
!		!dzeta_lat(2,IRANGE) = yscl**2 * zf1 / zf2
!		!dzeta(2) = sum(dzeta_lat(2, IRANGE))/dble(n1)/dble(n2)/dble(n3)
!		dzeta_lat(IRANGE) = yscl**2 * zf1 / zf2
!		dzeta = sum(dzeta_lat(IRANGE))/dble(n1)/dble(n2)/dble(n3)

!	end subroutine get_dzeta_smooth

	! Subroutine that takes in a field as calculated on the lattice and smoothing on some physical scale
	subroutine gauss_smooth_hub(f1, Fk1, Fk2, r_smooth, planf, planb)
		real(C_DOUBLE), pointer :: f1(:,:,:)								! Field to be smoothed
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:)		! Pointer for FFT of field to be smoothed
		real(dl) :: r_smooth																! number times comoving hubble scale on which to smooth
		real(dl) :: sig_s																		! Smoothing scale (comoving size)
		real(dl) :: ker																			! Gaussian kernal for smoothing
		integer :: i,j,k,ii,jj,kk														! lattice labels and wave numbers
		real(dl) :: rad2																		! radius squared in Fourier space
		real(dl) :: ker_norm																! normalization factor
		type(C_PTR) :: planf, planb													! FFTW forward and backward plans

		sig_s = r_smooth/(yscl*get_hubble())								! compute sig_s = r_smooth/(aH)
		call fftw_execute_dft_r2c(planf, f1, Fk1)						! compute forward FFT
		
		! Loop over lattice, multiply FFT by kernal
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					ker = exp(-0.5*rad2*sig_s**2)									! compute kernal, nvol factor cancelled with inverse transform
					!Fk1(LATIND) = Fk1(LATIND)*ker									! convolve and make complex
					Fk2(LATIND) = ker*(1._dl,0._dl)
				enddo
			enddo
		enddo

		Fk1 = Fk1*Fk2
		call fftw_execute_dft_c2r(planb, Fk2, f1)
		ker_norm = sum(f1)
		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1 = f1/(nvol*ker_norm)

	end subroutine gauss_smooth_hub

	! Subroutine to take the devivative of scale when performing gaussian smoothing
	subroutine dgauss_smooth_hub(f1, Fk1, Fk2, r_smooth, planf, planb)
		real(C_DOUBLE), pointer :: f1(:,:,:)								! Field to be smoothed
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:)		! Pointer for FFT of field to be smoothed
		real(dl) :: r_smooth																! number times comoving hubble scale on which to smooth
		real(dl) :: sig_s																		! Smoothing scale (comoving size)
		real(dl) :: ker																			! Gaussian kernal for smoothing
		integer :: i,j,k,ii,jj,kk														! lattice labels and wave numbers
		real(dl) :: rad2																		! radius squared in Fourier space
		real(dl) :: ker_norm																! normalization factor
		type(C_PTR) :: planf, planb													! FFTW forward and backward plans

		sig_s = r_smooth/(yscl*get_hubble())								! compute sig_s = r_smooth/(aH)
		call fftw_execute_dft_r2c(planf, f1, Fk1)						! compute forward FFT
		
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					ker = exp(-0.5*rad2*sig_s**2)									! compute kernal for un-normalized window
					Fk2(LATIND) = ker*(1._dl,0._dl)
					ker = (-rad2*r_smooth/(yscl*get_hubble())**2)*exp(-0.5*rad2*sig_s**2)
					Fk1(LATIND) = Fk1(LATIND)*ker									! convolve and make complex
				enddo
			enddo
		enddo
		
		call fftw_execute_dft_c2r(planb, Fk2, f1)
		ker_norm = sum(f1)
		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1 = f1/(nvol*ker_norm)

	end subroutine dgauss_smooth_hub

	! Subroutine that takes in a field as calculated on the lattice and smoothing on some physical scale
	subroutine sharp_k_smooth_hub(f1, Fk1, Fk2, k_smooth, planf, planb)
		real(C_DOUBLE), pointer :: f1(:,:,:)								! Field to be smoothed, pointer for windowing function normalization
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:)		! Pointer for FFT of field to be smoothed
		real(dl) :: k_smooth																! number times horizon to smooth
		real(dl) :: k_ker																		! sharp k smoothing scale
		integer :: i,j,k,ii,jj,kk														! lattice labels and wave numbers
		real(dl) :: rad2																		! radius squared in Fourier space
		real(dl) :: ker_norm																! normalization factor
		type(C_PTR) :: planf, planb													! FFTW forward and backward plans

		k_ker = (k_smooth*(yscl*get_hubble()))**2
		call fftw_execute_dft_r2c(planf, f1, Fk1)						! compute forward FFT
		
		Fk2 = (1._dl, 0._dl)
		! Loop over lattice, cut modes outside of sharp k-space filter
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = dble(ii**2 + jj**2 + kk**2)*dk**2
					if (rad2>k_ker) then
						Fk2(LATIND) = (0._dl, 0._dl)
					endif
				enddo
			enddo
		enddo
		Fk1 = Fk1*Fk2
		call fftw_execute_dft_c2r(planb, Fk2, f1)
		ker_norm = sum(f1)
		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1 = f1/(nvol*ker_norm)
	end subroutine sharp_k_smooth_hub

	! Subroutine that applies a sharp-k filter to an input field
	subroutine filter_sharp_k(dk, k2_l, k2_h, f1, Fk1, planf, planb)
		real(dl) :: dk																		! fundemental k mode
		real(dl) :: k2_l, k2_h														! lower and upper bounds of filter squared
		real(C_DOUBLE), pointer :: f1(:,:,:)							! field to be filtered
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)	
		type(C_PTR) :: planf, planb												! fftw plans

		integer :: i,j,k,ii,jj,kk													! lattice labels and wave numbers
		real(dl) :: rad2																	! k-space radius squared

		call fftw_execute_dft_r2c(planf, f1, Fk1)
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					if (rad2 < k2_l .or. rad2 > k2_h) then
						Fk1(LATIND) = (0._dl, 0._dl)
					endif
				enddo
			enddo
		enddo

		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1=f1/nvol
	end subroutine filter_sharp_k

	! Subroutine that applies a super-gaussian filter in k space to an input field
	! n.b. Uses functional exp(0.5*(x/x_0)^(2*n))
	subroutine filter_sgauss_k(dk, k_mean, k2_var, f1, Fk1, ord, planf, planb)
		real(dl) :: dk																		! fundemental k mode
		real(dl) :: k_mean, k2_var												! k mean and variance of gaussian filter
		real(C_DOUBLE), pointer :: f1(:,:,:)							! field to be filtered
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)	
		integer:: ord																			! order of super-gaussian (1=gaussian)
		type(C_PTR) :: planf, planb												! fftw plans

		integer :: i,j,k,ii,jj,kk													! lattice labels and wave numbers
		real(dl) :: rad2																	! k-space radius squared

		call fftw_execute_dft_r2c(planf, f1, Fk1)
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					!Fk1(LATIND) = exp(-0.5_dl*(sqrt(rad2) - k_mean)**2/k2_var)*Fk1(LATIND)
					Fk1(LATIND) = exp(-0.5_dl*((sqrt(rad2) - k_mean)**2/k2_var)**ord)*Fk1(LATIND)	! testing super gaussian filter
				enddo
			enddo
		enddo

		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1=f1/nvol
	end subroutine filter_sgauss_k

	! Subroutine that applies a super-gaussian filter in k space to an input fourier transformed filed
	! This subroutine takes as input the fft of a field and returns the filtered field in realspace
	! without overwriting the the input fft
	! n.b. I haven't tested this subroutine
	subroutine filter_sgauss_k_fft(dk, k_mean, k2_var, f1, Fk1, Fk2, ord, planf, planb)
		real(dl) :: dk																		! fundemental k mode
		real(dl) :: k_mean, k2_var												! k mean and variance of gaussian filter
		real(C_DOUBLE), pointer :: f1(:,:,:)							! output field
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:)	! input field and filtered field	
		integer:: ord																			! order of super-gaussian (1=gaussian)
		type(C_PTR) :: planf, planb												! fftw plans

		integer :: i,j,k,ii,jj,kk													! lattice labels and wave numbers
		real(dl) :: rad2																	! k-space radius squared

		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					Fk2(LATIND) = exp(-0.5_dl*((sqrt(rad2) - k_mean)**2/k2_var)**ord)*Fk1(LATIND)
				enddo
			enddo
		enddo

		call fftw_execute_dft_c2r(planb, Fk2, f1)
		f1=f1/nvol
	end subroutine filter_sgauss_k_fft

	! Subroutine with the aim of outputting the real-space and k-space filters
	! to do: write output for the real-space and k-space filters
	subroutine filter_sgauss_r(dk, k_mean, k2_var, f1, Fk1, ord, planf, planb)
		real(dl) :: dk																		! fundemental k mode
		real(dl) :: k_mean, k2_var												! k mean and variance of gaussian filter
		real(C_DOUBLE), pointer :: f1(:,:,:)							! field to be filtered
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)	
		integer:: ord																			! order of super-gaussian (1=gaussian)
		type(C_PTR) :: planf, planb												! fftw plans

		integer :: i,j,k,ii,jj,kk													! lattice labels and wave numbers
		real(dl) :: rad2																	! k-space radius squared

		! Initialize k-space filter
		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					rad2 = (ii**2 + jj**2 + kk**2)*dk**2
					Fk1(LATIND) = exp(-0.5_dl*((sqrt(rad2) - k_mean)**2/k2_var)**ord)*(1._dl, 0._dl)
				enddo
			enddo
		enddo
		! Calculate real-space filter
		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1=f1/nvol 
	end subroutine filter_sgauss_r

	! Subroutine to calculate the ind^th component of the gradient of f1
	subroutine gradient(f1, Fk1, planf, planb, ind)
		real(C_DOUBLE), pointer :: f1(:,:,:)							! pointer for field
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)	! pointer for FFT
		type(C_PTR) :: planf, planb												! fftw plan
		integer :: ind																		! component of gradient

    integer :: i,j,k,ii,jj,kk
		complex(C_DOUBLE_COMPLEX), parameter :: i_imag = (0._dl, 1._dl)

		call fftw_execute_dft_r2c(planf, f1, Fk1)

		do k=1,nz; if(k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
			do j=1,ny; if(j>nny) then; jj = j-ny-1; else; jj=j-1; endif
				do i=1,nnx; ii = i-1
					if (ind == 1) then
						Fk1(LATIND) = i_imag*dble(ii)*dk*Fk1(LATIND)! check sign
					elseif (ind == 2) then
						Fk1(LATIND) = i_imag*dble(jj)*dk*Fk1(LATIND)! check sign
					elseif (ind == 3) then
						Fk1(LATIND) = i_imag*dble(kk)*dk*Fk1(LATIND)! check sign
					endif
				enddo
			enddo
		enddo
		
		call fftw_execute_dft_c2r(planb, Fk1, f1)
		f1 = f1/nvol	

	end subroutine gradient

	! calculates grad(f1).grad(f2) and outputs in gg (f1 and f2 are also overwritten)
	! check if I can do this by updating a partial sum
	subroutine graddotgrad_test(nsize, dk, f1, f2, gg, Fk1, Fk2, Fk3, planf, planb)
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

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = f1*f2 !using i_imag

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = gg+f1*f2 !using i_imag
		
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = gg+f1*f2 !using i_imag
		gg = gg /dble(n1)**2 /dble(n2)**2 /dble(n3)**2

	end subroutine graddotgrad_test

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

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = f1*f2 !using i_imag

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = gg+f1*f2 !using i_imag
		
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = gg+f1*f2 !using i_imag
		gg = gg /dble(n1)**2 /dble(n2)**2 /dble(n3)**2

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
	! to do: use arrays and loops like a decent human being
	! n.b. I've updated the formating of the output
	! Formatting is in equal time data blocks, each line in the data block being a different scale
	subroutine write_zeta(time, n_file)
		real(dl) :: time
		!real(dl) :: zeta_moment1, zeta_moment2, zeta_moment3, zeta_moment4!, zeta_char
		!real(dl) :: zeta_s_moment1, zeta_s_moment2, zeta_s_moment3, zeta_s_moment4	! moments of smoothed zeta
		integer, dimension(1:3) :: n_size = (/nx,ny,nz/)
		integer :: n_file		! file number
		!real(dl) :: mean_phi, delta_phi2
		!real(dl) :: zeta_s_mean
		real(dl), dimension(0:n_scale, n_moment) :: zeta_moment
		integer :: i, j

		!n_size(1) = nx; n_size(2) = ny; n_size(3) = nz
		
		!zeta_moment1 = get_mean_3d(zeta_lat, n_size)
		!zeta_moment2 = get_moment_3d(zeta_lat, zeta_moment1, 2, n_size)
		!zeta_moment3 = get_moment_3d(zeta_lat, zeta_moment1, 3, n_size)
		!zeta_moment4 = get_moment_3d(zeta_lat, zeta_moment1, 4, n_size)
		!zeta_s_moment1 = get_mean_3d(zeta_smooth, n_size)
		!zeta_s_moment2 = get_moment_3d(zeta_smooth, zeta_s_moment1, 2, n_size)
		!zeta_s_moment3 = get_moment_3d(zeta_smooth, zeta_s_moment1, 3, n_size)
		!zeta_s_moment3 = get_mean_3d(zeta_smooth_pert, n_size)
		!zeta_s_moment4 = get_moment_3d(zeta_smooth_pert, zeta_s_moment3, 2, n_size)
		!mean_phi = get_mean_3d(fld(1,IRANGE), n_size)
		!delta_phi2 = get_moment_3d(fld(1,IRANGE), mean_phi, 2, n_size)
		!zeta_char = char_func(zeta_lat)
		zeta_moment(0,1) = get_mean_3d(zeta_lat, n_size)
		do j=1, n_moment
			zeta_moment(0,j) = get_moment_3d(zeta_lat, zeta_moment(0,1), j, n_size)
		enddo
		write(n_file,'(30(ES22.15,2x))') time, yscl, zeta_moment(0,:)
		do i=1, n_scale
			zeta_moment(i,1) = get_mean_3d(zeta_smooth(i,:,:,:), n_size)
			do j=1, n_moment
				zeta_moment(i,j) = get_moment_3d(zeta_smooth(i,:,:,:), zeta_moment(i,1), j, n_size)
			enddo
			write(n_file,'(30(ES22.15,2x))') time, yscl, zeta_moment(i,:)
		enddo
		write(n_file, *)
		!write(n_file,'(30(ES22.15,2x))') time, dzeta, zeta, zeta_moment1, zeta_moment2, zeta_moment3, zeta_moment4
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
	subroutine write_zeta_partial(time, kslice, n_file)
		real(dl), intent(in) :: time	! conformal time of output
		integer :: kslice							! the z index of the slice to be output
		integer :: n_file							! file number of output
		integer :: i,j								! lattice index (x-y slice)

		do i=1,nx; do j=1,ny
#ifdef ONEFLD
			write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, epsilon_part(1,i,j,kslice), dzeta_part(1,i,j,kslice),&
			dzeta_part(2,i,j,kslice), zeta_part(1,i,j,kslice), zeta_part(2,i,j,kslice), epsilon_part(2,i,j,kslice), dzeta_part(3,i,j,kslice),&
			dzeta_part(4,i,j,kslice), zeta_part(3,i,j,kslice), zeta_part(4,i,j,kslice)
#elif TWOFLD2
			write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, epsilon_part(1,i,j,kslice), dzeta_part(1,i,j,kslice),&
			dzeta_part(2,i,j,kslice), zeta_part(1,i,j,kslice), zeta_part(2,i,j,kslice), epsilon_part(2,i,j,kslice), dzeta_part(3,i,j,kslice),&
			dzeta_part(4,i,j,kslice), zeta_part(3,i,j,kslice), zeta_part(4,i,j,kslice)!, zeta_part(1,i,j,kslice)+zeta_part(2,i,j,kslice)
#endif
		enddo; enddo
		write(n_file,*)

	end subroutine write_zeta_partial

	! Subroutine to write a slice of an array of filters
	! n.b. I'm abandoning this subroutine as it is not needed here and can easily be reproduced in post processing
	! to do: initialize file elsewhere
	! to do: call function to calculate filters
	! to do: output time, or a, or aH or something
	! to do: format into data blocks, each block at a particular time step, rows=position on lattice, column=which filter and r/k-space
	! to do: declace variable that will hold filter slices
	! to do: fill in k2_var
	! to do: fill in ord
!	subroutine write_filter(time,k_slice,n_file)
!		real(dl), intent(in) :: time	! conformal time of output
!		integer :: k_slice
!		integer :: n_file
!		integer :: i,j
!		real(dl) :: filt(nx,ny,1) 

!		print*, 'writing filter'
!		call filter_sgauss_r(dk, 0._dl, 1._dl, laplace, Fk, 1, planf, planb)
!		filt(:,:,1) = laplace(:,:,k_slice)
		
!		do i=1,nx; do j=1,ny
!			write(n_file, '(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, k_slice, filt(i,j,:) 
!		enddo; enddo
!		write(n_file, *)
!	end subroutine write_filter

end module zeta_mod







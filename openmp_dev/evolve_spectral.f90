!
! Module to store the desired evolution scheme for the fields, etc.
!

!
! To do : Add packing and unpacking of phase space variables
!

!
! To do: Add the Minkowski evolution and fixed background evolution
!   make sure I have the macros in the right place
!

! Add some macros to allow more seamless switching between numbers of dimensions

!#define THREADS 1
#define OMP 1
!#define USEMPI

!#define VECTOR(A,B) (/A,B/)
! Define the model parameters
!#define DVDPHI sin(fld(l,IRANGE))
!#define POTENTIAL -cos(fld(IRANGE))
!#define DVDPHI(PHI,PSI) VECTOR((PHI)**3 + g2*(PHI)*(PSI)**2, g2*(PHI)**2)
!#define POTENTIAL 0.25*fld(1,IRANGE)**4 + 0.5*g2*fld(1,IRANGE)**2*fld(2,IRANGE)**2
!#define M2I(PHI,PSI) VECTOR((PHI)**2+g2*(PSI)**2, g2*(PHI)**2)
!#define NUMFIELD 2

! to do: figure out field outputing, should be done using a binary file
! to do: make a precompiler statement for choosing how rho is remormalized
! to do: make a precompiler for how fields are initialized
! to do: make a separate module for output subroutines
! to do: remove redundant FFTs from the zeta filtering computation
! to do: make zeta calculation method macro

program lattice
#include "macros.h"
  use fftw3
  use params
  use hamiltonian
  use analysis
#ifdef OMP
  use omp_lib
#endif
	use zeta_mod
	!use correlator_mod
	use correlator_int_mod
	use sample_sites_mod
!	use renorm_mod

  implicit none
! Time Stepping Properties
  integer :: i, j, jj
  integer, parameter :: nstep = 2**9 			! Determines total number of output steps
  integer, parameter :: stepsize = 2**5		! Determines number of integration steps between outputing data
	integer, parameter :: stepslice = 2**2	! Determines number of output steps between outputing a lattice slice
	integer, parameter :: zslice = nz/2 		! Determines which slice of the lattice is output
  real(dl), parameter :: tstep0 = 0.8_dl/(2.0_dl**14)!1.0_dl/(2.0_dl**16)!dx/10000.!tstep = dx/10. ! Base time step
	real(dl) :: tstep = tstep0 							! time step which is updated
	real(dl) :: t_cur = 0.0 								! tracks the elapsed conformal time
!  real(dl), parameter :: tstep = 0.05
	integer :: seed_in											! seed used for initializing field fluctuations
	real(dl), parameter :: efold_init = 4.0_dl					! number of efolds of linear integration

  integer :: terror

  integer(kind=8) :: ti1,ti2,clock_rate
  real*8 :: tr1,tr2

  complex(C_DOUBLE_COMPLEX), pointer :: Fk2(:,:,:)
	complex(C_DOUBLE_COMPLEX), pointer :: Fk4(:,:,:) ! Fk3 is defined in zeta_mod.f90

	character(len=*), parameter :: run_ident = "_ending"		! Identifier appended to output file names

!	real(C_DOUBLE), pointer :: ggrad1(:,:,:)
!	real(C_DOUBLE), pointer :: ggrad2(:,:,:)
!	real(C_DOUBLE), pointer :: ggrad3(:,:,:)
!	real(C_DOUBLE), pointer :: ggrad4(:,:,:)
!	real(C_DOUBLE), pointer :: zf1(:,:,:)
!	real(C_DOUBLE), pointer :: zf2(:,:,:)
!	real(C_DOUBLE), pointer :: zfp1(:,:,:)
!	real(C_DOUBLE), pointer :: zfp2(:,:,:)
!	real(C_DOUBLE), pointer :: zlap1(:,:,:)
!	real(C_DOUBLE), pointer :: zlap2(:,:,:)

! This is currently just for debugging
  real(C_DOUBLE), dimension(nx,ny,nz) :: grad_squared_s, grad_squared_d

  print*,"Running with ", nstep*stepsize," steps and ",nstep," output steps, symp8"
  print*,"Lattice size is ", nx,ny,nz

! Begin by initializing FFTW
#ifdef THREADS
  terror = fftw_init_threads()
  call fftw_plan_with_nthreads(4)
#endif
#ifdef OMP
!  omp_set_num_threads(4)
  terror = fftw_init_threads()
  print*,"Error code is ",terror
  terror = omp_get_max_threads()
  print*,"Num Threads = ",terror
  call fftw_plan_with_nthreads(omp_get_max_threads())
#endif

! initialize arrays for doing FFT
  call cpu_time(tr1); call system_clock(ti1)
  call initialize_arrays()
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"FFTW setup :", dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1); call system_clock(ti1)
  call init_output()

! initialize fields	
	call init_cor_rad_cosmic(efold_init, H0)										! provides initial power spectrum by integrating mode functions on oversampled, isotropic
	! This loop is for testing field initialization, ordinarily only initialize fileds once
	!call init_cor_rad_cosmic_z(5.0_dl, dphi0, H0)					! provides initial power spectrum for fields and zeta
	!call write_cor_lat()
	do seed_in=2,2
	print*, "nfld=", nfld
	call init_fields_cor_rad(cor_rad, seed_in)						! initialize fields from oversampled, isotropic power spectrum
	!call zeta_init_zero()																	! initialize zeta to 0 on the whole lattice
	call zeta_init_pert([nx,ny,nz],dk,Fk,Fk2,Fk3,planf,planb)																	! initialize zeta to one-field perterbation theory result

	!call init_fields_cor2(cor_lat,seed_in)								! initializes fields to match a supplied power spectrum on lattice
	!call init_fields_cor_conv_test(cor_lat, seed_in)			! initialization for convergence test
	!call init_fields_ind(seed_in)												! initialize fields indepentently
	!enddo
  !call init_fields(1.)	! initialize fields as in Minkowski space
	call init_sites()	! randomly sample lattice sites at which to track trajectories
!	call init_sites_hom()
!	call read_sites()
	call init_sample_zeta()
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
	if (run_hom==.FALSE.) then
		call output_homogeneous()
	elseif (run_hom==.TRUE.) then
		call init_homogeneous()
	endif
  print*,"Field Initialization time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1); call system_clock(ti1)
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"Initial Output Time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1); call system_clock(ti1)
	!initialize dzeta
	call get_dzeta([nx,ny,nz],dk,Fk,Fk2,Fk3,planf,planb)					! call for unsmoothed zeta
	!call get_dzeta_part([nx,ny,nz],dk,Fk,Fk2,Fk3,planf,planb)		! call for unsmoothed zeta split by fields
	!call get_dzeta_part_kgv([nx,ny,nz],dk,Fk,Fk2,Fk3,planf,planb)! call for unsmoothed zeta split by K, G, K+V, G+V
	!call get_dzeta_part_kgv_smooth([nx,ny,nz], dk, 1._dl, Fk, Fk2, Fk3, planf, planb)	! call for unsmoothed and smoothed zeta split by K+V, G
	!call get_dzeta_smooth([nx,ny,nz],dk,Fk,Fk2,Fk3,planf,planb)	! call for smoothed zeta
	do i=1,n_scale
		!call zeta_filter_sharp_k(dk, k_hor_smooth(i), Fk, 2*i-1,planf, planb)		! call for zeta_smooth calculated with a low-pass sharp k filter
		call zeta_filter_sgauss_k(dk, k_hor_smooth(i), Fk, ord(i), i, planf, planb)		! call for zeta_smooth calculated with a low-pass s-Gaussian k filter
	enddo

	t_cur = 0.0_dl
	call output_run_param(88)
	call make_output(0._dl, 0)
  do j=1,nstep
    print*,"step ", j
    !call symp6(tstep, stepsize)
		call symp8(tstep, stepsize)
    !step zeta
		call zeta_step(stepsize*tstep, [nx,ny,nz], dk, Fk, Fk2, Fk3, planf, planb)	! step zeta_lat without calculation of zeta partials
		!call zeta_step_kgv(stepsize*tstep, [nx,ny,nz], dk, Fk, Fk2, Fk3, planf, planb)		! Step zeta integration with K, G, K+V, G+V parts
		!call zeta_step_kgv_smooth(stepsize*tstep, [nx,ny,nz], dk, 1._dl, Fk, Fk2, Fk3, planf, planb)		! Step zeta integration with K+V, G, K+V smooth, G smooth parts

		!call zeta_kick(dk, k2_kick_low, 1.0_dl, Fk, planf, planb)	! call for zeta_smooth calculated with kicks from a sharp k filter
		!call zeta_kick_pert_rho(dk, k2_kick_low, 3.0_dl, Fk, Fk2, Fk3, planf, planb)
		!call zeta_kick_pert(dk, k2_kick_low_pert, 1.0_dl, Fk, Fk2, Fk3, planf, planb)

		! Calculate filtered zeta
		do i=1,n_scale
			!call zeta_filter_sharp_k(dk, k_hor_smooth(i), Fk, 2*i-1,planf, planb)		! call for zeta_smooth calculated with a low-pass sharp k filter
			call zeta_filter_sgauss_k(dk, k_hor_smooth(i), Fk, ord(i), i, planf, planb)		! call for zeta_smooth calculated with a low-pass s-Gaussian k filter
		enddo

		!adapt step size
		t_cur = t_cur + stepsize*tstep
		call tstep_adapt()
    !call make_output(j*stepsize*tstep)
		!if (modulo(j,2**5)==0) then
		call make_output(t_cur, j)
		!endif
  enddo
	enddo
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"Total Evolution Time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  contains

    subroutine initialize_threads()

    end subroutine initialize_threads

    subroutine initialize_arrays()
#ifdef THREEDIM
  call allocate_fftw_array(nx,ny,nz,laplace,Fk)
  planf = fftw_plan_dft_r2c_3d(nz, ny, nx, laplace, Fk, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
  planb = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, laplace, FFTW_MEASURE+FFTW_DESTROY_INPUT) 

  call allocate_3d_fourier_array(nx,ny,nz,Fk2)
	call allocate_3d_fourier_array(nx,ny,nz,Fk3)
	call allocate_3d_fourier_array(nx,ny,nz,Fk4)
	call allocate_3d_real_array(nx,ny,nz,ggrad1)
	call allocate_3d_real_array(nx,ny,nz,zlap1)
	call allocate_3d_real_array(nx,ny,nz,zlap2)
#endif

#ifdef TWODIM
  call allocate_fftw_array(nx,ny,laplace,Fk)
  planf = fftw_plan_dft_r2c_2d(ny,nx,laplace,Fk,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  planb = fftw_plan_dft_c2r_2d(ny,nx,Fk,laplace,FFTW_PATIENT+FFTW_DESTROY_INPUT)
#endif

#ifdef ONEDIM
  call allocate_fftw_array(nx,laplace,Fk)
  planf = fftw_plan_dft_r2c_1d(nx,laplace,Fk,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  planb = fftw_plan_dft_c2r_1d(nx,Fk,laplace,FFTW_PATIENT+FFTW_DESTROY_INPUT)
#endif
    end subroutine initialize_arrays

    subroutine symp8(dt, nsteps)
      real*8 :: dt
      integer :: nsteps

      integer :: j
      real*8, parameter :: w1 = 0.74167036435061295344822780
      real*8, parameter :: w2 = -0.40910082580003159399730010
      real*8, parameter :: w3 = 0.19075471029623837995387626
      real*8, parameter :: w4 = -0.57386247111608226665638733
      real*8, parameter :: w5 = 0.29906418130365592384446354
      real*8, parameter :: w6 = 0.33462491824529818378495798
      real*8, parameter :: w7 = 0.31529309239676659663205666
      real*8, parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3+w4+w5+w6+w7)

      call Hamiltonian_Split(w7*dt/2._dl,1)

      do j=1,nsteps
         call symp_o2step(dt, w7, w6)
         call symp_o2step(dt, w6, w5)
         call symp_o2step(dt, w5, w4)
         call symp_o2step(dt, w4, w3)
         call symp_o2step(dt, w3, w2)
         call symp_o2step(dt, w2, w1)
         call symp_o2step(dt, w1, w0)
         call symp_o2step(dt, w0, w1)
         call symp_o2step(dt, w1, w2)
         call symp_o2step(dt, w2, w3)
         call symp_o2step(dt, w3, w4)
         call symp_o2step(dt, w4, w5)
         call symp_o2step(dt, w5, w6)
         call symp_o2step(dt, w6, w7)
         if (j.eq.nsteps) then
            call symp_o2step(dt, w7, 0._dl)
         else
            call symp_o2step(dt, w7, w7)
         endif
      enddo

    end subroutine symp8

    subroutine symp6(dt, nsteps)
      real(dl) :: dt
      integer :: nsteps

      real(dl), parameter :: w3=0.784513610477560_dl
      real(dl), parameter :: w2=0.235573213359357_dl
      real(dl), parameter :: w1=-1.177679984177887_dl
      real(dl), parameter :: w0=1._dl - 2._dl*(w1+w2+w3)

      integer :: j

      call Hamiltonian_Split(w3*dt/2._dl,1)
      do j=1,nsteps
         call symp_o2step(dt,w3,w2)
         call symp_o2step(dt,w2,w1)
         call symp_o2step(dt,w1,w0)
         call symp_o2step(dt,w0,w1)
         call symp_o2step(dt,w1,w2)
         call symp_o2step(dt,w2,w3)
         if (j==nsteps) then
            call symp_o2step(dt,w3,0._dl)
         else
            call symp_o2step(dt,w3,w3)
         endif
      enddo
    end subroutine symp6

    subroutine symp4(dt, nsteps)
      real(dl) :: dt
      integer :: nsteps

      real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
      real(dl), parameter :: w0 = 1._dl - 2._dl*w1

      integer :: j
      
      call Hamiltonian_Split(w1*dt/2._dl,1)
      do j=1,nsteps
         call symp_o2step(dt, w1,w0)
         call symp_o2step(dt,w0,w1)
         if (j==nsteps) then
            call symp_o2step(dt,w1,0._dl)
         else
            call symp_o2step(dt,w1,w1)
         endif
      enddo
    end subroutine symp4

    subroutine symp2(dt,nsteps)
      real(dl) :: dt
      integer :: nsteps

      integer :: j
      call Hamiltonian_Split(dt/2._dl, 1)
      do j=1,nsteps-1
         call symp_o2step(dt,1._dl,1._dl)
      enddo
      call symp_o2step(dt, 1._dl,0._dl)
    end subroutine symp2

    subroutine write_fields(time)
      integer :: i
      real(dl) :: time

      integer :: j,k
      real(dl) :: lap(1:nfld), grad_sq(1:nfld)
      j=ny/2; k=nz/2

      laplace = fld(1,IRANGE)
      call laplacian_3d(nx, ny, nz,laplace,Fk,dk, planf, planb)
      laplace = fld(1,IRANGE)
      call gradient_squared_3d_spectral([nx,ny,nz],laplace,Fk,Fk2,grad_squared_s,dk,planf,planb)

!      do j=2, ny-1
      do i=2,nx-1
         k=nz/2; j=ny/2
         lap = (1./dx**2/cc)*(STENCIL(c,LAPLACIAN))
         grad_sq = (1./dx**2/cc)*(STENCIL(c,GRAD2))
#ifdef THREEDIM
!         write(99,*) i*dx, j*dx, fld(:,i,j,k), fldp(:,i,j,k), laplace(i,j,k), lap(:)
         write(99,'(30(ES22.15,2x))') i*dx, j*dx, fld(:,i,j,k), fldp(:,i,j,k), grad_squared_s(i,j,k), grad_sq(1)
#endif
#ifdef TWODIM
         write(99,*) time, i*dx, fld(:,i,j), fldp(:,i,j)
#endif
#ifdef ONEDIM
         write(99,*) time, i*dx, fld(:,i), fldp(:,i)
#endif
      enddo; !write(99,*); enddo
      write(99,*)

    end subroutine write_fields

    subroutine init_fields(vparam)
      real(dl) :: vparam

      integer :: i, j
      real(dl) :: x, x0, y, y0
      real(dl) :: ptemp, gamma
      integer :: nseed
      integer, allocatable :: seed(:)

!      gamma = 1. / (1.+vparam**2)**0.5
!      x0 = len/2.
!      y0 = len/2.

!      do i=1,nx
!         x=i*dx-x0
!         ptemp = vparam*cosh(gamma*x)
!#ifdef THREEDIM
!         fld(:,i,:,:) = 4.*atan(1./ptemp)
!#endif
!#ifdef TWODIM
!         fld(:,i,:) = 4.*atan(1./ptemp)
!#endif
!#ifdef ONEDIM
!         fld(:,i) = 4.*atan(1./ptemp)
!#endif
!      enddo

      call random_seed(SIZE=nseed)
			print*, 'nseed = ', nseed !Testing seed
      allocate(seed(nseed))
      seed = 29*(/ (i-1, i=1,nseed) /)!37*(/ (i-1, i=1,nseed) /)
			print*, 'seed = ', seed	!Testing seed
      call random_seed(PUT=seed)
			print*, 'seed = ', seed	!Testing seed
      deallocate(seed)
      
      !do j=1,nfld
      !    call sample(-0.25, 3.*fld0(1)**2)! fix this term for effective mass
      !  	 fld(j,IRANGE) = fld0(j) + laplace
      !  	 call sample(0.25, 3.*fld0(1)**2)
      !  	 fldp(j,IRANGE) = dfld0(j) + laplace
      !enddo
#ifdef TWOFLD2
			! This initializing routine uses m2eff calculated for chi0 = 0, if chi0 != 0 the mass matrix is not diagonal.
			! Should check that for the fluctuations produced using the mean field of chi is valid.
			!	Initialize phi field
			if (infl_option==1) then
      	call sample(-0.25, 3.*fld0(1)**2)
      	fld(1,IRANGE) = fld0(1) + laplace
      	call sample(0.25, 3.*fld0(1)**2)
      	fldp(1,IRANGE) = dfld0(1) + laplace
			elseif (infl_option==2) then
				if (potential_option==0) then
      		call sample(-0.25, m2)
      		fld(1,IRANGE) = fld0(1) + laplace
      		call sample(0.25, m2)
      		fldp(1,IRANGE) = dfld0(1) + laplace
				elseif (potential_option==1) then
      		call sample(-0.25, m2+g2*fld0(2))
      		fld(1,IRANGE) = fld0(1) + laplace !temporarily taken out fluctuations
      		call sample(0.25, m2+g2*fld0(2))
      		fldp(1,IRANGE) = dfld0(1) + laplace !temporarily taken out fluctuations
				elseif (potential_option==4 .OR. potential_option==5) then ! assumes fld0(2) = 0
      		call sample(-0.25, m2)
      		fld(1,IRANGE) = fld0(1) + laplace !temporarily taken out fluctuations
      		call sample(0.25, m2)
      		fldp(1,IRANGE) = dfld0(1) + laplace !temporarily taken out fluctuations
				endif
			endif
      ! Initialize chi field
			if (potential_option==1) then
				call sample(-0.25, g2*(fld0(1)-phi_p)**2-beta2)
      	fld(2,IRANGE) = fld0(2) + laplace !temporarily taken out fluctuations
      	call sample(0.25, g2*(fld0(1)-phi_p)**2-beta2)
      	fldp(2,IRANGE) = dfld0(2) + laplace !temporarily taken out fluctuations
			elseif (potential_option==4) then
				call sample(-0.25, m2_inf - 6.*(m2_inf-m2_p) / (5.+dcosh(2.*dsqrt(3.*g2)*(fld0(1)-phi_p)/dsqrt(m2_inf-m2_p))) )
      	fld(2,IRANGE) = fld0(2) + laplace !temporarily taken out fluctuations
      	call sample(0.25, m2_inf - 6.*(m2_inf-m2_p) / (5.+dcosh(2.*dsqrt(3.*g2)*(fld0(1)-phi_p)/dsqrt(m2_inf-m2_p))))
      	fldp(2,IRANGE) = dfld0(2) + laplace !temporarily taken out fluctuations
			elseif (potential_option==5) then
				call sample(-0.25, m2_inf)
      	fld(2,IRANGE) = fld0(2) + laplace !temporarily taken out fluctuations
      	call sample(0.25, m2_inf)
      	fldp(2,IRANGE) = dfld0(2) + laplace !temporarily taken out fluctuations			
			endif
#endif
      yscl = 1.
      call calc_metric()
    end subroutine init_fields

    subroutine dump_rho(time)
      real(dl) :: time
      integer :: l

      real(dl) :: GE, PE, KE, rho, mom
      real(dl) :: elap, lap(nfld)
      integer :: i,j,k
      real(dl) :: acur, fac1,fac2
      
      acur = get_scale_factor()

! This loop can be shortened, as I'm double summing a lot of things
!      GE = 0.
!      GE = sum( (yvec(2:nlat)-yvec(1:nlat-1))**2 )
!      GE = GE + (yvec(nlat)-yvec(1))**2
!      GE = 0.5*GE / dx**2 / nlat  ! This is 0.5 not 0.25 since there is an overall factor of 2 on all the above stuff

#ifdef SPECTRAL
      GE = 0._dl
      do l=1,nfld
         laplace = fld(l,IRANGE)
#ifdef THREEDIM
         call laplacian_spectral(nx, ny, nz,laplace,Fk,dk, planf, planb)
#endif
#ifdef TWODIM
         call laplacian_spectral(nx, ny, laplace, Fk, dk, planf, planb)
#endif
#ifdef ONEDIM
         call laplacian_spectral(nx, laplace, Fk, dk, planf, planb)
#endif

#ifdef VECTORIZE
         GE =  GE - sum(fld(l,IRANGE)*laplace(IRANGE))
#endif
#ifdef LOOPEVOLVE
         FLOOP
           GE = GE - fld(l,LATIND)*laplace(LATIND)
         FLOOPEND  
#endif
      enddo
      GE = 0.5_dl*GE / nvol
#endif

#ifdef DISCRETE
      call wrap_fields()
      elap = 0.5/dx**2/cc
      GE = 0._dl
!$OMP PARALLEL DO PRIVATE(lap) REDUCTION(+:GE)
      FLOOP
        lap = STENCIL(c,LAPLACIAN)
        GE = GE - sum(fld(:,i,j,k)*lap(:))
      FLOOPEND
!$OMP END PARALLEL DO
      GE = elap*GE / nvol
#endif

#ifdef VECTORIZE
#ifdef ONEFLD
      PE = sum(potential_test(fld(1,IRANGE)))!sum(potential(fld(1,IRANGE),fld(2,IRANGE))) !sum(potential(fld(1,IRANGE))) !
#elif TWOFLD2
			PE = sum(potential_test(fld(1,IRANGE),fld(2,IRANGE)))!sum(potential(fld(1,IRANGE),fld(2,IRANGE))) !sum(potential(fld(1,IRANGE))) !
#endif
      KE = sum(fldp(:,IRANGE)**2)
#endif
#ifdef LOOPEVOLVE
      PE=0.
!$OMP PARALLEL DO REDUCTION(+:PE)
      FLOOP
        PE = PE + potential(fld(1,LATIND),fld(2,LATIND)) !potential(fld(1,LATIND)) !
      FLOOPEND  
!$OMP END PARALLEL DO
#endif 

!      KE = 0._dl
!      FLOOP
!        fac1 = 1.+zeta*fld(1,LATIND)**2
!        fac2 = fac1/(1.+zetafac*fld(1,LATIND)**2)
!        KE = KE + fac1*(fac2*fldp(1,LATIND)**2 + fldp(2,LATIND)**2)
!      FLOOPEND
!      KE = 0.5*KE*kinetic_norm() / nvol

#ifdef LOOPEVOLVE
      KE=0.
!$OMP PARALLEL DO REDUCTION(+:KE)      
      FLOOP
        KE = KE + sum(fldp(:,LATIND)**2)
      FLOOPEND
!$OMP END PARALLEL DO
#endif
      KE = 0.5_dl*KE * kinetic_norm() / nvol
      PE = PE / nvol

      GE = GE / acur**2
      
      rho = KE + PE + GE
! Evaluate momentum
!      laplace = fld(1:nlat,1:nlat)
!      call derivative_2d(nlat, lattice, Fk, dk)
!      mom = sum(yvec(nlat+1:2*nlat)*lattice(1:nlat))
!      mom = mom / nlat

#ifdef THREEDIM
#ifdef ONEFLD
      write(98,'(30(ES22.15,2X))') time, acur, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, -ysclp**2/12._dl/acur**4, &
           sum(fld(1,IRANGE))/nvol, sum(fldp(1,IRANGE))/nvol
#elif TWOFLD2
			write(98,'(30(ES22.15,2X))') time, acur, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, -ysclp**2/12._dl/acur**4, &
           sum(fld(1,IRANGE))/nvol, sum(fldp(1,IRANGE))/nvol, sum(fld(2,IRANGE))/nvol, sum(fldp(2,IRANGE))/nvol
#endif
#endif
#ifdef TWODIM
      write(98,*) time, rho, KE, PE, GE, grav_energy(), sum(fld(1,:,:))/nvol, sum(fld(2,:,:))/nvol
#endif
#ifdef ONEDIM
      write(98,*) time, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, sum(fld(1,:))/nvol, sum(fld(2,:))/nvol
#endif

    end subroutine dump_rho

!
! Randomly sample a gaussian random field with the appropriate spectrum
!
#define KCUT_FAC 0.9_dl
    subroutine sample(gamma, m2eff, spec)
!      real(C_DOUBLE), pointer :: f(:,:,:)
      real(dl) :: gamma
      real(dl) :: m2eff
      real(dl), optional :: spec
      
      type(C_PTR) :: plan_sin
      type(C_PTR) :: plan1
      
      integer, parameter :: os = 16, nos = max(nx,ny,nz)*os**2
      real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = KCUT_FAC*min(nnx,nny,nnz)*dk/2.0
      complex, parameter :: w = (0.0, twopi)

      real(dl) :: ker(nos), a(nnx), p(nnx)
      integer, allocatable :: seed(:)
      integer :: nseed

      integer :: i, j, k, l; 
      real(dl) :: kk
#ifdef THREEDIM
      real(dl), parameter :: norm = 0.5/(nvol*(twopi*dk**3)**0.5*mpl)*(dkos/dxos)
#endif
#ifdef TWODIM
      real(dl), parameter :: norm = 0.5/nvol/mpl*(dkos/dxos)  ! Work this one out
#endif
#ifdef ONEDIM
      real(dl), parameter :: norm = 0.5/nvol/mpl*(dkos/dxos)  ! Work this one out
#endif
      
    ! calculate (oversampled) radial profile of convolution kernel
      do k = 1,nos; kk = (k-0.5)*dkos
         ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-(kk/kcut)**2)	! filter fluctuations for Niquist
				 !ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-(kk/KCUT_FAC*(H0*exp(H0*(phi_p-phi0)/dphi0)))**2)	 !filter fluctuations for horizon scale
				 !ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-0.5*(kk/(KCUT_FAC*H0*exp(H0*(phi_p-phi0)/dphi0)*exp(H0/sqrt(sqrt(1.e5)*abs(dphi0)))))**2)	 !filter fluctuations at the horizon sclae at the end of non-adiatic event !!!FIX THIS WHEN YOU GET THE CHANCE FILTERING SHOULD NOT DEPEND ON g2, I'VE HARD CODED A VALUE AS A TEMPORARY FIX.
      end do
      
!      Assign the kernel here
!      plan = fftw_plan_r2r_1d(nos, ker, ker, FFTW_RODFT10, FFTW_ESTIMATE)
!      call fftw_execute
      
      plan_sin = fftw_plan_r2r_1d(nos,ker,ker,FFTW_RODFT10,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      call fftw_execute_r2r(plan_sin, ker, ker)
      call fftw_destroy_plan(plan_sin)

      do k = 1,nos; ker(k) = norm * ker(k)/k; end do
       ! initialize 3D convolution kernel (using linear interpolation of radial profile)
         FLOOP
#ifdef THREEDIM
            kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os
#endif
#ifdef TWODIM
            kk = sqrt(dble(i-nn)**2 + dble(j-nn)**2)
#endif
#ifdef ONEDIM
            kk = sqrt(dble(i-nn)**2)
#endif
            l = floor(kk)
            
            if (l > 0) then
               laplace(LATIND) = ker(l) + (kk-l)*(ker(l+1)-ker(l))
            else
#ifdef THREEDIM
               laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0
#endif
#ifdef TWODIM
               laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0  ! check this one
#endif
#ifdef ONEDIM
               laplace(LATIND) =
#endif
            end if
      FLOOPEND

       ! convolve kernel with delta-correlated Gaussian noise
!         plan1 = fftw_plan_dft_r2c_3d(nz,ny,nx,f,Fk, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, laplace, Fk)
!         call fftw_destroy_plan(plan1)

      do k=1,nz; do j=1,ny
         call random_number(a); call random_number(p)
         Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
      enddo; enddo

!         plan1 = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_c2r(planb, Fk, laplace)
!         call fftw_destroy_plan(plan1)

    end subroutine sample

    subroutine init_output()
      !open(unit=99,file="field_values_spec.out")
			open(unit=98,file="energy_spec"//run_ident//".out")
			open(unit=97,file="spectrum"//run_ident//".out")
			open(unit=96,file="zeta"//run_ident//".out")
!			open(unit=94,file="lat_sample"//run_ident//".out")
			!open(unit=92,file="lat_dump.out",form="unformatted") !TESTING REQUIRED! binary file of field output
!			open(unit=93,file="sample_sites.out")
			open(unit=91,file="init_spectrum"//run_ident//".out")
			open(unit=90,file="lat_slice"//run_ident//".out")
			!open(unit=89,file="lat.out")
			open(unit=88,file="run_param"//run_ident//".out")
			!open(unit=87,file='zeta_smooth.out')
			open(unit=86,file="zeta_part"//run_ident//".out")
			open(unit=85,file="spectrum_zeta_part"//run_ident//".out")
			!open(unit=84,file='filter_test.out') ! done instead in post processing
#ifdef THREEDIM
      call init_spectrum_3d(nx,ny,nz)
#endif
#ifdef TWODIM
      call init_spectrum_2d(nx,ny)
#endif
#ifdef ONEDIM
      call init_spectrum_1d(nx)
#endif
    end subroutine init_output
		
		! Subroutine to output the run parameters in a human readable way
		! to do: make output for how zeta is calculated/smoothed
		! to do: make output identification for integration scheme
		subroutine output_run_param(n_file)
			integer :: n_file		! file number
			! write lattice parameters
			write(n_file, '(A, 2x, 3(I5,2x))'), "nx, ny, nz: ", nx, ny, nz
			write(n_file, '(A, 2x, ES22.15, 2x)'), "len: ", len
			write(n_file, '(A, 2x, ES22.15, 2x)'), "mpl: ", mpl
			! write run parameters
			write(n_file, '(A, 2x, 3(I5,2x))'), "nstep, stepsize, stepslice: ", nstep, stepsize, stepslice
			write(n_file, '(A, 2x, ES22.15, 2x)'), "tstep0: ", tstep0
			! write mean initial conditions
			if (nfld==1) then
				write(n_file, '(A, 2x, 3(ES22.15, 2x))'), "phi0, dphi0, H0: ", phi0, dphi0, H0
			elseif (nfld==2) then
				write(n_file, '(A, 2x, 5(ES22.15, 2x))'), "phi0, dphi0, chi0, dchi0, H0: ", phi0, dphi0, chi0, dchi0, H0
			endif
			! write fluctuation initialization parameters
			write(n_file, '(A, 2x, I5, 2x)'), "seed_in", seed_in
			write(n_file, '(A, 2x, 2(ES22.15, 2x))'), "efold_init, nstep_mode: ", efold_init, nstep_mode
			write(n_file, '(A, 2x, (ES22.15, 2x))'), "KCUT_FAC: ", KCUT_FAC 
			! write potential type and parameters
			write(n_file, '(A, 2x, I5,2x)'), "infl_option: ", infl_option
			if (nfld == 2) then
				write(n_file, '(A, 2x, I5,2x)'), "potential_option: ", potential_option
			endif
			if (infl_option == 2) then
				write(n_file, '(A, 2x, (ES22.15, 2x))'), "m2", m2
			endif
			if (potential_option == 5 .and. nfld == 2) then
				write(n_file, '(A, 2x, (ES22.15, 2x))'), "lambda_chi: ", lambda_chi 
				write(n_file, '(A, 2x, (ES22.15, 2x))'), "m2_inf: ", m2_inf
			elseif (potential_option == 6 .and. nfld == 2) then
				write(n_file, '(A, 2x, (ES22.15, 2x))'), "lambda_chi: ", lambda_chi
				write(n_file, '(A, 2x, 4(ES22.15, 2x))'), "phi_p, phi_w, m2_inf, m2_p: ", phi_p, phi_w, m2_inf, m2_p
			endif
			! write zeta calculation proceedure, ie what are the partial zetas

			! write zeta smoothing proccedure and scales
			write(n_file, '(A, 2x, 2(I5, 2x))'), "n_scale, n_moment: ", n_scale, n_moment
			write(n_file, '(A, 2x, 30(ES22.15, 2x))'), "k_hor_smooth: ", k_hor_smooth
			write(n_file, '(A, 2x, 30(I5, 2x))'), "ord: ", ord
			close(unit=n_file)
		end subroutine output_run_param

    subroutine make_output(time, step)
      real(dl) :: time
			integer :: step
      integer :: i
      !real(dl) :: spec(ns,2*nfld+2)
			!real(dl) :: spec(ns,(nfld*2)**2+4)  !to spectum and cross spec for both fields and zeta spectrum
			!real(dl) :: spec(ns,2*(4*nfld+4)) !to spectrum and cross spec for both fields and zeta and difference 

      call dump_rho(time)
			call write_zeta(time, 96)
			!call write_lat_sample(time)
			if (0 == mod(step, stepslice)) then
				call write_spec(97)
				call write_slice(time, zslice, 90)
				call write_zeta_partial(time, zslice, 86)
				call write_zeta_part_spec(85)
			endif
			!if (step==nstep .or. step==80) then
			if (step==nstep) then
				!call write_lat(time)
			endif
! if the WINT is not defined call lat_dump to output the lattice to a binary file 
!#ifndef WINT
!			call lat_dump()
!#endif 

!      laplace(IRANGE) = fld(2,IRANGE) !commented out due to being for two field model
!#ifdef THREEDIM
      !laplace(IRANGE) = fld(1,IRANGE)
      !call spectrum_3d(spec(:,1),laplace, Fk, planf)
      !Fk2=Fk
      !laplace(IRANGE) = fldp(1,IRANGE)
      !call spectrum_3d(spec(:,2), laplace, Fk, planf)
      !call crossspec_3d(Fk, Fk2, spec(:,3),spec(:,4)) !check ordering for real vs imaginary parts
! one field zeta spec output
			!laplace(IRANGE) = zeta_lat(IRANGE)
			!call spectrum_3d(spec(:,5), laplace, Fk, planf)
			!Fk2=Fk
			!laplace(IRANGE) = dzeta_lat(2,IRANGE)
			!call spectrum_3d(spec(:,6), laplace, Fk, planf)
			!call crossspec_3d(Fk, Fk2, spec(:,7), spec(:,8))
! two field model output
			!laplace(IRANGE) = fld(2,IRANGE)
			!call spectrum_3d(spec(:,5),laplace, Fk3, planf)
      !Fk4=Fk3
      !laplace(IRANGE) = fldp(2,IRANGE)
      !call spectrum_3d(spec(:,6), laplace, Fk3, planf)
			!call crossspec_3d(Fk3, Fk4, spec(:,7), spec(:,8))
! two field cross spectrum output
			! crossspec takes ft'ed inputs, can avoid doing the same fft twice
			! Fk - dphi
			! Fk2 - phi
			! Fk3 - dchi
			! Fk4 - chi
			!call crossspec_3d(Fk2, Fk4, spec(:,13), spec(:,14)) !phi-chi
			!call crossspec_3d(Fk2, Fk3, spec(:,15), spec(:,16)) !phi-dchi
			!call crossspec_3d(Fk, Fk4, spec(:,17), spec(:,18)) !dphi-chi
			!call crossspec_3d(Fk, Fk3, spec(:,19), spec(:,20)) !dphi-dchi
! zeta spec output
			!laplace(IRANGE) = zeta_lat(IRANGE)
			!call spectrum_3d(spec(:,9), laplace, Fk, planf)
			!Fk2=Fk
			!laplace(IRANGE) = dzeta_lat(2,IRANGE)
			!call spectrum_3d(spec(:,10), laplace, Fk, planf)
			!call crossspec_3d(Fk, Fk2, spec(:,11), spec(:,12))

! zeta spec output
!			laplace(IRANGE) = zeta_lat(IRANGE)
!			call spectrum_3d(spec(:,1), laplace, Fk, planf)
!			Fk2=Fk
!			laplace(IRANGE) = dzeta_lat(IRANGE)	!laplace(IRANGE) = dzeta_lat(2,IRANGE)
!			call spectrum_3d(spec(:,2), laplace, Fk, planf)
!			call crossspec_3d(Fk, Fk2, spec(:,3), spec(:,4))
! field spec output
!			do i=1,1
!      	laplace(IRANGE) = fld(i,IRANGE)
!      	call spectrum_3d(spec(:,4*i+1),laplace, Fk, planf)
!      	Fk2=Fk
!      	laplace(IRANGE) = fldp(i,IRANGE)
!      	call spectrum_3d(spec(:,4*i+2), laplace, Fk, planf)
!      	call crossspec_3d(Fk, Fk2, spec(:,4*i+3),spec(:,4*i+4)) !check ordering for real vs imaginary parts
!			enddo
!#ifdef TWOFLD2
! two field model output
!			laplace(IRANGE) = fld(2,IRANGE)
!			call spectrum_3d(spec(:,9),laplace, Fk3, planf)
!      Fk4=Fk3
!      laplace(IRANGE) = fldp(2,IRANGE)
!      call spectrum_3d(spec(:,10), laplace, Fk3, planf)
!			call crossspec_3d(Fk3, Fk4, spec(:,11), spec(:,12))
! two field cross spectrum output
			! crossspec takes ft'ed inputs, can avoid doing the same fft twice
			! Fk - dphi
			! Fk2 - phi
			! Fk3 - dchi
			! Fk4 - chi
!			call crossspec_3d(Fk2, Fk4, spec(:,13), spec(:,14)) !phi-chi
!			call crossspec_3d(Fk2, Fk3, spec(:,15), spec(:,16)) !phi-dchi
!			call crossspec_3d(Fk, Fk4, spec(:,17), spec(:,18)) !dphi-chi
!			call crossspec_3d(Fk, Fk3, spec(:,19), spec(:,20)) !dphi-dchi			 	
!#endif

! Spectrum of difference
! Untested, requires two modes of running, can have two verisions of this function that are toggled by either a precompiler
! statement, or from potential_mod potential_option parameter and use the opposite toggle for calling or not calling lat_dump.
			! call an instance of lat_read
!#ifdef WINT
!			call lat_read(f, fp, z, dz) ! pass varibles, need declaration/renaming

			! comupte spectra
!			laplace(IRANGE) = fld(1,IRANGE) - f(1,IRANGE)
!			call spectrum_3d(spec(:,13), laplace, Fk, planf)
!			Fk2=Fk
!      laplace(IRANGE) = fldp(1,IRANGE) - fp(1,IRANGE)
!      call spectrum_3d(spec(:,14), laplace, Fk, planf)
!      call crossspec_3d(Fk, Fk2, spec(:,15),spec(:,16))
! two field model output
!			laplace(IRANGE) = fld(2,IRANGE) - f(2,IRANGE)
!			call spectrum_3d(spec(:,17),laplace, Fk, planf)
!      Fk2=Fk
!      laplace(IRANGE) = fldp(2,IRANGE) - fp(2,IRANGE)
!      call spectrum_3d(spec(:,18), laplace, Fk, planf)
!			call crossspec_3d(Fk, Fk2, spec(:,19), spec(:,20))
! zeta spec output
!			laplace(IRANGE) = zeta_lat(IRANGE) - z(IRANGE)
!			call spectrum_3d(spec(:,21), laplace, Fk, planf)
!			Fk2=Fk
!			laplace(IRANGE) = dzeta_lat(2,IRANGE) - dz(IRANGE)
!			call spectrum_3d(spec(:,22), laplace, Fk, planf)
!			call crossspec_3d(Fk, Fk2, spec(:,23), spec(:,24))
!#endif
!#endif
!#ifdef TWODIM
!      call spectrum_2d(spec, laplace, Fk, planf)
!#endif
!#ifdef ONEDIM
!      call spectrum_1d(spec, laplace, Fk, planf)
!#endif
!      do i=1,ns
!         write(97,'(30(ES22.15,2x))') (i-1)*dk, spec(i,:)
!      enddo
!      write(97,*)
      
    end subroutine make_output

		! Subroutine that outputs the spectra and cross-spectra of the component parts of zeta
		! to do: calculate spectra
		! to do: calculate cross-spectra
		! to do: make output
		subroutine write_zeta_part_spec(n_file)
			integer :: n_file		! file number
			integer :: i,j
			real(dl) :: spec(ns,2*(nzpart) + 2)

			do i=1,nzpart
				laplace(IRANGE) = zeta_part(i,IRANGE)
				call spectrum_3d(spec(:,2*i-1), laplace, Fk, planf)
				laplace(IRANGE) = dzeta_part(i,IRANGE)	!dzeta_part(2,2*(i-1)+j,IRANGE)
				call spectrum_3d(spec(:,2*i), laplace, Fk, planf)
			enddo
			do i=1,2
				laplace(IRANGE) = epsilon_part(i,IRANGE)
				call spectrum_3d(spec(:,i + 8), laplace, Fk, planf)
			enddo

			do i=1,ns
         write(n_file,'(30(ES22.15,2x))') (i-1)*dk, spec(i,:)
      enddo
      write(n_file,*)	
		end subroutine write_zeta_part_spec

		! Subroutine to calculate and output the spectra and cross spectra of the fields and zeta
		! Formatting: 
		subroutine write_spec(n_file)
			integer :: n_file		! file number
			real(dl) :: spec(ns,(nfld*2)**2+4)  !to spectum and cross spec for both fields and zeta spectrum
			integer :: i				! indexing variable
#ifdef THREEDIM
! zeta spec output
			laplace(IRANGE) = zeta_lat(IRANGE)
			call spectrum_3d(spec(:,1), laplace, Fk, planf)
			Fk2=Fk
			laplace(IRANGE) = dzeta_lat(IRANGE)	!laplace(IRANGE) = dzeta_lat(2,IRANGE)
			call spectrum_3d(spec(:,2), laplace, Fk, planf)
			call crossspec_3d(Fk, Fk2, spec(:,3), spec(:,4))
			! one field spec output
			do i=1,1
      	laplace(IRANGE) = fld(i,IRANGE)
      	call spectrum_3d(spec(:,4*i+1),laplace, Fk, planf)
      	Fk2=Fk
      	laplace(IRANGE) = fldp(i,IRANGE)
      	call spectrum_3d(spec(:,4*i+2), laplace, Fk, planf)
      	call crossspec_3d(Fk, Fk2, spec(:,4*i+3),spec(:,4*i+4)) !check ordering for real vs imaginary parts
			enddo
#ifdef TWOFLD2
			! two field model output
			laplace(IRANGE) = fld(2,IRANGE)
			call spectrum_3d(spec(:,9),laplace, Fk3, planf)
      Fk4=Fk3
      laplace(IRANGE) = fldp(2,IRANGE)
      call spectrum_3d(spec(:,10), laplace, Fk3, planf)
			call crossspec_3d(Fk3, Fk4, spec(:,11), spec(:,12))
			! two field cross spectrum output
			! crossspec takes ft'ed inputs, can avoid doing the same fft twice
			! Fk - dphi, Fk2 - phi, Fk3 - dchi, Fk4 - chi
			call crossspec_3d(Fk2, Fk4, spec(:,13), spec(:,14)) !phi-chi
			call crossspec_3d(Fk2, Fk3, spec(:,15), spec(:,16)) !phi-dchi
			call crossspec_3d(Fk, Fk4, spec(:,17), spec(:,18)) !dphi-chi
			call crossspec_3d(Fk, Fk3, spec(:,19), spec(:,20)) !dphi-dchi			 	
#endif
#endif
#ifdef TWODIM
      call spectrum_2d(spec, laplace, Fk, planf)
#endif
#ifdef ONEDIM
      call spectrum_1d(spec, laplace, Fk, planf)
#endif
      do i=1,ns
         write(n_file,'(30(ES22.15,2x))') (i-1)*dk, spec(i,:)
      enddo
      write(n_file,*)
      
		end subroutine write_spec

		! Output homogeneous intitial conditions to file
		subroutine output_homogeneous()
			real(dl), dimension(nfld) :: fld_hom
			real(dl), dimension(nfld) :: fldp_hom
			integer :: i

			do i=1,nfld
				fld_hom(i) = sum(fld(i,IRANGE))/dble(nx)/dble(ny)/dble(nz)
				fldp_hom(i) = sum(fldp(i,IRANGE))/dble(nx)/dble(ny)/dble(nz)
			enddo
			open(unit=95,file="hom_ic.out")
			do i=1,nfld
				write(95,'(30(ES22.15,2x))') fld_hom(i)
				write(95,'(30(ES22.15,2x))') fldp_hom(i)
			enddo
			close(95)
		end subroutine output_homogeneous

		! Initialize homogeneous fields from a file
		subroutine init_homogeneous()
			real(dl), dimension(nfld) :: fld_hom
			real(dl), dimension(nfld) :: fldp_hom
			integer :: i

			open(unit=95, file="hom_ic.out")
			do i=1,nfld
				read(95,*) fld_hom(i)
				read(95,*) fldp_hom(i)
			enddo
			close(95)
			do i=1,nfld
				fld(i,IRANGE) = fld_hom(i)
				fldp(i,IRANGE) = fldp_hom(i)
			enddo
		end subroutine init_homogeneous

	! Write lattice to an unformatted file
	! Still untested
	subroutine lat_dump()
		integer :: i
		do i=1, nfld
			! write field and field momentum n the whole lattice to a file
			write(92) fld(i,IRANGE)
			write(92) fldp(i,IRANGE)
		enddo
		! write zeta and dzeta/dtau to a file for whole lattice
		write(92) zeta_lat(IRANGE)
		write(92) dzeta_lat(IRANGE)	!dzeta_lat(2,IRANGE) ! the index 2 is because dzeta_lat(1,IRANGE) hold dzeta_lat from the previous time step
	end subroutine lat_dump

	! Read in lattice from previous run and store value in some arrays
	! Still untested
	! Need to check if * in formating will correctly identify format from binary
	subroutine lat_read(f, fp, z, dz)
		! Declare variables
		real(dl), dimension(nfld,IRANGE) :: f, fp
		real(dl), dimension(IRANGE) :: z, dz
		integer :: i
		
		do i=1,nfld
			read(92,*) f(i,IRANGE)
			read(92,*) fp(i,IRANGE)
		enddo
		read(92,*) z(IRANGE)
		read(92,*) dz(IRANGE)
	end subroutine lat_read

	! Subroutine to write a slice of the lattice to an output file
	subroutine write_slice(time, kslice, n_file)
		real(dl), intent(in) :: time	! conformal time of output
		integer :: kslice							! the z index of the slice to be output
		integer :: i,j								! lattice index (x-y slice)
		integer :: n_file							! file number

		do i=1,nx; do j=1,ny
#ifdef ONEFLD
			!write(90,'(30(ES22.15,2x))') time, i, j, kslice, zeta_lat(i,j,kslice), dzeta_lat(2,i,j,kslice), fld(1,i,j,kslice), fldp(1,i,j,kslice)
			write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, zeta_lat(i,j,kslice), dzeta_lat(i,j,kslice), fld(1,i,j,kslice), fldp(1,i,j,kslice), zeta_smooth(1,i,j,kslice),  dzeta_smooth(i,j,kslice), zeta_smooth_pert(i,j,kslice), zeta_pert(i,j,kslice), delta_alpha(i,j,kslice)
#elif TWOFLD2
			!write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, zeta_lat(i,j,kslice), dzeta_lat(i,j,kslice), fld(1,i,j,kslice), fldp(1,i,j,kslice), fld(2,i,j,kslice), fldp(2,i,j,kslice), zeta_smooth(i,j,kslice),  dzeta_smooth(i,j,kslice)
			write(n_file,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, kslice, zeta_lat(i,j,kslice), dzeta_lat(i,j,kslice), fld(1,i,j,kslice), fldp(1,i,j,kslice), fld(2,i,j,kslice), fldp(2,i,j,kslice), zeta_smooth(:,i,j,kslice)
#endif
		enddo; enddo
		write(90,*)
	end subroutine write_slice

	! Subroutine to write lattice to an output file
	subroutine write_lat(time)
		real(dl), intent(in) :: time	! conformal time of output
		integer :: i,j,k							! lattice index
		
		do i=1,nx; do j=1,ny; do k=1,nz
#ifdef ONEFLD
			!write(90,'(30(ES22.15,2x))') time, i, j, kslice, zeta_lat(i,j,kslice), dzeta_lat(2,i,j,kslice), fld(1,i,j,kslice), fldp(1,i,j,kslice)
			write(89,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, i, j, k, zeta_lat(i,j,k), dzeta_lat(i,j,k), fld(1,i,j,k), fldp(1,i,j,k),&
			zeta_smooth(1,i,j,k),  dzeta_smooth(i,j,k)
#elif TWOFLD2
			write(89,'(ES22.15,2x,3(I5,2x),30(ES22.15,2x))') time, zeta_lat(i,j,k), dzeta_lat(i,j,k), zeta_part(1,i,j,k), zeta_part(2,i,j,k),&
			fld(1,i,j,k), fldp(1,i,j,k), fld(2,i,j,k), fldp(2,i,j,k), zeta_smooth(1,i,j,k),  dzeta_smooth(i,j,k)
#endif
		enddo; enddo; enddo
		write(89,*)		

	end subroutine write_lat

	subroutine tstep_adapt()
		tstep = tstep0/yscl
	end subroutine tstep_adapt

	! Subroutine to take a radial profile of power spectrum and interpolate and use it to interpolate
	! the power for the modes on the lattice
!	subroutine interp_cor(cor)

!	end subroutine interp_cor

	! Subroutine which takes the 2-point correlation functions of fields and field momenta as input
	! and generates random fields to match the 2-point correlation functions given.
	! This will replace the subroutine init_fields
	! The 2-point correlators need to be calculated separately, eg by integrating the mode functions
	! to do: option for correct real-space/ Fourier-space 2-point correlation
	! to do: set os (oversampling factor) as a precompiler so it can easily be passed around. Moved to correlator_mod
	! to do: apply filter for Niquist frequency
	! to do: the interpolation is interpolating the 2-point correlation, should maybe not interpolate on the square
	subroutine init_fields_cor(cor)
		!complex(dl), dimension(2*nfld,2*nfld,:) :: corr					! matrix of power spectrum for each (oversampled) mode
		!complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld,:) :: cor(:,:,:)
		real(dl), dimension(2*nfld,2*nfld,:) :: cor(:,:,:)

    integer, parameter :: os = 16, nos = max(nx,ny,nz)*os**2	! over-sampling factor and new number of samples
    real, parameter :: dxos = dx/os, dkos = dk/(2*os)					! oversampled lattice spacing and mode spacing
    complex, parameter :: w = (0._dl, twopi)									! i*2*pi, used for random phase

    !real(dl) :: ker(nos), a(nnx), p(nnx)											! kernal, amplitude, phase
		real(dl) :: ker(nos), a(2*nfld), p(2*nfld)								! kernal, amplitude, phase
    integer, allocatable :: seed(:)														! rng variable
    integer :: nseed																					! rng variable

    integer :: i, j, k, n, m, l																! lattice location indicies/ matrix indicies, rounded radial distance
    real(dl) :: kk																						! radial distance (real/momentum space)
#ifdef THREEDIM
    !real(dl), parameter :: norm = 0.5/(nvol*(twopi*dk**3)**0.5*mpl)*(dkos/dxos) ! double check that this is the correct value
		real(dl), parameter :: norm = 2._dl!*twopi*dkos / (dxos**2*nos*nvol) ! this is wrong
#endif
		type(C_PTR) :: plan_sin																		! FFT plan
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld,IRANGE) :: cor_interp							! Power spec interpolated to each Fourier mode on the lattice
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: cor_fac												! factored power at a particular k mode							
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv																		! vector of Gaussian random variables
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: conv																	! vector of random variables matching given power spectrum
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: Fk_cor											! vector to hold DFT of fld1, fldp1, fld2, fldp2, ...

		!!!!! Matching the power spectum !!!!!
		! Interpolate correlation matrix
		! make sine transform plan (used in interpolation)
		plan_sin = fftw_plan_r2r_1d(nos,ker,ker,FFTW_RODFT10,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))

		do n = 1,2*nfld; do m = n,2*nfld ! loop over correlation matrix entries
			do k = 1,nos; kk = (dble(k)-0.5)*dkos	! the -0.5 is due to the definition of the sine transform (FFTW_RODFT10)
				ker(k) = kk*cor(n,m,k) * exp(-0.5*(kk/(KCUT_FAC*H0))**2)	! filter here, knows nothing on Niquist frequency
			enddo
			call fftw_execute_r2r(plan_sin, ker, ker)	! transform to real space (radial)
			do k = 1,nos; ker(k) = norm * ker(k)/k; enddo	! normalize transform and insert 1/r factor

		  FLOOP
#ifdef THREEDIM
      	kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os	! calculate radial distance to the point (i,j,k)
#endif
#ifdef TWODIM
        kk = sqrt(dble(i-nn)**2 + dble(j-nn)**2)
#endif
#ifdef ONEDIM
        kk = sqrt(dble(i-nn)**2)
#endif
        l = floor(kk)	! round down radial distance to nearest sampled point
            
        if (l > 0) then
					!corr_interp(n,m,LATIND) = ker(l) + (kk-l)*(ker(l+1)-ker(l))	! linear interpolation between nearest sample points
        	laplace(LATIND) = ker(l) + (kk-l)*(ker(l+1)-ker(l))	! linear interpolation between nearest sample points
        else
#ifdef THREEDIM
					!corr_interp(n,n,LATIND) = (4.0*ker(1)-ker(2))/3.0
          laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0
#endif
#ifdef TWODIM
          laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0  ! check this one
#endif
#ifdef ONEDIM
          laplace(LATIND) =
#endif
            end if
      FLOOPEND
			call fftw_execute_dft_r2c(planf, laplace, Fk)
			cor_interp(n,m,:,:,:) = Fk(:,:,:)
		enddo; enddo
		
		! Initialize random number generator
		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
    seed = 29*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)

		! Cholesky factorization of corr_interp
		do i=1,nx; do j=1,ny; do k=1,nnz	! check ordering
			do m = 1,2*nfld; do n = 1,2*nfld
				if (n>=m) then
					cor_fac(m,n) = cor_interp(m,n,LATIND)
				else
					cor_fac(m,n) = CONJG(cor_interp(n,m,LATIND))! Take Hermitian conjugate here
				endif
			enddo; enddo
			print*, "i,j,k = ", i, j, k	! testing
			do m=1,2*nfld	! testing
				print*, cor_fac(m,:)
			enddo
			call zpotrf('L',2*nfld,cor_fac,2*nfld,l)
			if (l /= 0) then
				print*, "Factorization warning: l = ", l
			endif
			! Initialize GRV vector
			call random_number(a); call random_number(p)	! a and p are length 2*nfld
			grv(:) =	sqrt(-2.0*log(a)) * exp(w*p)
			call zgemv('C', 2*nfld, 2*nfld, (1._dl,0._dl), cor_fac, 2*nfld, grv, 1, (0._dl,0._dl), conv, 1)
			Fk_cor(:,LATIND) = conv(:)
		enddo; enddo; enddo
		
		! Peform FFT to calculate fluctuations in real space
		! Add fluctuations to homogeneous fields
		! Check if /nvol is required
		do m=1, nfld
			Fk(:,:,:) = Fk_cor(2*m-1,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)!!!!! fft to real space
			fld(m,IRANGE) = fld0(m) + laplace(IRANGE)
			Fk(:,:,:) = Fk_cor(2*m,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)
			fldp(m,IRANGE) = dfld0(m) + laplace(IRANGE)
		enddo
			
		call fftw_destroy_plan(plan_sin)	! destroy FFT plan		

		!!!!!!!! momentum space 2-point !!!!!!!!

		! make FFT plan
!		plan_sin = fftw_plan_r2r_1d(nos,ker,ker,FFTW_RODFT10,ior(FFTW_ESTIMATE,FFTW_UNALIGNED)) ! this is named plan_sin since it is a sine transform for real space 2-point correlator matching
		! Perform interpolation on correlation matrix
!		do n = 1,2*nfld; do m = n,2*nfld ! loop over correlation matrix entries
!			do k = 1,nos
!				ker(k) = corr(n,m,k)
!			enddo
			! Transform to real space and interpolate
!    	call fftw_execute_r2r(plan_sin, ker, ker)
!			ker = ker/nvol	! normalization factor for FFT
			! Here is the interpolation part
!			FLOOP
!#ifdef THREEDIM
!      	kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os	! calculate radial distance in units of dxos, this misplaces the center.
!#endif
!				l = floor(kk)	! calculate radial distance to next lowest quantity for which ker is calculated
!				if (l > 0) then
!        	laplace(LATIND) = ker(l) + (kk-l)*(ker(l+1)-ker(l))	! linear interpolation
!        else
!#ifdef THREEDIM
          !laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0	! I don't know where this part comes from, need to figure that out and if it is correct here
!					laplace(LATIND) = ker(1) - (kk-l)*(ker(2)-ker(1)) ! this is a linear extrapolation
!#endif
!			FLOOPEND
			! Interpolation is done, real space result stored in laplace.

			! Calculate Fourier space result
!			call fftw_execute_dft_r2c(planf, laplace, Fk)
			! Figure a way to store Fk so that the next correlation matrix entry can be calculated
!			corr_interp(n,m,:,:,:) = Fk(:,:,:)
			! Do the same using for Hermetian conjugate
			
!		enddo; enddo
		! Remember there is a particular way entries are stored in r2c transforms

		! Initialize random number generator
!		call random_seed(SIZE=nseed)
!		print*, 'nseed = ', nseed
!    allocate(seed(nseed))
!    seed = 29*(/ (i-1, i=1,nseed) /)
!    call random_seed(PUT=seed)
!    deallocate(seed)

		! looping over each mode
!		do k=1,nz; do j=1,ny; do i=1,nnx
				! Decompose correlation matrix for each mode (need to declare decomposition matrix)

				! Initialize GRV vector
!				call random_number(a); call random_number(p)	! a and p are length 2*nfld
!				grv(:) =	sqrt(-2.0*log(a)) * exp(w*p)
				! Multiply GRV vector by decomposed correlation matrix mode by mode
				! Result is a vector of length 2*nfld with Fourier space information for a single mode
				! Store result for FFT
!				Fk_cor(:,LATIND) = grv(:)
!    enddo; enddo; enddo
		
!		do n=1, 2*nfld, 2
			! Perform FFT on transformed GRV and add results to background field
!			Fk(:,:,:) = Fk_cor(n,:,:,:)
!			call fftw_execute_dft_c2r(planb, Fk, laplace)
!			fld(n,IRANGE) = fld0(n) + laplace/nvol
			! Repeat for momentum
!			Fk(:,:,:) = Fk_cor(n+1,:,:,:)
!			call fftw_execute_dft_c2r(planb, Fk, laplace)
!			fldp(n,IRANGE) = dfld0(n) + laplace/nvol
!		enddo	
		
		! loop over field indicies
		!do n=1,2*nfld
			! loop over modes
		!	do j=1,ny; do k=1,nz
		!		call random_number(a); call random_number(p) ! a and p are length nnx, so random number generates an array of length nnx
				! Check normalization on Gaussian
		!		norm_grv(n,:,j,k) =	sqrt(-2.0*log(a)) * exp(w*p)! the i component only goes to i=nnx do to the r2c transform reality condition
		!	enddo; enddo
		!enddo

		! matrix multiplication on norm_grv ()
	end subroutine init_fields_cor

	! Subroutine which takes the 2-point correlation functions of fields and field momenta as input
	! and generates random fields to match the 2-point correlation functions given.
	! This will replace the subroutine init_fields
	! The 2-point correlators need to be calculated separately, eg by integrating the mode functions
	! to do: option for correct real-space/ Fourier-space 2-point correlation
	! to do: 
	! to do: should make exception for k=0 mode to avoid factorization error
#define KCUT_FAC_H 16._dl
	subroutine init_fields_cor2(cor_in, seed_in)
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld,nnx,ny,nz) :: cor_in			! Input power/cross spectrum for each k on the lattice
		integer :: seed_in

    complex, parameter :: w = (0._dl, twopi)									! i*2*pi, used for random phase
		real(dl) :: a(2*nfld), p(2*nfld)													! amplitude, phase
    integer, allocatable :: seed(:)														! rng variable
    integer :: nseed																					! rng variable

    integer :: i, j, k, n, m, l																! lattice location indicies/ matrix indicies, check integer
		integer :: ii,jj,kk																				! mode numbers used for filtering
		real(dl), parameter :: kcut = KCUT_FAC*(min(nnx,nny,nnz)-1)!sqrt(2.*KCUT_FAC*(min(nnx,nny,nnz)-1))	!	frequency filtering at Nyquist
		real(dl), parameter :: kcut_h = KCUT_FAC_H*H0									! frequency filtering at horizon

#ifdef THREEDIM
    !real(dl), parameter :: norm = 0.5/(nvol*(twopi*dk**3)**0.5*mpl)*(dkos/dxos) ! double check that this is the correct value
		!real(dl), parameter :: norm = 2._dl*twopi*dkos / (dxos**2*nos*nvol)
#endif
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: cor_fac												! factored power at a particular k mode		
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv																		! vector of Gaussian random variables
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: Fk_cor											! vector to hold DFT of fld1, fldp1, fld2, fldp2, ...

		print*, "kcut: ", kcut
		
		! Initialize random number generator
		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
    !seed = 27*(/ (i-1, i=1,nseed) /)
		seed = seed_in*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)

		! Cholesky factorization of power
		do k=1,nz; do j=1,ny; do i=1,nnx
			if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
			if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			ii = i-1
			do m = 1,2*nfld; do n = 1,2*nfld
				cor_fac(m,n) = cor_in(m,n,LATIND)
			enddo; enddo
			!call dpotrf('L',2*nfld,cor_fac,2*nfld,l)
			call zpotrf('L',2*nfld,cor_fac,2*nfld,l)
			if (l /= 0) then
				print*, "Factorization warning: l = ", l
			endif
			! Initialize GRV vector
			call random_number(a); call random_number(p)	! a and p are dimension 2*nfld
			if (i==1 .and. j==1 .and. k==1) then
				grv(:) = (0._dl,0._dl)
			elseif (i==1 .or. i==nnx) then
				grv(:) = (1._dl,0._dl)*sign(sqrt(-2._dl*log(a)),p-0.5_dl)!/sqrt(2._dl*dble(nfld))	! Normalize grv accounting for number of fields. Is this one right?
			else
				!grv(:) =	sqrt(-1._dl*log(a)) * exp(w*p)			! Removed factor of sqrt(2) from Box Mueller transform to normalize for complex grv. That was wrong needs normalization based on nfld
				grv(:) =	sqrt(-1._dl*log(a)) * exp(w*p)!/sqrt(2._dl*dble(nfld))	! Normalize grv accounting for number of fields
			endif
			!print*, "grv norm= ", sqrt((grv*conjg(grv)))
			do m = 1,2*nfld; do n = 1,2*nfld	! can do this with a where statement
				if (n<m) then
					cor_fac(n,m) = (0._dl,0._dl)
				endif
			enddo; enddo
			call ztrmv('L','N','N', 2*nfld, cor_fac, 2*nfld, grv,  1)	! multiply random vector by transposed lower triangular factor of matrix
			!Fk_cor(:,LATIND) = grv(:)/nvol*exp(-0.5*dble(ii**2+jj*2+kk**2)/(kcut**2))	! filtering at Nyquist
			Fk_cor(:,LATIND) = grv(:)/nvol!*exp(-0.5*dble(ii**2+jj**2+kk**2)/(kcut**2))*exp(-0.5*dble(ii**2+jj**2+kk**2)*dk**2/(kcut_h**2))	! fitering at horizon, modified Nyq filtering
			Fk_cor(:,LATIND) = grv(:)/nvol*exp(-0.5*dble(ii**2+jj**2+kk**2)*dk**2/(kcut_h**2))	! fitering at horizon, modified Nyq filtering
			! Apply tophat filter to remove Nyquist frequency
			if (dble(ii**2+jj**2+kk**2) - kcut**2>0) then
				Fk_cor(:,LATIND) = grv(:)*0.0_dl
			endif
		enddo; enddo; enddo
		
		! Peform FFT to calculate fluctuations in real space
		! Add fluctuations to homogeneous fields
		! Check if /nvol is required, nvol is used when initializing Fk_cor
		do m=1, nfld
			Fk(:,:,:) = Fk_cor(2*m-1,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)!!!!! fft to real space
			fld(m,IRANGE) = fld0(m) + laplace(IRANGE)
#ifdef GHOST
			ghst(m,IRANGE) = laplace(IRANGE)
#endif
			Fk(:,:,:) = Fk_cor(2*m,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)
			fldp(m,IRANGE) = dfld0(m) + laplace(IRANGE)
#ifdef GHOST
			ghstp(m,IRANGE) = laplace(IRANGE)
#endif
		enddo

		yscl = 1._dl
    call calc_metric()
	end subroutine init_fields_cor2

	! subroutine to initialize fields given a power spectrum
	! n.b. cor_in is provided as
	subroutine init_fields_cor_rad(cor_in, seed_in)
		real(dl), dimension(2*nfld,2*nfld,nos_mode) :: cor_in						! power on some oversampled radial profile
		integer :: seed_in
		
		complex, parameter :: w = (0._dl, twopi)									! i*2*pi, used for random phase
		real(dl) :: a(2*nfld), p(2*nfld)													! amplitude, phase
    integer, allocatable :: seed(:)														! rng variable
    integer :: nseed																					! rng variable

		integer :: i, j, k, n, m, l																! lattice location indicies/ matrix indicies, check integer
		integer :: ii,jj,kk																				! mode numbers used for filtering
		real(dl) :: rad																						! radius squared of wavenumber
		real(dl), parameter :: kcut = KCUT_FAC*(min(nnx,nny,nnz)-1)!sqrt(2.*KCUT_FAC*(min(nnx,nny,nnz)-1))	!	frequency filtering at Nyquist
		real(dl), parameter :: kcut_h = KCUT_FAC_H*H0									! frequency filtering at horizon

		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: cor_fac												! factored power at a particular k mode		
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv																		! vector of Gaussian random variables
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: Fk_cor											! vector to hold DFT of fld1, fldp1, fld2, fldp2, ...

		! Initialize random number generator
		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
		seed = seed_in*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)

		! Loop over wavenumbers on lattice
		do k=1,nz; do j=1,ny; do i=1,nnx
			if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
			if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			ii = i-1
			! Interpolate wavenumber
			rad = sqrt(dble(ii**2 + jj**2 + kk**2))
			l = floor(rad*os_mode)
			do m = 1,2*nfld; do n = 1,2*nfld
				cor_fac(m,n) = (1._dl,0._dl)*(rad-dble(l)/dble(os_mode))*cor_in(m,n,l+1) + (1._dl,0._dl)*(1._dl - rad + dble(l)/dble(os_mode))*cor_in(m,n,l)	
			enddo; enddo
			!if (400<=ii**2+jj**2+kk**2) then
			!	print*, "ii, jj, kk, rad, l = ", ii, jj, kk, rad, l
			!	print*, cor_fac/(nvol**2)
			!	print*
			!endif
			! Choleski factorization
			call zpotrf('L',2*nfld,cor_fac,2*nfld,l)
			if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1)) then
				print*, "Factorization warning: l = ", l
				print*, "ii, jj, kk, rad, l = ", ii, jj, kk, rad, l
				print*, cor_fac
				print*
				!print*, "i, j, k = ", i ,j ,k
			endif
			! Initialize GRV vector
			call random_number(a); call random_number(p)	! a and p are dimension 2*nfld
			if (i==1 .and. j==1 .and. k==1) then
				grv(:) = (0._dl,0._dl)
			elseif (i==1 .or. i==nnx) then
				grv(:) = (1._dl,0._dl)*sign(sqrt(-1._dl*log(a)),p-0.5_dl)!(1._dl,0._dl)*sign(sqrt(-2._dl*log(a)),p-0.5_dl)
			else
				grv(:) =	sqrt(-1._dl*log(a)) * exp(w*p) !sqrt(-0.5_dl*log(a)) * exp(w*p)!
			endif
			!print*, "grv norm= ", sqrt((grv*conjg(grv)))
			do m = 1,2*nfld; do n = 1,2*nfld	! can do this with a where statement
				if (n<m) then
					cor_fac(n,m) = (0._dl,0._dl)
				endif
			enddo; enddo
			call ztrmv('L','N','N', 2*nfld, cor_fac, 2*nfld, grv,  1)	! multiply random vector by transposed lower triangular factor of matrix
			!Fk_cor(:,LATIND) = grv(:)/nvol*exp(-0.5*dble(ii**2+jj**2+kk**2)*dk**2/(kcut_h**2))	! fitering at horizon
			Fk_cor(:,LATIND) = grv(:)/nvol	! no filtering at horizon
			! Apply tophat filter to remove Nyquist frequency
			if (dble(ii**2+jj**2+kk**2) - kcut**2 > 0) then
				Fk_cor(:,LATIND) = grv(:)*0.0_dl
			endif
		enddo; enddo; enddo
			
		! Peform FFT to calculate fluctuations in real space
		! Add fluctuations to homogeneous fields
		! Check if /nvol is required, nvol is used when initializing Fk_cor
		do m=1, nfld
			Fk(:,:,:) = Fk_cor(2*m-1,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)!!!!! fft to real space
			fld(m,IRANGE) = fld0(m) + laplace(IRANGE)
			Fk(:,:,:) = Fk_cor(2*m,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)
			fldp(m,IRANGE) = dfld0(m) + laplace(IRANGE)
		enddo

		yscl = 1._dl
    call calc_metric()
	end subroutine init_fields_cor_rad

	! subroutine to initialize fields given a correlation matrix and zeta consistant with the fields up to linear order
	! to do: Giving factorization warning from 0 matrix
	! to do:
	subroutine init_fields_cor_z_rad(cor_in, z_amp_in, z_phase_in, seed_in)
		real(dl), dimension(2*nfld,2*nfld,nos_mode) :: cor_in						! correlation matrix on some oversampled radial profile
		real(dl), dimension(nos_mode) :: z_amp_in												! Amplitude of zeta modes on oversampled radial profile
		real(dl), dimension(nos_mode) :: z_phase_in											! relative phase of zeta compared to phi on radial profile
		integer :: seed_in																							! random number generator seed
		
		complex, parameter :: w = (0._dl, twopi)									! i*2*pi, used for random phase
		real(dl) :: a(2*nfld), p(2*nfld)													! amplitude, phase
    integer, allocatable :: seed(:)														! rng variable
    integer :: nseed																					! rng variable

		integer :: i, j, k, n, m, l																! lattice location indicies/ matrix indicies, check integer
		integer :: ii,jj,kk																				! mode numbers used for filtering
		real(dl) :: rad																						! radius squared of wavenumber
		real(dl), parameter :: kcut = KCUT_FAC*(min(nnx,nny,nnz)-1)!sqrt(2.*KCUT_FAC*(min(nnx,nny,nnz)-1))	!	frequency filtering at Nyquist
		real(dl), parameter :: kcut_h = KCUT_FAC_H*H0									! frequency filtering at horizon

		real(dl) :: z_amp
		real(dl) :: z_phase
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: cor_fac												! factored power at a particular k mode		
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv																		! vector of Gaussian random variables
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: Fk_cor											! vector to hold DFT of fld1, fldp1, fld2, fldp2, ...

		! Initialize random number generator
		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
		seed = seed_in*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)

		!print*, cor_in

		! Loop over wavenumbers on lattice
		do k=1,nz; do j=1,ny; do i=1,nnx
			if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
			if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			ii = i-1
			! Interpolate wavenumber
			rad = sqrt(dble(ii**2 + jj**2 + kk**2))
			l = floor(rad*os_mode)
			z_amp = (rad-dble(l)/dble(os_mode))*z_amp_in(l+1) + (1._dl - rad + dble(l)/dble(os_mode))*z_amp_in(l)
			z_phase = (rad-dble(l)/dble(os_mode))*z_phase_in(l+1) + (1._dl - rad + dble(l)/dble(os_mode))*z_phase_in(l)
			do m = 1,2*nfld; do n = 1,2*nfld
				cor_fac(m,n) = (1._dl,0._dl)*(rad-dble(l)/dble(os_mode))*cor_in(m,n,l+1) + (1._dl,0._dl)*(1._dl - rad + dble(l)/dble(os_mode))*cor_in(m,n,l)	
			enddo; enddo
			! Choleski factorization
			!print*, cor_fac
			call zpotrf('L',2*nfld,cor_fac,2*nfld,l)
			if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1)) then
				print*, "Factorization warning: l = ", l
				print*, "ii, jj, kk, rad, l = ", ii, jj, kk, rad, l
				print*, cor_fac
				print*
			endif
			! Initialize GRV vector
			call random_number(a); call random_number(p)	! a and p are dimension 2*nfld							
			if (i==1 .and. j==1 .and. k==1) then
				grv(:) = (0._dl,0._dl)
			elseif (i==1 .or. i==nnx) then
				grv(:) = (1._dl,0._dl)*sign(sqrt(-1._dl*log(a)),p-0.5_dl)!(1._dl,0._dl)*sign(sqrt(-2._dl*log(a)),p-0.5_dl)
			else
				grv(:) =	sqrt(-1._dl*log(a)) * exp(w*p) !sqrt(-0.5_dl*log(a)) * exp(w*p)!
			endif
			! Set Fourier transform of zeta
			!print*, "a,p = ", a, p
			!print*, "grv(:) = ", grv(:)
			!print*, "grv(1) = ", grv(1)
			!print*, "z_amp = ", z_amp
			! move z_amp = line to after modes are multipled by matrix
			!z_amp = z_amp * abs(grv(1))																															! find amplitude of zeta mode
			!print*, "z_amp = ", z_amp	
			!print*, "abs(grv(1)) = ", abs(grv(1))	
			!print*, "AIMAG(grv(1)) = ", AIMAG(grv(1))	
			!print*, "REAL(grv(1)) = ", REAL(grv(1))					
			!print*, "atan2(AIMAG(grv(1)), REAL(grv(1)))", atan2(AIMAG(grv(1)), REAL(grv(1)))
			!z_phase = z_phase - atan2(AIMAG(grv(1)), REAL(grv(1)))					! find phase of zeta mode
			!Fk(LATIND) = z_amp * exp((0._dl,-1._dl)*z_phase) / nvol
			do m = 1,2*nfld; do n = 1,2*nfld	! can do this with a where statement
				if (n<m) then
					cor_fac(n,m) = (0._dl,0._dl)
				endif
			enddo; enddo
			call ztrmv('L','N','N', 2*nfld, cor_fac, 2*nfld, grv,  1)	! multiply random vector by transposed lower triangular factor of matrix
			z_amp = z_amp * abs(grv(1))
			z_phase = z_phase - atan2(AIMAG(grv(1)), REAL(grv(1)))					! find phase of zeta mode
			Fk(LATIND) = z_amp * exp((0._dl,-1._dl)*z_phase) / nvol
			! Set Fourier transform of fields
			!Fk_cor(:,LATIND) = grv(:)/nvol*exp(-0.5*dble(ii**2+jj**2+kk**2)*dk**2/(kcut_h**2))	! fitering at horizon
			Fk_cor(:,LATIND) = grv(:)/nvol	! no filtering at horizon
			! Apply tophat filter to remove Nyquist frequency
			if (dble(ii**2+jj**2+kk**2) - kcut**2 > 0) then
				Fk_cor(:,LATIND) = grv(:)*0.0_dl
				Fk(LATIND) = (0.0_dl, 0.0_dl)
			endif
		enddo; enddo; enddo
			
		! Peform FFT to calculate fluctuations in real space
		! Add fluctuations to homogeneous fields
		call fftw_execute_dft_c2r(planb, Fk, laplace)
		zeta_lat(IRANGE) = laplace(IRANGE)
		!print*, zeta_lat
		do m=1, nfld
			Fk(:,:,:) = Fk_cor(2*m-1,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)!!!!! fft to real space
			fld(m,IRANGE) = fld0(m) + laplace(IRANGE)
			Fk(:,:,:) = Fk_cor(2*m,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)
			fldp(m,IRANGE) = dfld0(m) + laplace(IRANGE)
		enddo

		yscl = 1._dl
    call calc_metric()
	end subroutine init_fields_cor_z_rad

	! Alternative field initialization routine to use for convergence tests. Initializes fields using the Fourier modes of a lattice with double the grid spacing.
	! to do: order in which modes are initialized must match that used in init_fields_cor2
	! to do: Nyquist frequency (and above) must be removed, so apply a tophat filter
	! to do: the filter must remove exactly those modes which are not found on the smaller lattice
	! to do: those modes which are found on the smaller lattice must be initialized in exactly the same way
	! to do: check if the lattices match at the even, or odd points
	subroutine init_fields_cor_conv_test(cor_in, seed_in)
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld,nnx,ny,nz) :: cor_in			! Input power/cross spectrum for each k on the lattice
		integer :: seed_in

    complex, parameter :: w = (0._dl, twopi)									! i*2*pi, used for random phase
		real(dl) :: a(2*nfld), p(2*nfld)													! amplitude, phase
    integer, allocatable :: seed(:)														! rng variable
    integer :: nseed																					! rng variable

    integer :: i, j, k, n, m, l																! lattice location indicies/ matrix indicies, check integer
		integer :: ii,jj,kk																				! mode numbers used for filtering
		real(dl), parameter :: kcut = KCUT_FAC*(min(nx/4,ny/4,nz/4))!sqrt(2.*KCUT_FAC*(min(nnx,nny,nnz)-1))	!	frequency filtering at Nyquist
		real(dl), parameter :: kcut_h = KCUT_FAC_H*H0									! frequency filtering at horizon

		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: cor_fac												! factored power at a particular k mode		
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv																		! vector of Gaussian random variables
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: Fk_cor											! vector to hold DFT of fld1, fldp1, fld2, fldp2, ...
		
		integer, dimension(nx/2) :: k_conv

		! Ordering of Fourier modes for lattice with double grid spacing
		k_conv = (/ (k, k=1,nx/2) /)
		where (k_conv > nx/4+1)
			k_conv = k_conv+nx/2
		elsewhere
			k_conv = k_conv
		end where 

		print*, "kcut: ", kcut
		
		! Initialize random number generator
		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
    !seed = 27*(/ (i-1, i=1,nseed) /)
		seed = seed_in*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)

		! Put all modes to zero
		Fk_cor(:,:,:,:) = 0._dl

		! Cholesky factorization of power
		do k=1,nz/2; do j=1,ny/2; do i=1,nx/4+1
			if (k>nz/4+1) then; kk = k_conv(k)-nz-1; else; kk=k-1; endif
			if (j>ny/4+1) then; jj = k_conv(j)-ny-1; else; jj=j-1; endif
			ii = i-1
			do m = 1,2*nfld; do n = 1,2*nfld
				cor_fac(m,n) = cor_in(m,n,i,k_conv(j),k_conv(k))
			enddo; enddo
			call zpotrf('L',2*nfld,cor_fac,2*nfld,l)	! Factorize matrix
			if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1)) then
				print*, "Factorization warning: l = ", l
			endif
			! Initialize GRV vector
			call random_number(a); call random_number(p)	! a and p are dimension 2*nfld
			if (i==1 .and. j==1 .and. k==1) then
				grv(:) = (0._dl,0._dl)
			elseif (i==1 .or. i==nx/4+1) then
				grv(:) = (1._dl,0._dl)*sign(sqrt(-2._dl*log(a)),p-0.5_dl)	! Box Mueller transform
			else
				grv(:) =	sqrt(-1._dl*log(a)) * exp(w*p)
			endif
			!print*, "grv norm= ", sqrt((grv*conjg(grv)))
			do m = 1,2*nfld; do n = 1,2*nfld	! can do this with a where statement
				if (n<m) then
					cor_fac(n,m) = (0._dl,0._dl)
				endif
			enddo; enddo
			call ztrmv('L','N','N', 2*nfld, cor_fac, 2*nfld, grv,  1)	! multiply random vector by transposed lower triangular factor of matrix
			Fk_cor(:,i,k_conv(j),k_conv(k)) = grv(:)/nvol!*exp(-0.5*dble(ii**2+jj**2+kk**2)/(kcut**2))*exp(-0.5*dble(ii**2+jj**2+kk**2)*dk**2/(kcut_h**2))	! fitering at horizon, modified Nyq filtering
			! Apply tophat filter to remove Nyquist frequency
			if (dble(ii**2+jj**2+kk**2) - kcut**2>0) then
				Fk_cor(:,LATIND) = grv(:)*0.0_dl
			endif
		enddo; enddo; enddo
		
		! Peform FFT to calculate fluctuations in real space
		! Add fluctuations to homogeneous fields
		do m=1, nfld
			Fk(:,:,:) = Fk_cor(2*m-1,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)!!!!! fft to real space
			fld(m,IRANGE) = fld0(m) + laplace(IRANGE)
#ifdef GHOST
			ghst(m,IRANGE) = laplace(IRANGE)
#endif
			Fk(:,:,:) = Fk_cor(2*m,:,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)
			fldp(m,IRANGE) = dfld0(m) + laplace(IRANGE)
#ifdef GHOST
			ghstp(m,IRANGE) = laplace(IRANGE)
#endif
		enddo

		yscl = 1._dl
    call calc_metric()
	end subroutine init_fields_cor_conv_test

	! subroutine to initialize fields independently for testing purposes
	! to do: the reality condition inforces that some of the entries to Fk_cor must be real, do this.
	! to do: since this is the power spectrum of the fluctuations the zero k-mode must be set to zero
	subroutine init_fields_ind(seed_in)
		integer :: seed_in

    complex, parameter :: w = (0._dl, twopi)									! i*2*pi, used for random phase
		real(dl) :: a(2*nfld), p(2*nfld)													! amplitude, phase
    integer, allocatable :: seed(:)														! rng variable
    integer :: nseed																					! rng variable

    integer :: i, j, k, n, m, l																! lattice location indicies/ matrix indicies, check integer

		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: cor_fac												! factored power at a particular k mode		
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv																		! vector of Gaussian random variables
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: Fk_cor											! vector to hold DFT of fld1, fldp1, fld2, fldp2, ...
		
		! Initialize random number generator
		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
    !seed = 27*(/ (i-1, i=1,nseed) /)
		seed = seed_in*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)

		! Initialize power spectrum (independent fields, uniform power)
		do k=1,nz; do j=1,ny; do i=1,nnx
			! Initialize GRV vector
			call random_number(a); call random_number(p)	! a and p are dimension 2*nfld
			!grv(:) =	sqrt(-1._dl*log(a)) * exp(w*p)			! Removed factor of sqrt(2) from Box Mueller transform to normalize for complex grv
			if (k==1 .and. j==1 .and. i==1) then
				grv(:) = (0._dl, 0._dl)
			elseif (i==1 .or. i==nnx) then
				grv(:) = (1._dl, 0._dl)
			else
				grv(:) =	(0._dl, 1._dl)!* exp(w*p)																! Amplitudes set to 1 for testing, phase set to 0
			endif
			Fk_cor(:,LATIND) = grv(:)/nvol			
		enddo; enddo; enddo
			!Fk_cor(:,:,:,:) = (1._dl,0._dl)/nvol

		
		! Peform FFT to calculate fluctuations in real space
		! Add fluctuations to homogeneous fields
		! Check if /nvol is required, nvol is used when initializing Fk_cor
		do m=1, nfld
			Fk(:,:,:) = Fk_cor(2*m-1,:,:,:)
			!Fk(:,:,:) = (1._dl, 1._dl)/nvol
			!print*,"before:",Fk(3,:,:)
			call fftw_execute_dft_c2r(planb, Fk, laplace)!!!!! fft to real space
			fld(m,IRANGE) = fld0(m) + laplace(IRANGE)
			call fftw_execute_dft_r2c(planf, laplace, Fk)
			!print*,"after:", Fk(3,:,:)
			Fk(:,:,:) = Fk_cor(2*m,:,:,:)
			!Fk(:,:,:) = (0._dl, 1._dl)/nvol
			call fftw_execute_dft_c2r(planb, Fk, laplace)
			fldp(m,IRANGE) = dfld0(m) + laplace(IRANGE)
		enddo
		
		yscl = 1.
    call calc_metric()
	end subroutine init_fields_ind

  end program lattice

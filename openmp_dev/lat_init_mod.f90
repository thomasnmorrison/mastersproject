! lat_init_mod.f90

! Module for subroutines relating to initializing the lattice.

! to do: write subroutine that initializes fields from a supplied correlation
!        matrix and self consitantly initializes metric perturbations from
!        the same random variables
! to do: write a subroutine that will initialize zeta from the metric
!        perturbations
! to do: write a subroutine to call the evolution of the background, field
!        perturbations, and metric pertubations
! to do: decouple this module from hamiltonian so it can be tested and used independently
! to do: need to initialize yscl separately, since it is in the hamiltonian module which is not used

#define LATINITOPT 1
! LATINITOPT 0 for homogeneous initialization
! LATINITOPT 1 for spectral initialization
! LATINITOPT 2 for convolution initialization
#define ZETAINIT 0
! ZETAINIT 0 to initialize zeta to 0
! ZETAINIT 1 to initialize zeta from constraint equation
!          assuming longitudinal gauge

module lat_init_mod
#include "macros.h"
  use grv_mod
  use fftw3
  use omp_lib
  use params
  use potential_mod
  use analysis
  !use bg_cosmo_mod
  !use pert_cosmo_mod
  !use field_pert_mod
  !use metric_pert_mod
  !use Hamiltonian

implicit none

  !complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_fp_mat  

contains

  ! Subroutine that wraps various options for lattice initialization
  ! to do: update inputs
  !subroutine lat_init(corr_in, k_cut, seed_in)
  subroutine lat_init(corr_in, kos, k_cut, seed_in, f0_in, df0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
    !real(dl) :: corr_in(:,:,:)!corr_in(2*nfld,2*nfld,fp_nos)
    !real(dl), intent(in) :: k_cut
    !integer, intent(in) :: seed_in
    real(dl) :: corr_in(:,:,:)      ! input correlation !corr_in(2*nfld,2*nfld,cp_nos)
    integer, intent(in) :: kos     ! over sampling of k in corr_in
    real(dl), intent(in) :: k_cut   ! value of k cutoff
    integer, intent(in) :: seed_in  ! seed for rng
    real(dl), intent(in) :: f0_in(:)  ! mean fields inputs
    real(dl), intent(in) :: df0_in(:) ! mean momenta inputs
    real(dl) :: f_init(:,:,:,:)     ! array of fields to be initialized
    real(dl) :: df_init(:,:,:,:)    ! array of momenta to be initialized
    real(dl) :: zeta_init(:,:,:)    ! zeta to be initialized
    real(C_DOUBLE), pointer :: f(:,:,:)								! pointer field for lattice
		complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)	  ! pointer for fft
    type(C_PTR) :: planb
    real(dl), intent(in) :: norm_in		

#if (LATINITOPT == 0)
    call lat_init_homo()
#elif (LATINITOPT == 1)
    call lat_init_spec(corr_in, kos, k_cut, seed_in, f0_in, df0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
#elif (LATINITOPT == 2)
    call lat_init_convolution(corr_in, k_cut, seed_in)
#endif

  end subroutine lat_init

  ! Subroutine to initialize homogeneous fields on the lattice from input. 
  subroutine lat_init_homo(f0_in, df0_in,f_init, df_init, zeta_init)
    real(dl), intent(in) :: f0_in(:)
    real(dl), intent(in) :: df0_in(:)
    real(dl) :: f_init(:,:,:,:)
    real(dl) :: df_init(:,:,:,:)
    real(dl) :: zeta_init(:,:,:)

    integer :: n

    do n=1,nfld
      f_init(n,:,:,:) = f0_in(n)
      df_init(n,:,:,:) = df0_in(n)
    enddo
    zeta_init(:,:,:) = 0._dl
  end subroutine lat_init_homo

  ! to do: should probably delete corr_Fk after calling
  ! to do: make sure my index ordering is consistant here and in pert_cosmo_test.f90
  ! n.b. zeta mode realization assumes yscl=1
  ! n.b. wrap fields after initializing
  ! to do: Figure out what is happening with passing arrays of bounds starting at 0 into this
  !        subroutine. may be indexing from 1 internally.
  ! to do: off by a factor of sqrt(2)
  ! to do: check storing in fftw
  subroutine lat_init_spec(corr_in, kos, k_cut, seed_in, f0_in, df0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
    real(dl) :: corr_in(:,:,:)      ! input correlation !corr_in(2*nfld,2*nfld,cp_nos)
    integer, intent(in) :: kos     ! over sampling of k in corr_in
    real(dl), intent(in) :: k_cut   ! value of k cutoff
    integer, intent(in) :: seed_in  ! seed for rng
    real(dl), intent(in) :: f0_in(:)  ! mean fields inputs
    real(dl), intent(in) :: df0_in(:) ! mean momenta inputs
    real(dl) :: f_init(:,:,:,:)     ! array of fields to be initialized
    real(dl) :: df_init(:,:,:,:)    ! array of momenta to be initialized
    real(dl) :: zeta_init(:,:,:)    ! zeta to be initialized
    real(C_DOUBLE), pointer :: f(:,:,:)								! pointer field for lattice
		complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)	  ! pointer for fft
    type(C_PTR) :: planb												      ! fft c2r plan
    real(dl), intent(in) :: norm_in

    integer :: i, j, k
    integer :: ii, jj, kk, l
    real(dl) :: rad
    integer :: n, m

    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: corr_fac
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: corr_Fk

    complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv
   
    ! initialize rng
    call init_rng(seed_in)    

    ! loop over wavenumbers
    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
      do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif 
        do i=1,nnx; ii=i-1
          rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))
          if (rad*dk .GE. k_cut) then
            corr_Fk(:,LATIND) = 0._dl*grv(:) ! cut field fluctuations
          else
            ! initialize corr_fac
            l = floor(rad*kos)
            if (l>0) then
              corr_fac(:,:) = (1._dl,0._dl)*((1._dl + l - rad*kos)*corr_in(:,:,l) & !((1._dl-l+rad*kos)*corr_in(:,:,l) &
                                             + (rad*kos - l)*corr_in(:,:,l+1))
            else
              corr_fac(:,:) = (1._dl,0._dl)*((rad-1._dl)*corr_in(:,:,2) &
                                             + (2._dl-rad)*corr_in(:,:,1))
            endif
            ! Choleski factorization of corr_fac
            call zpotrf('L',2*nfld,corr_fac,2*nfld,l)
            ! print warning if factorization failed
            if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1)) then
              print*, "Factorization Warning: l = ", l
              print*, "ii, jj, kk = ", ii, jj, kk
              print*, "corr_fac = ", corr_fac
            endif
            ! call Gaussian random variables
            if (i==1 .and. j==1 .and. k==1) then
              grv(:) = (0._dl, 0._dl)
            !elseif (i==1 .or. i==nx) then
            !  grv(:) = (1._dl, 0._dl)*get_grv_real(2*nfld)
            else
              grv(:) = get_grv_complex(2*nfld)
            endif
            ! enforce triangular form on correlation matrix
            do m=1,2*nfld; do n=1,2*nfld
              if (n<m) then
                corr_fac(n,m) = (0._dl, 0._dl)
              endif
            enddo; enddo
            ! realize field modes
            !print*, 'corr_fac(2,2): ', corr_fac(2,2)
            !print*, 'grv(:): ', grv(:)
            call ztrmv('L','N','N', 2*nfld, corr_fac, 2*nfld, grv, 1)
            corr_Fk(:,LATIND) = grv(:)/nvol*norm_in
          endif
          ! cut modes above kcut
          !if (rad .GE. k_cut) then
          !  corr_Fk(:,LATIND) = 0._dl*grv(:) ! cut field fluctuations
          !endif
#if (ZETAINIT == 1)
          ! Realize modes for zeta
          Fk(LATIND) = (0._dl,0._dl)
          do n=1,nfld
            Fk(LATIND) = Fk(LATIND) + ((bg_dv(fld0,n) + 3._dl*H0*dfld(n))/(sum(dfld0**2) - 2._dl*rad) &
                                       + (bg_dv(fld0,n))/(3._dl*sum(dfld0**2)))*corr_Fk(2*n-1,LATIND)
            Fk(LATIND) = Fk(LATIND) + ((dfld0(n))/(sum(dfld0**2) - 2._dl*rad) &
                                       + (dfld0(n))/(3._dl*sum(dfld0**2)))*corr_Fk(2*n,LATIND)
          enddo
#endif
        enddo
      enddo
    enddo
#if (ZETAINIT == 0)
    zeta_init(:,:,:) = 0._dl
#elif (ZETAINIT == 1)
    ! Invert FFT to initialize zeta
    call fftw_execute_dft_c2r(planb, Fk, f)
    zeta_init(:,:,:) = f(:,:,:)
#endif

    ! initialize homogeneous fields
    call lat_init_homo(f0_in,df0_in,f_init, df_init, zeta_init)
    ! invert FFT, add fluctuations to homogeneous fields
    do m=1,nfld
      Fk(:,:,:) = corr_Fk(2*m-1,:,:,:)
      call fftw_execute_dft_c2r(planb, Fk, f)
      f_init(m,IRANGE) = f_init(m,IRANGE) + f(IRANGE)
      Fk(:,:,:) = corr_Fk(2*m,:,:,:)
      !print*, 'Fk: ', Fk
      call fftw_execute_dft_c2r(planb, Fk, f)
      df_init(m,IRANGE) = df_init(m,IRANGE) + f(IRANGE)
    enddo
    ! set scale factor
    !yscl = 1._dl
  end subroutine lat_init_spec


  ! Subroutine to initialize the fields on the lattice to match Fourier space
  ! correlations using spectral initialization. Initializes zeta from metric 
  ! perturbations initialized  consitantly with the field fluctuations.
  ! to do: initialize metric perturbations in line
  ! to do: multiply mp_fp_interp and mp_fms_interp
  ! to do: multiply product by grv(1:2)
  ! to do: initialize zeta fft inline do this by filling Fk
  ! to do: check if offset is required for wavenumber from corr_in
  ! to do: matrix multiplication alpha check data type
!  subroutine lat_init_spec(corr_in, k_cut, seed_in)
!    real(dl) :: corr_in(2*nfld,2*nfld,fp_nos)
!    real(dl), intent(in) :: k_cut
!    integer, intent(in) :: seed_in
!
!    integer :: i, j, k
!    integer :: ii, jj, kk, l
!    real(dl) :: rad
!    integer :: n, m
!
!    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: corr_fac
!    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: corr_Fk
!
!    complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_fp_interp
!    complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_fms_interp
!    complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_init
!
!    complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv
!    complex(C_DOUBLE_COMPLEX), dimension(2) :: mp
!   
!    ! initialize rng
!    call init_rng(seed_in)    
!
!    ! loop over wavenumbers
!    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else kk=k-1; endif
!      do j=1,ny; if (j>nny) then; jj = ny+1-j; else jj=j-1; endif 
!        do i=1,nnx; ii=i-1
!          rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))
!          ! initialize corr_fac
!          l = floor(rad*fp_kos)
!          corr_fac(:,:)=(1._dl,0._dl)*((1._dl-l+rad*fp_kos)*corr_in(:,:,l)&
!                                       + (rad*fp_kos-l)*corr_in(:,:,l+1))
!          !do m=1,2*nfld; do n=1,2*nfld
!          !  corr_fac(m,n) = (1._dl,0._dl)*()
!          !enddo; enddo
!          l = floor(rad*mp_kos)
!          mp_fp_interp(:,:)=(1._dl,0._dl)*((1._dl-l+rad*mp_kos)*mp_fp_mat(:,:,l)&
!                                             + (rad*mp_kos-l)*mp_fp_mat(:,:,l+1))
!          mp_fms_interp(:,:)=(1._dl,0._dl)*((1._dl-l+rad*mp_kos)*mp_fms(:,:,l) &
!                                             + (rad*mp_kos-l)*mp_fms(:,:,l+1))
!          ! Multiply mp_fms by mp_fp
!          call zgemm('N','N',2,2,2,(1._dl,0._dl),mp_fms_interp,2, &
!                     mp_fp_interp,2,(0._dl,0._dl),mp_init,2)
!          ! Choleski factorization of corr_fac
!          call zpotrf('L',2*nfld,corr_fac,2*nfld,l)
!          ! print warning if factorization failed
!          if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1))
!            print*, "Factorization Warning: l = ", l
!            print*, "ii, jj, kk = ", ii, jj, kk
!            print*, "corr_fac = ", corr_fac
!          endif
!          ! call Gaussian random variables
!          if (i==1 .and. j==1 .and. k==1) then
!            grv(:) = (0._dl, 0._dl)
!          elseif (i==1 .or. i==nx) then
!            grv(:) = (1._dl, 0._dl)*get_grv_real(2*nfld)
!          else
!            grv(:) = get_grv_complex(2*nfld)
!          endif
!          ! realize metric perturpation modes
!          call zgemv('N',2,2,(1._dl,0._dl),mp_init,2,grv(1:2),1, &
!                     (0._dl,0._dl),mp,1)
!          ! enforce triangular form on correlation matrix
!          do m=1,2*nfld; do n=1,2*nfld
!            if (n<m) then
!              corr_fac(n,m) = (0._dl, 0._dl)
!            endif
!          enddo
!          ! realize field modes
!          call ztrmv('L','N','N', 2*nfld, corr_fac, n*nfld, grv, 1)
!          corr_Fk(:,LATIND) = grv(:)/nvol
!          ! realize zeta modes
!          call set_bg_cosmo_ic(a_bg_0,hub_bg_0,f_bg_0,df_bg_0)
!          Fk(LATIND) = (2._dl*(mp(2)/hub_bg_0 + mp(1)) &
!                       /(3._dl*(1._dl+get_bg_w_conf())) + mp(1))/nvol
!          ! cut modes above kcut
!          if (rad .GE. k_cut) then
!            corr_Fk(:,LATIND) = 0._dl*grv(:) ! cut field fluctuations
!            Fk(LATIND) = (0._dl,0._dl) ! cut zeta fluctuations
!          endif
!        enddo
!      enddo
!    enddo

!    ! invert FFT to initialize zeta
!    call fftw_execute_dtf_c2r(planb, Fk, laplace)
!    zeta_lat(:,:,:) = laplace(:,:,:)
!
!    ! initialize homogeneous fields
!    call lat_init_homo()
!    ! invert FFT, add fluctuations to homogeneous fields
!    do m=1,nfld
!      Fk(:,:,:) = corr_Fk(2*m-1,:,:,:)
!      call fftw_execute_dft_c2r(planb, Fk, laplace)
!      fld(m,:,:,:) = fld(m,:,:,:) + laplace(:,:,:)
!      Fk(:,:,:) = corr_Fk(2*m,:,:,:)
!      call fftw_execute_dft_c2r(planb, Fk, laplace)
!      fldp(m,:,:,:) = fldp(m,:,:,:) + laplace(:,:,:)
!    enddo
!    ! set scale factor
!    yscl = 1._dl
!  end subroutine lat_init_spec


  ! Subroutine to initialize the fields on the lattice to match real space
  ! correlations using convolution initialization. Initializes zeta from metric
  ! perturbations initialized consitantly with the field fluctuations.
  ! to do: figure out how to 
  subroutine lat_init_conv()

  end subroutine lat_init_conv







end module lat_init_mod

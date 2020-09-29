! pert_cosmo_mod.f90

! Module for computing 2pt correlations of cosmological perturbations for a
! scalar field. Calculation done in conformal time and in longitudinal gauge.

! Calculation is done at linear level for an oversampled radial profile of
! modes.

! n.b. careful of field and momentum ordering in cp_corr
! n.b. 

! to do: apply cp_norm
! to do: clean up the ic setting situation
! to do: write initialization routine for next order in adiabatic approximation
! to do: init_cp_corr for cosmic time equations 


#define MODE_MP 0
! MODE_MP 0 for no metric perturbations in evolution
! MODE_PM 1 for linear metric perturbations in evolution
#define OMEGOPT 3
! OMEGOPT 1 for omega^2 = k^2 + a^2*m^2
! OMEGOPT 2 for omega^2 = k^2 + a^2*m^2 - 2(a'/a)^2
! OMEGOPT 3 for omega^2 = k^2 + a^2*m^2 - a''/a
#define INITOPT -1
! INITOPT -1 to call whatever the test initialization is
! INITOPT 1 to enforce ICs with (aL)'=0
! INITOPT 2 to set ICs with linear metric perturbations
! INITOPT 3 to set on the massless de Sitter solution
! INITOPT 4 to set with an adiabatic approximation
#define CPVAROPT 1
! CPVAROPT 0 to use the phi_A = L_AB chi_B set of variables
! CPVAROPT 1 to use the (a phi_A) = L_AB chi_B set of variables
! CPVAROPT 2 to use sqrt(2k)(a phi_A) = L_AB chi_B, k**2 delta_AB -2(a'/a)**2*cp_c_dxdx = <chi_A',chi_B'> set of variables
! CPVAROPT 3 to use sqrt(2k)(a phi_A) = L_AB chi_B, k**2/a**2 delta_AB + u_AB = <chi_A.,chi_B.> in cosmic time

module pert_cosmo_mod

  use params
  use potential_mod
  use gl_integrator_mod
  use bg_cosmo_mod

implicit none

  integer, parameter :: cp_sys_dim = bg_sys_dim+nfld*(2*nfld+1)! dimension of cosmo pert system
  integer, parameter :: cp_kos = 2**3 ! oversampling of k space radial profile 
  integer, parameter :: cp_nos = nn*cp_kos!max(nx,ny,nz)*cp_kos ! number of oversampled modes
  real(dl), parameter :: cp_dkos = dk/cp_kos!dk*cp_kos
  
  real(dl), parameter :: cp_norm = sqrt(nvol/dx**3/mpl**2) ! normalization for modes
  
  ! Derived quantities
  !real(dl), dimension(2,2) :: cp_constraint  ! mp-fp constraint matrix
  real(dl), dimension(2,2*nfld) :: cp_constraint  ! mp-fp constraint matrix
  real(dl), dimension(nfld,nfld) :: cp_m2    ! effective mass matrix
  real(dl), dimension(nfld,nfld) :: cp_a     ! antisym A matrix
  real(dl), dimension(nfld,nfld) :: cp_b     ! B matrix
  real(dl), dimension(nfld,nfld) :: cp_c     ! C matrix
  real(dl), dimension(nfld,nfld) :: cp_temp, cp_temp2  ! temp matrix
  ! Matrix quantities
  real(dl), dimension(nfld,nfld) :: cp_l_m
  real(dl), dimension(nfld,nfld) :: cp_dl_m
  real(dl), dimension(nfld,nfld) :: cp_c_xdx_m
  real(dl), dimension(nfld,nfld) :: cp_c_dxdx_m
  ! Evolved variables
  real(dl), target, dimension(nfld*(nfld+1)/2) :: cp_l
  real(dl), target, dimension(nfld*(nfld+1)/2) :: cp_dl
  real(dl), target, dimension((nfld-1)*nfld/2) :: cp_c_xdx
  real(dl), target, dimension(nfld*(nfld+1)/2) :: cp_c_dxdx

  !integer, target, dimension(cp_sys_dim) :: cp_key ! I think this is now extraneous so shouldn't need
  real(dl) :: cp_k2  ! square wavenumber

  real(dl) :: cp_sys_temp(cp_sys_dim)
  real(dl), allocatable :: cp_sys_g(:,:,:) ! GL integrator g
  type(ode_sys), dimension(cp_sys_dim) :: cp_ode

  real(dl), dimension(2*nfld,2*nfld,cp_nos) :: cp_corr ! 2pt corr mode by mode 
  real(dl) :: cp_corr_muksas          ! 2pt of MS variable

contains

  ! Subroutine to initialize pointers for the cp system
  ! to do: point to ode_param
  subroutine init_cp(da, dh, df, ddf, dl, ddl, dxdx, ddxdx)
    procedure(d_template) :: da, dh, df, ddf, dl, ddl, dxdx, ddxdx

    integer :: i,j,k,n1,n2

    n1 = nfld*(nfld+1)/2; n2 = (nfld-1)*nfld/2
    !print*, 'n1, n2 = ', n1, n2

    if(.NOT. allocated(cp_sys_g)) then
      call gl_init_g(cp_sys_g, cp_sys_dim)
    endif

    ! set pointers to derivative functions
    cp_ode(1)%f => da
    cp_ode(2)%f => dh
    do i=1,nfld
      cp_ode(2+i)%f => df
      cp_ode(2+nfld+i)%f => ddf
    enddo
    do k=1,n1
      cp_ode(bg_sys_dim+k)%f => dl
      cp_ode(bg_sys_dim+n1+k)%f => ddl
      cp_ode(bg_sys_dim+2*n1+n2+k)%f => ddxdx
    enddo
    do k=1,n2
      cp_ode(bg_sys_dim+2*n1+k)%f => dxdx
    enddo

    ! set pointers to system variables
    cp_ode(1)%y => a_bg
    cp_ode(2)%y => hub_bg
    do i=1,nfld
      cp_ode(2+i)%y => f_bg(i)
      cp_ode(2+nfld+i)%y => df_bg(i)
    enddo
    ! set pointer for matrix indices
    do k=1,n1
      cp_ode(bg_sys_dim+k)%y => cp_l(k)
      cp_ode(bg_sys_dim+n1+k)%y => cp_dl(k)
      cp_ode(bg_sys_dim+2*n1+n2+k)%y => cp_c_dxdx(k)
    enddo
    do k=1,n2
      cp_ode(bg_sys_dim+2*n1+k)%y => cp_c_xdx(k)
    enddo

    ! set keys
    do i=2+1,bg_sys_dim
      if(.NOT. allocated(cp_ode(i)%key)) then
        allocate(cp_ode(i)%key(1))
      endif
    enddo
    do i=1,nfld
      cp_ode(2+i)%key = (/i/)
      cp_ode(2+nfld+i)%key = (/i/)
    enddo
    ! set cp_l and cp_dl keys
    do i=bg_sys_dim+1,cp_sys_dim
      if(.NOT. allocated(cp_ode(i)%key)) then
        allocate(cp_ode(i)%key(2))
      endif
    enddo
    k = 1
    do i=1,nfld
      do j=1,i
        cp_ode(bg_sys_dim+k)%key = (/i,j/)
        cp_ode(bg_sys_dim+n1+k)%key = (/i,j/)
        cp_ode(bg_sys_dim+2*n1+n2+k)%key = (/i,j/)  ! check the +k here
        k = k+1 
      enddo
    enddo
    k = 1
    do i=2,nfld
      do j=1,i-1
        cp_ode(bg_sys_dim+2*n1+k)%key = (/i,j/)
        k = k+1
      enddo
    enddo

    ! set params
    do i=1, cp_sys_dim
      cp_ode(i)%ode_param => null()
    enddo
  end subroutine init_cp

  ! Subroutine to set cp_sys_temp from cp_ode
  subroutine set_cp_sys_temp()
    integer :: i
    do i=1, cp_sys_dim
      cp_sys_temp(i) = cp_ode(i)%y
    enddo
  end subroutine set_cp_sys_temp

  ! Subroutine to set ICs for the system of cp 2pt corr
  ! n.b. ICs for the bg part of the system need to be set separately
  ! n.b. eigenvectors are stored in columns
  ! to do: test this subroutine
  ! to do: check if LAPACK is storing things in the upper triangle of lower
  !        triangle matrices
#if (CPVAROPT == 0)
  subroutine set_cp_ic(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    
    integer :: i,j,k
    integer :: a_i, hub_i
    real(dl), dimension(nfld,nfld) :: u_m  ! matrix of eigenvectors
    real(dl), dimension(nfld) :: diag         ! array of eigenvalues
    real(dl), dimension((nfld+2)*nfld) :: work ! array for dsyev subroutine
#if (INITOPT == 1)
    a_i = 1; hub_i = 2
    cp_l = 0._dl; cp_dl = 0._dl; cp_c_xdx = 0._dl; cp_c_dxdx = 0._dl
    call get_cp_m2_conf(y, k2_in)
    !print*, 'cp_m2 = ', cp_m2(:,:)
    ! Diagonalize mass matrix
#if (OMEGOPT == 1)
    u_m = cp_m2
#elif (OMEGOPT == 2)
    u_m = cp_m2
    do i=1,nfld
      u_m(i,i) = u_m(i,i) - 2._dl*y(hub_i)**2
    enddo
#elif (OMEGOPT == 3)
    u_m = cp_m2
    do i=1,nfld
      u_m(i,i) = u_m(i,i) - (diff_bg_hub_conf(y) + y(hub_i)**2) ! k^2 + a^2V_{,ff} -a''/a
    enddo
#endif
    call dsyev('V','L',nfld,u_m,nfld,diag,work,(nfld+2)*2,i)
    !print*, 'eigenvalues cp_m2 = ', diag
    ! Compute amplitudes
    cp_m2 = 0._dl
    do i=1,nfld
      cp_m2(i,i) = 1._dl/(2._dl*y(a_i)**2*sqrt(diag(i)))
    enddo
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_m2,nfld,u_m,nfld, & 
               0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,u_m,nfld,cp_temp,nfld, & 
               0._dl,cp_l_m,nfld)
    call dpotrf('L',nfld,cp_l_m,nfld,i)
!#if (INITOPT == 1)
    cp_dl_m = -y(hub_i)*cp_l_m
!#endif
    ! Fill variables
    k = 1
    do i=1,nfld
      do j=1,i
        cp_l(k) = cp_l_m(i,j)
        cp_dl(k) = cp_dl_m(i,j)
        k = k+1
      enddo 
    enddo
    ! Compute correlations
    cp_m2 = 0._dl
    do i=1,nfld
      cp_m2(i,i) = sqrt(diag(i))/(2._dl*y(a_i)**2)
    enddo
    call dtrtri('L','N',nfld,cp_l_m,nfld,i)
    !print*, 'cp_l_m = ', cp_l_m
    !print*, 'u_m = ', u_m
    call dgemm('T','T',nfld,nfld,nfld,1._dl,u_m,nfld,cp_l_m,nfld, & 
               0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_m2,nfld,cp_temp,nfld, &
               0._dl,cp_temp2,nfld)
    call dgemm('T','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_temp2,nfld, &
               0._dl,cp_c_dxdx_m,nfld)
    !call dgemm('T','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp2,nfld, &
    !           0._dl,cp_c_dxdx_m,nfld)  ! I think this line has a bug should multiply cp_temp^T cp_temp2, otherwise miss the change of basis
    !print*, 'cp_c_dxdx_m = ', cp_c_dxdx_m
    ! Fill variables
    k = 0
    do i=1,nfld
      do j=1,i  
        k = k+1
        !cp_c_dxdx(k) = cp_c_dxdx_m(i,j) - y(hub_i)**2
        cp_c_dxdx(k) = cp_c_dxdx_m(i,j)
      enddo
    enddo
#endif

#if (INITOPT == -1)
    call set_cp_ic_test(y, k2_in)
#elif (INITOPT == 2)
    call set_cp_ic_mp(y,k2_in)
#elif(INITOPT == 3)
    call set_cp_ic_ds(y,k2_in)
#endif
  end subroutine set_cp_ic
#endif

  ! Subroutine to set ICs for the system of cp 2pt corr with the (a delta phi) equations
  ! n.b. ICs for the bg part of the system need to be set separately
  ! n.b. eigenvectors are stored in columns
  ! to do: extend to case of interacting fields
  ! to do: extend to the case of m^2 varying
  ! to do: extend to higher oreder in adiabatic approximation
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
  subroutine set_cp_ic(y, k2_in)
  !subroutine set_cp_ic_ascl(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    
    integer :: i,j,k
    integer :: a_i, h_i
    real(dl) :: ddh
    real(dl), dimension(nfld,nfld) :: u_m  ! matrix of eigenvectors
    real(dl), dimension(nfld) :: diag         ! array of eigenvalues
    real(dl), dimension((nfld+2)*nfld) :: work ! array for dsyev subroutine

    a_i = 1; h_i = 2
    cp_l = 0._dl; cp_dl = 0._dl; cp_c_xdx = 0._dl; cp_c_dxdx = 0._dl
    cp_l_m = 0._dl; cp_dl_m = 0._dl; cp_c_dxdx_m = 0._dl
    call get_cp_m2_conf(y, k2_in)!call get_cp_m2_ascl_conf(y, k2_in)
#if (INITOPT == -1)
    ddh = 0._dl
    do i=1,nfld
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
      ddh = ddh + y(2+nfld+i)*diff_bg_df_conf(y,key=(/i/))
#elif (CPVAROPT == 3)
      ddh = ddh - y(2+nfld+i)*diff_bg_df_cosm(y,key=(/i/)) ! d^2/dt^2(H)
#endif
    enddo
    do i=1,nfld
      cp_l_m(i,i) = 1._dl/sqrt(2._dl*sqrt(cp_m2(i,i)))
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
      cp_dl_m(i,i) = -(2._dl*y(h_i)*y(a_i)**2*bg_ddv(y(3:2+nfld), i, i) + ddh - 4._dl*y(h_i)*diff_bg_hub_conf(y))*cp_l_m(i,i)/4._dl
#elif (CPVAROPT == 3)
      cp_dl_m(i,i) = -y(a_i)**2*(2._dl*y(h_i)*bg_ddv(y(3:2+nfld), i, i) - 6._dl*diff_bg_hub_cosm(y)*y(h_i) &
                                  - 4._dl*y(h_i)**3 - ddh))*cp_l_m(i,i)/4._dl
#endif
      cp_c_dxdx_m(i,i) = cp_m2(i,i)
    enddo
#endif
#if (CPVAROPT == 2)
    cp_l_m = sqrt(2._dl*sqrt(k2_in))*cp_l_m
    cp_dl_m = sqrt(2._dl*sqrt(k2_in))*cp_dl_m
    do i=1,nfld
      cp_c_dxdx_m(i,i) = cp_c_dxdx_m(i,i) - k2_in
    enddo
    cp_c_dxdx_m = -cp_c_dxdx_m/(2._dl*y(h_i)**2)
#elif (CPVAROPT == 3)
    cp_l_m = sqrt(2._dl*sqrt(k2_in))*cp_l_m
    cp_dl_m = sqrt(2._dl*sqrt(k2_in))*cp_dl_m
    cp_c_dxdx_m = cp_c_dxdx_m/y(a_i)**2
    do i=1,nfld
      cp_c_dxdx_m(i,i) = cp_c_dxdx_m(i,i) - k2_in/y(a_i)**
    enddo
#endif
    k = 1
    do i=1,nfld
      do j=1,i
        cp_l(k) = cp_l_m(i,j)
        cp_dl(k) = cp_dl_m(i,j)
        cp_c_dxdx(k) = cp_c_dxdx_m(i,j)
        k = k+1
      enddo 
    enddo
  !end subroutine set_cp_ic_ascl
  end subroutine set_cp_ic

#elif (CPVAROPT == 3)
  subroutine set_cp_ic(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    
    integer :: i,j,k
    integer :: a_i, h_i
    real(dl) :: ddh
    real(dl), dimension(nfld,nfld) :: u_m  ! matrix of eigenvectors
    real(dl), dimension(nfld) :: diag         ! array of eigenvalues
    real(dl), dimension((nfld+2)*nfld) :: work ! array for dsyev subroutine

    a_i = 1; h_i = 2
    cp_l = 0._dl; cp_dl = 0._dl; cp_c_xdx = 0._dl; cp_c_dxdx = 0._dl
    cp_l_m = 0._dl; cp_dl_m = 0._dl; cp_c_dxdx_m = 0._dl
    call get_cp_m2_conf(y, k2_in)
#if (INITOPT == -1)
    ddh = 0._dl
    do i=1,nfld
      ddh = ddh + y(2+nfld+i)*diff_bg_df_cosm(y,key=(/i/))
    enddo
    do i=1,nfld
      cp_l_m(i,i) = 1._dl/sqrt(2._dl*sqrt(cp_m2(i,i)))
      cp_dl_m(i,i) = -y(a_i)**2*(2._dl*y(h_i)*bg_ddv(y(3:2+nfld),i,i) + ddh - 6._dl*diff_bg_hub_cosm(y)*y(h_i) - 4._dl*y(h_i)**3) &
                      /(cp_m2(i,i))*cp_l_m(i,i)/4._dl
      cp_c_dxdx_m(i,i) = bg_ddv(y(3:2+nfld),i,i) - diff_bg_hub_cosm(y)!cp_m2(i,i)
      !print*, 'cp_m2(i,i) = ', cp_m2(i,i)
    enddo
#endif
    cp_l_m = sqrt(2._dl*sqrt(k2_in))*cp_l_m
    cp_dl_m = sqrt(2._dl*sqrt(k2_in))*cp_dl_m
    !cp_c_dxdx_m = cp_c_dxdx_m/y(a_i)**2
    !print*, 'cp_c_dxdx_m(1,1) = ', cp_c_dxdx_m(1,1)
    !do i=1,nfld
    !  cp_c_dxdx_m(i,i) = cp_c_dxdx_m(i,i) - k2_in/y(a_i)**2 + 2._dl*y(h_i)**2
    !  print*, 'cp_c_dxdx_m(i,i) = ', cp_c_dxdx_m(i,i)
    !  print*, '- k2_in/y(a_i)**2 + 2._dl*y(h_i)**2 = ', - k2_in/y(a_i)**2 + 2._dl*y(h_i)**2
    !enddo
    k = 1
    do i=1,nfld
      do j=1,i
        cp_l(k) = cp_l_m(i,j)
        cp_dl(k) = cp_dl_m(i,j)
        cp_c_dxdx(k) = cp_c_dxdx_m(i,j)
        k = k+1
      enddo 
    enddo
  end subroutine set_cp_ic
#endif

  ! Subroutine to set initial conditions using the Mukhanov-Sasaki variable
  ! quantization for field 1.
  ! n.b. This method assumes the longitudinal DoF is purely the index 1 field
  !      when ICs are being set. This requires df_bg and dV to have 1 as the 
  !      only nonzero components.
  ! n.b. This method also assumes that ddV is diagonal when ICs are being set.
  ! to do: move the matrix stuff to another subroutine
  subroutine set_cp_ic_mp(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    
    integer :: i,j,k
    integer :: a_i, hub_i, df_i

    real(dl) :: dh, ddf

    real(dl), dimension(2,2) :: g
    real(dl), dimension(3,7) :: tr_a
    real(dl), dimension(7,3) :: tr_b
    real(dl), dimension(3,3) :: tr_c
    real(dl), dimension(3,1) :: corr_u
    real(dl), dimension(3,1) :: corr_v
    integer, dimension(3) :: piv
    integer, dimension(nfld) :: v_temp

    print*, 'sefault check 1'
    a_i = 1; hub_i = 2; df_i = 2+nfld+1
    cp_l = 0._dl; cp_dl = 0._dl; cp_c_xdx = 0._dl; cp_c_dxdx = 0._dl

    dh = diff_bg_hub_conf(y)
    ddf = diff_bg_df_conf(y, key=(/1/))  ! I require a key here

    print*, 'sefault check 2'
    g = 0._dl; tr_a = 0._dl; tr_b = 0._dl; tr_c = 0._dl
    g(1,1) = y(a_i)**2*bg_dv(y(3:2+nfld),1) + 3._dl*y(hub_i)*y(df_i)
    g(1,2) = y(df_i)
    g(2,1) = -y(hub_i)*y(a_i)**2*bg_dv(y(3:2+nfld),1) - (y(a_i)**2*bg_v(y(3:2+nfld))-k2_in)*y(df_i)
    g(2,2) = -y(hub_i)*y(df_i)
    g = g/(y(df_i)**2-2._dl*k2_in)
    print*, 'sefault check 3'
    tr_a(1,1) = 2._dl*y(df_i)/y(hub_i)
    tr_a(1,5) = (y(df_i)/y(hub_i))**2
    tr_a(2,1) = ddf/y(hub_i) - dh*y(df_i)/y(hub_i)**2 + 2._dl*y(df_i)
    tr_a(2,2) = y(df_i)/y(hub_i)
    tr_a(2,3) = y(df_i)/y(hub_i)
    tr_a(2,5) = (ddf/y(hub_i) - dh*y(df_i)/y(hub_i)**2 + y(df_i))*y(df_i)/y(hub_i)
    tr_a(2,6) = (y(df_i)/y(hub_i))**2
    tr_a(3,1) = 2._dl*(ddf-dh/y(hub_i)+y(hub_i)*y(df_i))
    tr_a(3,2) = 2._dl*y(hub_i)
    tr_a(3,3) = 2._dl*(ddf/y(hub_i) - dh*y(df_i)/y(hub_i)**2 + y(df_i))
    tr_a(3,4) = 2._dl*y(df_i)/y(hub_i)
    tr_a(3,5) = (ddf/y(hub_i) - dh*y(df_i)/y(hub_i)**2 + y(df_i))**2
    tr_a(3,6) = (ddf/y(hub_i) - dh*y(df_i)/y(hub_i)**2 + y(df_i))*y(df_i)/y(hub_i)
    tr_a(3,7) = (y(df_i)/y(hub_i))**2
    print*, 'sefault check 4'
    tr_b(1,1) = g(1,1); tr_b(1,2) = g(1,1)
    tr_b(2,1) = g(2,1); tr_b(2,2) = g(2,2)
    tr_b(3,2) = g(1,1); tr_b(3,3) = g(1,2)
    tr_b(4,2) = g(2,1); tr_b(4,3) = g(2,2)
    tr_b(5,1) = g(1,1)**2; tr_b(5,2) = 2._dl*g(1,1)*g(1,2); tr_b(5,3) = g(1,2)**2
    tr_b(6,1) = g(1,1)*g(2,1); tr_b(6,2) = g(1,1)*g(2,2)+g(1,2)*g(2,1); tr_b(6,3) = g(1,2)*g(2,2)
    tr_b(7,1) = g(2,1)**2; tr_b(7,2) = 2._dl*g(2,1)*g(2,2); tr_b(7,3) = g(2,2)**2
    print*, 'sefault check 5'
    tr_c(1,1) = 1._dl
    tr_c(2,1) = y(hub_i); tr_c(2,2) = 1._dl
    tr_c(3,1) = y(hub_i)**2; tr_c(3,2) = 2._dl*y(hub_i); tr_c(3,3) = 1._dl

    print*, 'sefault check 6'
    call dgemm('N','N',3,3,7,1._dl,tr_a,3,tr_b,7,1._dl,tr_c,3)
    ! Find <uu>, <uu'>, <u'u'>
    call get_cp_m2_ms_conf(y, k2_in)
    print*, "m2 with mp = ", cp_m2
    corr_v(1,1) = 1._dl/(2._dl*sqrt(cp_m2(1,1)))
    corr_v(2,1) = 0._dl
    corr_v(3,1) = sqrt(cp_m2(1,1))/2._dl
    corr_u = corr_v
    call dgesv(3,1,tr_c,3,piv,corr_u,3,i)  ! solve for <uu>, <uu'>, <u'u'>

    print*, 'sefault check 7'
    ! Set L
    cp_l_m = 0._dl    
    cp_l_m(1,1) = corr_u(1,1)
    do i=2,nfld
      cp_l_m(i,i) = 1._dl/(2._dl*sqrt(cp_m2(i,i) - (diff_bg_hub_conf(y) + y(hub_i)**2)))
    enddo
    call dpotrf('L',nfld,cp_l_m,nfld,i)
    !cp_l_m = cp_l_m/y(a_i)
    
    print*, 'sefault check 8'
    ! Set dL
    cp_dl_m = 0._dl
    cp_dl_m(1,1) = corr_u(2,1)
    do i=2,nfld
      cp_dl_m(i,i) = -y(hub_i)/(2._dl*sqrt(cp_m2(i,i) - (diff_bg_hub_conf(y) + y(hub_i)**2)))
    enddo
    call dgesv(nfld,nfld,cp_l_m,nfld,v_temp,cp_dl_m,nfld,i)
    cp_temp = cp_dl_m
    do i=1,nfld
      do j=1,nfld
        cp_dl_m(i,j) = cp_temp(j,i)
      enddo
    enddo

    print*, 'sefault check 9'
    ! Set c_dxdx
    cp_c_dxdx_m = 0._dl
    cp_c_dxdx_m(1,1) = corr_u(3,1)
    do i=2,nfld
      cp_c_dxdx_m(i,i) = y(hub_i)**2/(2._dl*sqrt(cp_m2(i,i) - (diff_bg_hub_conf(y) + y(hub_i)**2))) &
                         + sqrt(cp_m2(i,i) - (diff_bg_hub_conf(y) + y(hub_i)**2))/2._dl
    enddo
    call dgemm('N','T',nfld,nfld,nfld,-1._dl,cp_dl_m,nfld,cp_dl_m,nfld,1._dl,cp_c_dxdx_m,nfld)
    cp_temp = cp_l_m
    call dtrtri('L','N',nfld,cp_temp,nfld,i)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_dxdx_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_temp2,nfld,0._dl,cp_c_dxdx_m,nfld)

    print*, 'sefault check 10'
    ! Fill variables
    k = 1
    do i=1,nfld
      do j=1,i
        cp_l(k) = cp_l_m(i,j)
        cp_dl(k) = cp_dl_m(i,j)
        cp_c_dxdx(k) = cp_c_dxdx_m(i,j)
        k = k+1
      enddo 
    enddo
    print*, 'sefault check 11'
    ! Factors of a
    cp_l = cp_l/y(a_i); cp_dl = cp_dl/y(a_i)
  end subroutine set_cp_ic_mp

  ! Subroutine to set initial conditions to match the analytic solution for massless
  ! fields in de Sitter space.
  subroutine set_cp_ic_ds(y,k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    
    integer :: i,j,k
    integer :: a_i, h_i
    real(dl), dimension(nfld) :: diag

    a_i = 1; h_i = 2
    cp_l = 0._dl; cp_dl = 0._dl; cp_c_xdx = 0._dl; cp_c_dxdx = 0._dl
    cp_l_m = 0._dl; cp_dl_m = 0._dl; cp_c_dxdx_m = 0._dl

    do i=1,nfld
      cp_l_m(i,i) = sqrt((1._dl+y(h_i)**2/k2_in)/(2._dl*sqrt(k2_in)))/y(a_i)
      cp_dl_m(i,i) = -y(h_i)/sqrt(2._dl*sqrt(k2_in)*(1._dl+y(h_i)**2/k2_in))/y(a_i)
      cp_c_dxdx_m(i,i) = k2_in/(1._dl+y(h_i)**2/k2_in) - y(h_i)**2/(1._dl+y(h_i)**2/k2_in)**2
    enddo

    ! Fill variables
    k = 1
    do i=1,nfld
      do j=1,i
        cp_l(k) = cp_l_m(i,j)
        cp_dl(k) = cp_dl_m(i,j)
        cp_c_dxdx(k) = cp_c_dxdx_m(i,j)
        k = k+1
      enddo 
    enddo
  end subroutine set_cp_ic_ds

  ! Subroutine to test another way of setting initial conditions for the linear calculation.
  ! The way this goes is the adiabatic approximation is used to get a functional form for
  ! L, then that functional form is differentiated to set L'. We impose a condition on
  ! phase space volume to get c_dxdx.
  ! Set c_ff functional form, implies c_fdf, then phase space volume implies c_dfdf
  ! n.b. Starting off with just massless de Sitter before generalizing
  ! to do: 
  subroutine set_cp_ic_test(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    
    integer :: i,j,k
    integer :: a_i, h_i
    real(dl), dimension(nfld) :: diag

    a_i = 1; h_i = 2
    cp_l = 0._dl; cp_dl = 0._dl; cp_c_xdx = 0._dl; cp_c_dxdx = 0._dl
    cp_l_m = 0._dl; cp_dl_m = 0._dl; cp_c_dxdx_m = 0._dl

    do i=1,nfld
      cp_l_m(i,i) = 1._dl/sqrt(2._dl*sqrt(k2_in - 2._dl*y(h_i)**2))/y(a_i)
      cp_dl_m(i,i) = y(h_i)*((3._dl*y(h_i)**2 - k2_in)/(k2_in - 2._dl*y(h_i)**2))*cp_l_m(i,i)
      cp_c_dxdx_m(i,i) = k2_in - 2._dl*y(h_i)**2
    enddo

    ! Fill variables
    k = 1
    do i=1,nfld
      do j=1,i
        cp_l(k) = cp_l_m(i,j)
        cp_dl(k) = cp_dl_m(i,j)
        cp_c_dxdx(k) = cp_c_dxdx_m(i,j)
        k = k+1
      enddo 
    enddo

  end subroutine set_cp_ic_test
  
  ! Function to calculate conformal time derivative of cp_l. key input is used
  ! to determine which component is calculated.
!  function diff_cp_l_conf(y, ode_param, key)
!    real(dl) :: diff_cp_l_conf
!    real(dl) :: y(:)
!    real(dl), pointer, optional :: ode_param(:)
!    integer, optional :: key(:)

!    integer :: i,j
!    i = key(1); j = key(2)
!    call set_cp_dl_m(y)
!    diff_cp_l_conf = cp_dl_m(i,j)
!  end function diff_cp_l_conf

  ! Function to calculate the derivative of cp_l for the (a delta phi)=L chi equations
  ! in conformal time.
  ! Function to calculate conformal time derivative of cp_l. key input is used
  ! to determine which component is calculated.
  function diff_cp_l_conf(y, ode_param, key)
  !function diff_cp_l_ascl_conf(y, ode_param, key)
    real(dl) :: diff_cp_l_conf  !diff_cp_l_ascl_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    call set_cp_dl_m(y)
    diff_cp_l_conf = cp_dl_m(key(1),key(2))  !diff_cp_l_ascl_conf = cp_dl_m(key(1),key(2))
  !end function diff_cp_l_ascl_conf
  end function diff_cp_l_conf

  ! Function to calculate conformal time derivative of cp_dl. key input is used
  ! to determine which component is calculated.
#if (CPVAROPT == 0)
  function diff_cp_dl_conf(y, ode_param, key)
    real(dl) :: diff_cp_dl_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i,j
    real(dl) :: cp_ddl(nfld,nfld)

    i = key(1); j = key(2)
    cp_ddl = 0._dl
    call set_cp_constraint(y, cp_k2)
    call set_cp_l_m(y)
    call set_cp_dl_m(y)
    call set_cp_c_xdx_m(y)
    call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y)
    call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_a,nfld,0._dl,cp_ddl,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_ddl,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_c_dxdx_m,nfld,1._dl,cp_ddl,nfld)
    !print*, 'cp_ddl(i,j) = ', cp_ddl(i,j)
    !print*, 'cp_b(i,j) = ', cp_b(i,j)
    cp_ddl = cp_ddl - cp_b
    diff_cp_dl_conf = cp_ddl(i,j)
  end function diff_cp_dl_conf
#endif

  ! Function to calculate the derivative of cp_dl for the (a delta phi)=L chi equations
  ! in conformal time.
#if (CPVAROPT == 1)
  function diff_cp_dl_conf(y, ode_param, key)
  !function diff_cp_dl_ascl_conf(y, ode_param, key)
    real(dl) :: diff_cp_dl_conf!diff_cp_dl_ascl_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

#if (MODE_MP == 1)
    call set_cp_constraint_ascl(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)!call get_cp_m2_ascl_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)!call get_cp_b_ascl_conf(y); call get_cp_c_ascl_conf(y)
    call get_cp_a_conf(y)!call get_cp_a_ascl_conf(y)

    cp_temp = -cp_b
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_c_dxdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    diff_cp_dl_conf = cp_temp(key(1),key(2))!diff_cp_dl_ascl_conf = cp_temp(key(1),key(2))
  !end function diff_cp_dl_ascl_conf
  end function diff_cp_dl_conf
#endif

  ! Function to calculate the derivative of cp_dl for the (a delta phi)=L chi equations
  ! with c_dxdx change of variables in conformal time.
#if (CPVAROPT == 2)
  function diff_cp_dl_conf(y, ode_param, key)
  !function diff_cp_dl_ascl_conf(y, ode_param, key)
    real(dl) :: diff_cp_dl_conf!diff_cp_dl_ascl_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: h_i

    h_i = 2
#if (MODE_MP == 1)
    call set_cp_constraint_ascl(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    cp_temp = -cp_b
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,-2._dl*y(h_i)**2,cp_l_m,nfld,cp_c_dxdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    diff_cp_dl_conf = cp_temp(key(1),key(2))!diff_cp_dl_ascl_conf = cp_temp(key(1),key(2))
  !end function diff_cp_dl_ascl_conf
  end function diff_cp_dl_conf
#endif

#if (CPVAROPT == 3)
  function diff_cp_dl_conf(y, ode_param, key)
    real(dl) :: diff_cp_dl_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

#if (MODE_MP == 1)
    call set_cp_constraint_ascl(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    !call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    cp_temp = -cp_b
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_c_dxdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    diff_cp_dl_conf = cp_temp(key(1),key(2))
  end function diff_cp_dl_conf
#endif

  ! Function to calculate the conformal time derivative of the correlation of
  ! xdx.
!  function diff_cp_c_xdx_conf(y, ode_param, key)
!    real(dl) :: diff_cp_c_xdx_conf
!    real(dl) :: y(:)
!    real(dl), pointer, optional :: ode_param(:)
!    integer, optional :: key(:)
!
!    integer :: i,j
!
!    i = key(1); j = key(2)
!    call set_cp_constraint(y, cp_k2)
!    call set_cp_l_m(y)
!    call set_cp_dl_m(y)
!    call set_cp_c_xdx_m(y)
!    call set_cp_c_dxdx_m(y)
!    call get_cp_m2_conf(y, cp_k2)
!    call get_cp_b_conf(y)
!    call get_cp_c_conf(y)
!    call get_cp_a_conf(y)
!
!    diff_cp_c_xdx_conf = cp_a(i,j)
!  end function diff_cp_c_xdx_conf

  ! Function to calculate the conformal time derivative of the correlation of
  ! xdx for the (a delta phi)=L chi equations in conformal time.
  ! Function to calculate the conformal time derivative of the correlation of
  ! xdx.
  function diff_cp_c_xdx_conf(y, ode_param, key)
    real(dl) :: diff_cp_c_xdx_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

#if (MODE_MP == 1)
    call set_cp_constraint(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    diff_cp_c_xdx_conf = cp_a(key(1),key(2))
  end function diff_cp_c_xdx_conf

  ! Function to calculate the conformal time derivative of the correlation of
  ! dxdx.
  ! n.b. I've taken a transpose inside Sym[] in this evaluation.
#if (CPVAROPT == 0)
  function diff_cp_c_dxdx_conf(y, ode_param, key)
    real(dl) :: diff_cp_c_dxdx_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i,j,l
    real(dl) :: cp_dc_dxdx(nfld,nfld)

    i = key(1); j = key(2)
    cp_dc_dxdx = 0._dl
    call set_cp_constraint(y, cp_k2)
    call set_cp_l_m(y)
    call set_cp_dl_m(y)
    call set_cp_c_xdx_m(y)
    call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y)
    call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    !call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx_m,nfld,0._dl,cp_dc_dxdx,nfld)
    !cp_dc_dxdx = cp_dc_dxdx + cp_c_dxdx_m
    !call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_dc_dxdx,nfld,0._dl,cp_dc_dxdx,nfld)  !!
    !call dtrtri('L','N',nfld,cp_l_m,nfld,l)
    !call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_l_m,nfld,cp_dc_dxdx,nfld,0._dl,cp_dc_dxdx,nfld) !!
    !call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_dxdx_m,nfld,1._dl,cp_dc_dxdx,nfld) !!
    !call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_c_xdx_m,nfld,cp_a,nfld,1._dl,cp_dc_dxdx,nfld)
    !diff_cp_c_dxdx_conf = cp_dc_dxdx(i,j) + cp_dc_dxdx(j,i)

    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx_m,nfld,0._dl,cp_temp,nfld)
    cp_temp = cp_temp + cp_c_dxdx_m
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)  !!
    call dtrtri('L','N',nfld,cp_l_m,nfld,l)
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_l_m,nfld,cp_temp2,nfld,0._dl,cp_temp,nfld) !!
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_dxdx_m,nfld,1._dl,cp_temp,nfld) !!
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_c_xdx_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    diff_cp_c_dxdx_conf = cp_temp(i,j) + cp_temp(j,i)
  end function diff_cp_c_dxdx_conf
#endif

  ! Function to calculate the derivative of cp_c_dxdx for the (a delta phi)=L chi equations
  ! in conformal time.
#if (CPVAROPT == 1)
  function diff_cp_c_dxdx_conf(y, ode_param, key)
  !function diff_cp_c_dxdx_ascl_conf(y, ode_param, key)
    real(dl) :: diff_cp_c_dxdx_conf!diff_cp_c_dxdx_ascl_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i,j,l

#if (MODE_MP == 1)
    call set_cp_constraint(y, cp_k2)!call set_cp_constraint_ascl(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)!call get_cp_m2_ascl_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)!call get_cp_b_ascl_conf(y); call get_cp_c_ascl_conf(y)
    call get_cp_a_conf(y)!call get_cp_a_ascl_conf(y)

    cp_temp = -cp_c_dxdx_m
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dtrtri('L','N',nfld,cp_l_m,nfld,i)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp2,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    
    diff_cp_c_dxdx_conf = cp_temp(key(1),key(2)) + cp_temp(key(2),key(1))!diff_cp_c_dxdx_ascl_conf = cp_temp(key(1),key(2)) + cp_temp(key(2),key(1))
  !end function diff_cp_c_dxdx_ascl_conf
  end function diff_cp_c_dxdx_conf
#endif

  ! Function to calculate the derivative of cp_c_dxdx for the (a delta phi)=L chi equations
  ! and c_dxdx change of variables in conformal time.
#if (CPVAROPT == 2)
  function diff_cp_c_dxdx_conf(y, ode_param, key)
    real(dl) :: diff_cp_c_dxdx_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i,j,l
    integer :: h_i

    h_i = 2
#if (MODE_MP == 1)
    call set_cp_constraint(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    cp_temp = 2._dl*y(h_i)**2*cp_c_dxdx_m
    do i=1,nfld
      cp_temp(i,i) = cp_temp(i,i) - cp_k2
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dtrtri('L','N',nfld,cp_l_m,nfld,i)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp2,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    
    diff_cp_c_dxdx_conf = -(cp_temp(key(1),key(2)) + cp_temp(key(2),key(1)))/(2._dl*y(h_i)**2) &
                          - 2._dl*diff_bg_hub_conf(y)*cp_c_dxdx_m(key(1),key(2))/y(h_i)
  end function diff_cp_c_dxdx_conf
#endif

#if (CPVAROPT == 3)
  function diff_cp_c_dxdx_conf(y, ode_param, key)
    real(dl) :: diff_cp_c_dxdx_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i,j,l
    integer :: a_i, h_i

    a_i = 1; h_i = 2
#if (MODE_MP == 1)
    call set_cp_constraint(y, cp_k2)
#endif
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    !call get_cp_m2_conf(y, cp_k2)
    call get_cp_b_conf(y); call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    cp_temp = -cp_c_dxdx_m
    do i=1,nfld
      cp_temp(i,i) = cp_temp(i,i) - cp_k2/y(a_i)**2 + 2._dl*y(h_i)**2
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dtrtri('L','N',nfld,cp_l_m,nfld,i)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp2,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_a,nfld,1._dl,cp_temp,nfld)
    
    cp_temp2 = 0._dl
    do i=1,nfld
      cp_temp2(i,i) = 2._dl*y(h_i)*cp_k2/y(a_i)**2 + 4._dl*diff_bg_hub_cosm(y)*y(h_i)
    enddo
    diff_cp_c_dxdx_conf = (cp_temp(key(1),key(2)) + cp_temp(key(2),key(1))) &
                          - cp_c_dxdx_m(key(1),key(2)) + cp_temp2(key(1),key(2))
  end function diff_cp_c_dxdx_conf
#endif

  ! Subroutine to perform a single time step of the cp integration
  subroutine step_cp(dt)
    real(dl) :: dt

    call gl_integrate(cp_sys_temp, cp_ode, cp_sys_g, dt, cp_sys_dim)
  end subroutine step_cp

  ! to do: add key reference
  subroutine evolve_cp(y_end, ind, dt)
    real(dl) :: y_end
    integer :: ind
    real(dl) :: dt

    integer :: i
    real(dl) :: diff
 
    !print*, 'in evolve_cp 0'
    call set_cp_sys_temp()
    diff = cp_ode(ind)%f(cp_sys_temp, ode_param=cp_ode(ind)%ode_param, key=cp_ode(ind)%key)
    !print*, 'in evolve_cp 1'
    if (sign(1._dl,diff)*(cp_ode(ind)%y-y_end) .LE. 0._dl) then
      !print*, 'in evolve_cp 2'
      do while (sign(1._dl,diff)*(cp_ode(ind)%y-y_end) .LE. 0._dl)
        !print*, 'in evolve_cp 2a'
        call step_cp(dt)
      enddo
    else
      !print*, 'in evolve_cp 3'
      do while (sign(1._dl,diff)*(cp_ode(ind)%y-y_end) .GE. 0._dl)
        !print*, 'In backstep'
        call step_cp(-dt)
      enddo
    endif
    !print*, 'in evolve_cp 4'
    call gl_newton_root(cp_sys_temp, cp_ode, cp_sys_g, dt, cp_sys_dim, y_end, ind)
  end subroutine evolve_cp

  ! Subroutine to set cp_corr for one mode
#if (CPVAROPT == 0)
  subroutine set_cp_corr(y, k2_in, k_ind)
    real(dl) :: y(:)
    real(dl) :: k2_in
    integer :: k_ind  ! oversampled mode number

    integer :: i,j

    call set_cp_constraint(y, k2_in)
    call set_cp_l_m(y)
    call set_cp_dl_m(y)
    call set_cp_c_xdx_m(y)
    call set_cp_c_dxdx_m(y)
    call get_cp_m2_conf(y, k2_in)
    call get_cp_b_conf(y)
    call get_cp_c_conf(y)
    call get_cp_a_conf(y)

    cp_temp = cp_l_m
    call dtrmm('R','L','T','N',nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld)
    cp_corr(1:nfld,1:nfld,k_ind) = cp_temp

    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_l_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_dl_m,nfld,1._dl,cp_temp2,nfld)
    cp_corr(1:nfld,nfld+1:2*nfld,k_ind) = cp_temp2

    do i=1,nfld
      do j=1,nfld
        cp_corr(nfld+i,j,k_ind) = cp_corr(j,nfld+i,k_ind)
      enddo    
    enddo

    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_dxdx_m,nfld,cp_l_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_temp2
    call dgemm('N','T',nfld,nfld,nfld,-1._dl,cp_c_xdx_m,nfld,cp_dl_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) + cp_temp2
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_l_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_dl_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) + cp_temp2
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_dl_m,nfld,cp_dl_m,nfld,0._dl,cp_temp,nfld)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) + cp_temp
  end subroutine set_cp_corr
#endif

  ! Subroutine to 2pt functions c_ff, c_fdf, c_dfdf for the (a delta phi)=L chi equations
  ! in conformal time, optionally including the change of variables to c_dxdx
#if (CPVAROPT == 1 .OR. CPVAROPT == 2 .OR. CPVAROPT == 3)
  subroutine set_cp_corr(y, k2_in, k_ind)
  !subroutine set_cp_corr_ascl(y, k2_in, k_ind)
    real(dl) :: y(:)
    real(dl) :: k2_in
    integer :: k_ind  ! oversampled mode number

    integer :: i
    integer :: a_i, h_i

    a_i = 1; h_i = 2
    call set_cp_l_m(y); call set_cp_dl_m(y); call set_cp_c_xdx_m(y); call set_cp_c_dxdx_m(y)
    
#if (CPVAROPT == 2)
    ! modify variables for use in proceedure below
    cp_l_m = cp_l_m/sqrt(2._dl*sqrt(k2_in))
    cp_dl_m = cp_dl_m/sqrt(2._dl*sqrt(k2_in))
    cp_c_dxdx_m = -2._dl*y(h_i)**2*cp_c_dxdx_m
    do i=1,nfld
      cp_c_dxdx_m(i,i) = cp_c_dxdx_m(i,i) + k2_in
    enddo
#elif (CPVAROPT == 3)
    ! modify variables for use in proceedure below
    cp_l_m = cp_l_m/sqrt(2._dl*sqrt(k2_in))
    cp_dl_m = y(a_i)*cp_dl_m/sqrt(2._dl*sqrt(k2_in))
    cp_c_dxdx_m = y(a_i)**2*cp_c_dxdx_m
    do i=1,nfld
      cp_c_dxdx_m(i,i) = cp_c_dxdx_m(i,i) + k2_in - 2._dl*y(a_i)**2*y(h_i)**2
    enddo
#endif

    cp_temp = cp_l_m
    call dtrmm('R','L','T','N',nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld)
    cp_corr(1:nfld,1:nfld,k_ind) = cp_temp

    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_l_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_dl_m,nfld,1._dl,cp_temp2,nfld)
    cp_corr(1:nfld,nfld+1:2*nfld,k_ind) = cp_temp2 - y(h_i)*cp_corr(1:nfld,1:nfld,k_ind)
    
    cp_temp = cp_corr(1:nfld,nfld+1:2*nfld,k_ind)
    cp_corr(nfld+1:2*nfld,1:nfld,k_ind) = transpose(cp_temp)

    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_dl_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_l_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_dl_m,nfld,cp_dl_m,nfld,1._dl,cp_temp2,nfld)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_temp2
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_xdx_m,nfld,cp_l_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_dl_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) + cp_temp2
    call dgemm('N','T',nfld,nfld,nfld,1._dl,cp_c_dxdx_m,nfld,cp_l_m,nfld,0._dl,cp_temp,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_l_m,nfld,cp_temp,nfld,0._dl,cp_temp2,nfld)
    cp_temp = cp_corr(1:nfld,nfld+1:2*nfld,k_ind)
    cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) = cp_corr(nfld+1:2*nfld,nfld+1:2*nfld,k_ind) + cp_temp2 &
                                                 - y(h_i)**2*cp_corr(1:nfld,1:nfld,k_ind) &
                                                 - y(h_i)*cp_temp - y(h_i)*transpose(cp_temp)
    
    cp_corr = cp_corr/y(a_i)**2
  !end subroutine set_cp_corr_ascl
  end subroutine set_cp_corr
#endif

  ! Subroutine to calculate the variance of the Mukhanov-Sasaki variable
  ! n.b. This is currently only for the 1field definintion of the Mukhanov-
  !      Sasaki variable. I'm using it for debugging, but will need to
  !      extentend to the multifield case.
  ! to do: extend to the multifield case
  subroutine set_cp_corr_muksas(y, k2_in, k_ind)
    real(dl) :: y(:)
    real(dl) :: k2_in
    integer :: k_ind

    integer :: i,j
    integer :: a_i, hub_i, df_i
 
    a_i = 1; hub_i = 2; df_i = 2+nfld+1

    call set_cp_corr(y, k2_in, k_ind)

    cp_corr_muksas = y(a_i)**2*&
                     ((1._dl+y(df_i)*cp_constraint(1,1)/y(hub_i))**2*cp_corr(1,1,k_ind)&
                      +2._dl*y(df_i)*cp_constraint(1,nfld+1)/y(hub_i)*(1._dl+y(df_i)*cp_constraint(1,1)/y(hub_i))*cp_corr(1,nfld+1,k_ind) &
                      +(y(df_i)*cp_constraint(1,nfld+1)/y(hub_i))**2*cp_corr(nfld+1,nfld+1,k_ind))
  end subroutine set_cp_corr_muksas

  ! Subroutine to initialize cp_corr
  subroutine init_cp_corr_old(a_in, hub_in, f_in, df_in, k_phys, dt)
    real(dl), intent(in) :: a_in, hub_in
    real(dl), intent(in) :: f_in(:), df_in(:)
    real(dl) :: k_phys  ! k/aH at mode initialization
    real(dl) :: dt      ! time step

    integer :: i, k
    integer :: a_i, h_i, f_i

    a_i = 1; h_i = 2; f_i = 3
    do k=1,cp_nos
      cp_k2 = (k*cp_dkos)**2
      call set_bg_cosmo_ic(a_in, hub_in, f_in, df_in)
      call evolve_bg_cosmo(k_phys, h_i, dt)
      call set_cp_ic((/a_bg/),cp_k2)
      call evolve_cp(f_in(1), f_i, dt)  ! ending on phi value
      do i=1,cp_sys_dim
        cp_sys_temp(i) = cp_ode(i)%y
      enddo
      call set_cp_corr(cp_sys_temp, cp_k2, k)
    enddo
  end subroutine init_cp_corr_old

  ! Subroutine to initialize cp corr
  subroutine init_cp_corr(y0_in, kah_fac, dt_in, end_ind)
    real(dl), intent(in) :: y0_in(:)  ! background variables at lattice initialization
    real(dl) :: kah_fac              ! 
    real(dl), intent(in) :: dt_in    ! time step
    integer :: end_ind               ! index of variable used to end evolution

    integer :: k
    integer :: h_i

    h_i = 2
    ! Loop over sampled modes
    do k=1,cp_nos
      print*, 'init_cp_corr mode: '
      cp_k2 = (k*cp_dkos)**2
      print*, 'cp_k2 = ', cp_k2
      ! evolve bg to kah_fac = k/aH
      call evolve_bg_cosmo(k*cp_dkos/kah_fac, h_i, dt_in)
      ! set cp ICs
      call set_bg_sys_temp()
      print*, 'background evolution'
      call set_cp_ic(bg_sys_temp, cp_k2)
      print*, 'perturbation evolution'
      print*, 'y0_in = ', y0_in
      print*, 'y0_in(end_ind) = ', y0_in(end_ind)
      call evolve_cp(y0_in(end_ind), end_ind, dt_in)
      print*, 'perturbation evolution complete'
      call set_cp_sys_temp()
      call set_cp_corr(cp_sys_temp, cp_k2, k)
    enddo
  end subroutine init_cp_corr

  ! Subroutine to set the effetive mass matrix in conformal time
  ! n.b. This is setting k^2 + a^2V_,ff
#if (CPVAROPT == 0)
  subroutine get_cp_m2_conf(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in

    integer :: i,j
    integer :: a_i

    a_i = 1
    do i=1,nfld
      do j=1,nfld
        !cp_m2(i,j) = y(a_i)**2*pot_ddv(y(3),0._dl,i,j)
        cp_m2(i,j) = y(a_i)**2*bg_ddv(y(3:2+nfld),i,j)
      enddo 
    enddo
    do i=1,nfld
      cp_m2(i,i) = cp_m2(i,i) + k2_in 
    enddo
  end subroutine get_cp_m2_conf
#endif

  ! Subroutine to set the effective mass matrix in conformal time
  ! n.b. This is setting k^2 + a^2V_,ff -a''/a
#if (CPVAROPT == 1 .OR. CPVAROPT == 2 .OR. CPVAROPT == 3)
  subroutine get_cp_m2_conf(y, k2_in)
  !subroutine get_cp_m2_ascl_conf(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in

    integer :: i,j
    integer :: a_i, h_i

    a_i = 1; h_i = 2
    do i=1,nfld
      do j=1,nfld
        cp_m2(i,j) = y(a_i)**2*bg_ddv(y(3:2+nfld),i,j)
      enddo
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
      cp_m2(i,i) = cp_m2(i,i) + k2_in - y(h_i)**2 - diff_bg_hub_conf(y)
#elif (CPVAROPT == 3)
      cp_m2(i,i) = cp_m2(i,i) + k2_in - y(a_i)**2*(2._dl*y(h_i)**2 + diff_bg_hub_cosm(y))
#endif
    enddo
  !end subroutine get_cp_m2_ascl_conf
  end subroutine get_cp_m2_conf
#endif

  ! Subroutine to set the effective mass matrix in conformal time with metric
  ! perturbations and the longitudinal DoF purely in the 1 index direction.
  subroutine get_cp_m2_ms_conf(y, k2_in)
    real(dl) :: y(:)
    real(dl) :: k2_in
    real(dl) :: dh

    integer :: i,j
    integer :: a_i, hub_i, df_i

    a_i = 1; hub_i = 2; df_i = 2+nfld+1
    call get_cp_m2_conf(y, k2_in)
    dh = diff_bg_hub_conf(y)
    ! -z''/z
    cp_m2(1,1) = cp_m2(1,1) + 3._dl*(y(hub_i)**2-dh) &
                 - 2._dl*(dh**2/y(hub_i)**2 + y(a_i)**2*bg_dv(y(3:2+nfld),1)*y(df_i)/y(hub_i))
  end subroutine get_cp_m2_ms_conf

  ! Function to get the anti-symmetric matrix defined in the conformal time cp
  ! 2pt corr system
#if (CPVAROPT == 0)
  subroutine get_cp_a_conf(y)
    real(dl) :: y(:)
    integer :: i,j,l

    cp_temp = cp_l_m
    call dtrtri('L','N',nfld,cp_temp,nfld,l)
    
    !call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_c,nfld,cp_c_xdx_m,nfld,0._dl,cp_temp2,nfld)
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_c,nfld,cp_c_xdx_m,nfld,0._dl,cp_temp2,nfld)
    cp_temp2 = cp_temp2 + cp_b
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_temp2,nfld,0._dl,cp_a,nfld)
    cp_a = cp_a - cp_c_dxdx_m

    do i=1,nfld
      do j=1,i-1
        cp_a(i,j) = -cp_a(j,i)
      enddo
      cp_a(i,i) = 0._dl
    enddo
  end subroutine get_cp_a_conf
#endif

  ! Function to get the anti-symmetric matrix defined in the conformal time cp
  ! 2pt corr system using (a delta phi) = L chi.
  ! n.b. Call get_cp_b_ascl_conf, and get_cp_c_ascl_conf first.
#if (CPVAROPT == 1)
  subroutine get_cp_a_conf(y)
  !subroutine get_cp_a_ascl_conf(y)
    real(dl) :: y(:)

    integer :: i,j

    cp_a = -cp_c_dxdx_m
    cp_temp = cp_l_m
    call dtrtri('L','N',nfld,cp_temp,nfld,i)
    cp_temp2 = cp_b
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp2,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_temp2,nfld,1._dl,cp_a,nfld)
    
    do i=1,nfld
      do j=1,i-1
        cp_a(i,j) = -cp_a(j,i)
      enddo
      cp_a(i,i) = 0._dl
    enddo
  !end subroutine get_cp_a_ascl_conf
  end subroutine get_cp_a_conf
#endif

  ! Function to get the anti-symmetric matrix defined in the conformal time cp
  ! 2pt corr system using (a delta phi) = L chi and c_dxdx variable change
  ! n.b. Call get_cp_b_ascl_conf, and get_cp_c_ascl_conf first.
#if (CPVAROPT == 2)
  subroutine get_cp_a_conf(y)
  !subroutine get_cp_a_ascl_conf(y)
    real(dl) :: y(:)

    integer :: i,j
    integer :: h_i

    h_i = 2
    cp_a = 2._dl*y(h_i)**2*cp_c_dxdx_m
    cp_temp = cp_l_m
    call dtrtri('L','N',nfld,cp_temp,nfld,i)
    cp_temp2 = cp_b
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp2,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_temp2,nfld,1._dl,cp_a,nfld)
    
    do i=1,nfld
      do j=1,i-1
        cp_a(i,j) = -cp_a(j,i)
      enddo
      cp_a(i,i) = 0._dl
    enddo
  !end subroutine get_cp_a_ascl_conf
  end subroutine get_cp_a_conf
#endif

#if (CPVAROPT == 3)
  subroutine get_cp_a_conf(y)
    real(dl) :: y(:)

    integer :: i,j

    cp_a = - cp_c_dxdx_m
    cp_temp = cp_l_m
    call dtrtri('L','N',nfld,cp_temp,nfld,i)
    cp_temp2 = cp_b
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cp_c,nfld,cp_c_xdx_m,nfld,1._dl,cp_temp2,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_temp2,nfld,1._dl,cp_a,nfld)
    
    do i=1,nfld
      do j=1,i-1
        cp_a(i,j) = -cp_a(j,i)
      enddo
      cp_a(i,i) = 0._dl
    enddo
  end subroutine get_cp_a_conf
#endif

  ! Function to get the B matrix defined in the conformal time cp 2pt corr
  ! system
  ! n.b. First call set_cp_constraint, set_cp_l_m, set_cp_dl_m, set_cp_m2
#if (CPVAROPT == 0)
  subroutine get_cp_b_conf(y)
    real(dl) :: y(:)
    !real(dl) :: k2_in

    integer :: i
    integer :: a_i, h_i!, f_i, df_i
    real(dl), dimension(nfld,nfld) :: temp_m

    a_i = 1; h_i = 2!; f_i = 3; df_i = 4
    cp_b = 2._dl*y(h_i)*cp_dl_m
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_m2,nfld,cp_l_m,nfld,1._dl, &
               cp_b,nfld)
#if (MODE_MP == 1)
    do i=1,nfld
      temp_m(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,1:nfld) &
                    +2._dl*y(a_i)**2*bg_dv(y(3:2+nfld),i)*cp_constraint(1,1:nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp_m,nfld,cp_l_m,nfld,1._dl, &
               cp_b,nfld)
    do i=1,nfld
      temp_m(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,nfld+1:2*nfld) &
                    +2._dl*y(a_i)**2*bg_dv(y(3:2+nfld),i)*cp_constraint(1,nfld+1:2*nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp_m,nfld,cp_dl_m,nfld,1._dl, &
               cp_b,nfld)
    !do i=1,nfld
    !  cp_b(i,:) = cp_b(i,:)&
    !              +(-4._dl*y(df_i)*cp_constraint(2,1)&
    !                +2._dl*y(a_i)**2*bg_dv(y(f_i))*cp_constraint(1,1)&
    !               )*cp_l_m(1,:)&
    !              +(2._dl*y(a_i)**2*bg_dv(y(f_i))*cp_constraint(1,2))*cp_dl_m(1,:)
    !enddo
#endif
  end subroutine get_cp_b_conf
#endif

  ! Function to get the B matrix defined in the conformal time cp 2pt corr
  ! system using (a delta phi) = L chi.
  ! n.b. Call subroutines to set cp_l_m, cp_dl_m, cp_xdx_m, cp_dxdx_m, and set_cp_m2 first.
  ! n.b. If using metric perturbations also call set_cp_constraint_ascl first.
  ! to do: take k2 as an input
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
  subroutine get_cp_b_conf(y)
  !subroutine get_cp_b_ascl_conf(y)
    real(dl) :: y(:)

    integer :: i
    integer :: a_i, h_i

    a_i = 1; h_i = 2
    cp_temp = cp_m2
#if (CPVAROPT == 2)
    do i=1,nfld
      cp_temp(i,i) = cp_temp(i,i) - cp_k2
    enddo
#endif
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_l_m,nfld,0._dl,cp_b,nfld)
#if (MODE_MP == 1)
    do i=1,nfld
      cp_temp(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,1:nfld) &
                     + (2._dl*y(a_i)**2*bg_dv(y(3:2+nfld),i) + 4._dl*y(h_i)*y(2+nfld+i))*cp_constraint(1,1:nfld)
      cp_temp2(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,nfld+1:2*nfld) &
                     + (2._dl*y(a_i)**2*bg_dv(y(3:2+nfld),i) + 4._dl*y(h_i)*y(2+nfld+i))*cp_constraint(1,nfld+1:2*nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_l_m,nfld,1._dl,cp_b,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp2,nfld,cp_dl_m,nfld,1._dl,cp_b,nfld)
#endif
  !end subroutine get_cp_b_ascl_conf
  end subroutine get_cp_b_conf
#endif

#if (CPVAROPT == 3)
  subroutine get_cp_b_conf(y)
    real(dl) :: y(:)

    integer :: i,j
    integer :: a_i, h_i

    a_i = 1; h_i = 2
    do i=1,nfld
      do j=1,nfld
        cp_temp(i,j) = bg_ddv(y(3:2+nfld),i,j)
      enddo
      cp_temp(i,i) = cp_temp(i,i) - diff_bg_hub_cosm(y)
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_l_m,nfld,0._dl,cp_b,nfld)
    cp_b = cp_b - y(h_i)*cp_dl_m
#if (MODE_MP == 1)
    do i=1,nfld
      cp_temp(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,1:nfld) &
                     + (2._dl*bg_dv(y(3:2+nfld),i) + 4._dl*y(h_i)*y(2+nfld+i))*cp_constraint(1,1:nfld)
      cp_temp2(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,nfld+1:2*nfld) &
                     + (2._dl*bg_dv(y(3:2+nfld),i) + 4._dl*y(h_i)*y(2+nfld+i))*cp_constraint(1,nfld+1:2*nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp,nfld,cp_l_m,nfld,1._dl,cp_b,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cp_temp2,nfld,cp_dl_m,nfld,1._dl,cp_b,nfld)
#endif
  !end subroutine get_cp_b_ascl_conf
  end subroutine get_cp_b_conf
#endif

  ! Function to get the C matrix defined in the conformal time cp 2pt corr
  ! system
  ! n.b. Frist call set_cp_constraint, set_cp_l_m, set_cp_dl_m
#if (CPVAROPT == 0)
  subroutine get_cp_c_conf(y)
    real(dl) :: y(:)

    integer :: i
    integer :: a_i, h_i!, f_i, df_i
    real(dl), dimension(nfld,nfld) :: temp_m

    a_i = 1; h_i = 2!; f_i = 3; df_i = 4
    cp_c = 2._dl*(cp_dl_m + y(h_i)*cp_l_m)
#if (MODE_MP == 1)
    do i=1,nfld
      temp_m(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,nfld+1:2*nfld) &
                    +2._dl*y(a_i)**2*bg_dv(y(3:2+nfld),i)*cp_constraint(1,nfld+1:2*nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp_m,nfld,cp_l_m,nfld,1._dl, &
               cp_c,nfld)
    !do i=1,nfld
    !  cp_c(i,:) = cp_c(i,:)&
    !              +2._dl*y(a_i)**2*bg_dv(y(f_i))*cp_constraint(1,2)*cp_l_m(1,:)
    !enddo
#endif
  end subroutine get_cp_c_conf
#endif

  ! Function to get the C matrix defined in the conformal time cp 2pt corr
  ! system using (a delta phi) = L chi.
  ! n.b. Call subroutines to set cp_l_m, cp_dl_m, cp_xdx_m, cp_dxdx_m first.
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
  subroutine get_cp_c_conf(y)
  !subroutine get_cp_c_ascl_conf(y)
    real(dl) :: y(:)

    integer :: i
    integer :: a_i, h_i

    a_i = 1; h_i = 2
    cp_c = 2._dl*cp_dl_m
#if (MODE_MP == 1)
    do i=1,nfld
      cp_temp(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,nfld+1:2*nfld) &
                     + (2._dl*y(a_i)**2*bg_dv(y(3:2+nfld),i) + 4._dl*y(h_i)*y(2+nfld+i))*cp_constraint(1,nfld+1:2*nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nlfd,1._dl,cp_temp,nfld,cp_l_m,nfld,1._dl,cp_c,nfld)
#endif
  !end subroutine get_cp_c_ascl_conf
  end subroutine get_cp_c_conf
#endif

#if (CPVAROPT == 3)
  subroutine get_cp_c_conf(y)
    real(dl) :: y(:)

    integer :: i
    integer :: h_i

    h_i = 2
    cp_c = 2._dl*cp_dl_m + y(h_i)*cp_l_m
#if (MODE_MP == 1)
    do i=1,nfld
      cp_temp(i,:) = -4._dl*y(2+nfld+i)*cp_constraint(2,nfld+1:2*nfld) &
                     + (2._dl*bg_dv(y(3:2+nfld),i) + 4._dl*y(h_i)*y(2+nfld+i))*cp_constraint(1,nfld+1:2*nfld)
    enddo
    call dgemm('N','N',nfld,nfld,nlfd,1._dl,cp_temp,nfld,cp_l_m,nfld,1._dl,cp_c,nfld)
#endif
  end subroutine get_cp_c_conf
#endif

  ! Subroutine to calculate the constraint matrix to calculate metric
  ! perturbations from field perturbations in conformal time and longitudinal
  ! gauge.
  ! n.b. This is for conformal time
  ! n.b. LAPACK use (rows,columns) convension
#if (CPVAROPT == 0)
  subroutine set_cp_constraint(y, k2_in)  
    real(dl) :: y(:)  ! array of background variables
    real(dl) :: k2_in

    integer :: i
    integer :: a_i, h_i, f_i, df_i
    real(dl), dimension(2,2) :: temp_m
    real(dl), dimension(2,2*nfld) :: temp_constraint
    real(dl) :: rho_kin

    a_i = 1; h_i = 2!; f_i = 3; df_i = 4

    ! set 2x2 components
    temp_m(1,1) = 1._dl; temp_m(1,2) = 3._dl*y(h_i)
    temp_m(2,1) = -y(h_i); temp_m(2,2) = -y(a_i)**2*bg_v(y(3:2+nfld)) - k2_in

    ! set 2x2*nfld components
    do i=1,nfld
      cp_constraint(1,i) = y(a_i)**2*bg_dv(y(3:2+nfld),i)
      cp_constraint(1,nfld+i) = y(2+nfld+i)
      cp_constraint(2,i) = y(2+nfld+i)
      cp_constraint(2,nfld+i) = 0._dl
    enddo

    rho_kin = 0._dl
    do i=1,nfld
      rho_kin = rho_kin + 0.5_dl*y(2+nfld+i)**2
    enddo

    ! compute 2x2 * 2x2*nfld
    call dgemm('N','N',2,2*nfld,2,.5_dl,temp_m,2,cp_constraint,2,0._dl, &
               temp_constraint,2)
    !print*, 'temp_constraint(1,2) (no norm) = ', temp_constraint(1,2)
    !print*, 'k2_in = ', k2_in
    cp_constraint = temp_constraint &
                    / (rho_kin - k2_in)
    !print*, 'cp_constraint(1,2)', cp_constraint(1,2)
    !print*, 'cp_constraint(2,1)', cp_constraint(2,1)
                    !/ (y(a_i)**2*bg_v(y(3:2+nfld)) - k2_in + 3._dl*y(h_i)**2)

    !cp_constraint(1,1) = y(a_i)**2*bg_dv(y(f_i)) + y(h_i)*y(df_i)
    !cp_constraint(1,2) = y(df_i)
    !cp_constraint(2,1) = (-k2_in-(y(h_i)**2-0.5_dl*y(df_i)**2))*y(df_i)
    !cp_constraint(1,2) = 0._dl
    ! divide by prefactor
    !cp_constraint = cp_constraint/(-k2_in+0.5_dl*y(df_i)**2)
  end subroutine set_cp_constraint
#endif

  ! Subroutine to calculate the constraint matrix to calculate metric
  ! perturbations from field perturbations in conformal time and longitudinal
  ! gauge. Uses (a Phi), (a Phi)', (a delta phi), (a delta phi)'.
  ! n.b. This is for conformal time
!#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
!  subroutine set_cp_constraint(y, k2_in)
!  !subroutine set_cp_constraint_ascl(y, k2_in)
!    real(dl) :: y(:)  ! array of background variables
!    real(dl) :: k2_in
!
!    integer :: i
!    integer :: a_i, h_i
!    real :: horizon
!
!    a_i = 1; h_i = 2
!    cp_constraint = 0._dl
!    horizon = y(h_i)**2 - diff_bg_hub_conf(y) - k2_in
!    do i=1,nfld
!      cp_constraint(1,i) = -diff_bg_df_conf(y, key=(/i/))
!      cp_constraint(1,nfld+1) = y(2+nfld+i)
!      cp_constraint(2,i) = horizon*y(2+nfld+i)
!    enddo
!    
!     cp_constraint = cp_constraint/(2._dl*horizon)
!  !end subroutine set_cp_constraint_ascl
!  end subroutine set_cp_constraint
!#endif
#if (CPVAROPT == 1 .OR. CPVAROPT == 2 .OR. CPVAROPT == 3)
  subroutine set_cp_constraint(y, k2_in)
    real(dl) :: y(:)  ! array of background variables
    real(dl) :: k2_in

    integer :: i
    integer :: a_i, h_i

    a_i = 1; h_i = 2
    cp_constraint = 0._dl
    do i=1,nfld
#if (CPVAROPT == 1 .OR. CPVAROPT == 2)
      cp_constraint(1,i) = -diff_bg_df_conf(y, key=(/i/))/(y(h_i)**2 - diff_bg_hub_conf(y) - k2_in)
      cp_constraint(1,nfld+1) = y(2+nfld+i)/(y(h_i)**2 - diff_bg_hub_conf(y) - k2_in)
      cp_constraint(2,i) = y(2+nfld+i)
#elif (CPVAROPT == 3)
      cp_constraint(1,i) = (y(h_i)*y(2+nfld+i) + diff_bg_df_cosm(y, key=(/i/)))/(diff_bg_hub_cosm(y) + k2_in/y(a_i)**2)
      cp_constraint(1,nfld+1) = -y(2+nfld+i)/(diff_bg_hub_cosm(y) + k2_in/y(a_i)**2)
      cp_constraint(2,i) = y(2+nfld+i)
#endif
    enddo 
    cp_constraint = cp_constraint/2._dl
  end subroutine set_cp_constraint
#endif

  ! Subroutine to set cp_l_m matrix from an input array of cp system variables
  subroutine set_cp_l_m(cp_sys_in)
    real(dl), dimension(:) :: cp_sys_in

    integer :: i,j,k

    cp_l_m = 0._dl
    k = 1
    do i=1,nfld
      do j=1,i
        !cp_l_m(i,j) = cp_sys_in(4+k)
        cp_l_m(i,j) = cp_sys_in(bg_sys_dim+k)
        k = k+1
      enddo
    enddo
  end subroutine set_cp_l_m
 
  ! Subroutine to set cp_dl_m matrix from an input array of cp system variables
  subroutine set_cp_dl_m(cp_sys_in)
    real(dl), dimension(:) :: cp_sys_in

    integer :: i,j,k,n1

    n1 = nfld*(nfld+1)/2
    cp_dl_m = 0._dl
    k = 1
    do i=1,nfld
      do j=1,i
        !cp_dl_m(i,j) = cp_sys_in(4+n1+k)
        cp_dl_m(i,j) = cp_sys_in(bg_sys_dim+n1+k)
        k = k+1
      enddo
    enddo
  end subroutine set_cp_dl_m

  ! Subroutine to set the cp_c_xdx_m matrix from an input array of cp system
  ! variable 
  subroutine set_cp_c_xdx_m(cp_sys_in)
    real(dl), dimension(:) :: cp_sys_in

    integer :: i,j,k,n1

    n1 = nfld*(nfld+1)/2
    cp_c_xdx_m = 0._dl
    k = 1
    do i=2,nfld
      do j=1,i-1
        !cp_c_xdx_m(i,j) = cp_sys_in(4+2*n1+k)
        cp_c_xdx_m(i,j) = cp_sys_in(bg_sys_dim+2*n1+k)
        k = k+1
      enddo
      do j=i+1,nfld
        cp_c_xdx_m(i,j) = -cp_c_xdx_m(j,i) 
      enddo
    enddo
  end subroutine set_cp_c_xdx_m

  ! Subroutine to set the cp_c_dxdx_m matrix from an input array of cp system
  ! variable
  subroutine set_cp_c_dxdx_m(cp_sys_in)
    real(dl), dimension(:) :: cp_sys_in

    integer :: i,j,k,n1,n2

    n1 = nfld*(nfld+1)/2; n2 = (nfld-1)*nfld/2
    cp_c_dxdx_m = 0._dl
    k = 1
    do i=1,nfld
      do j=1,i
        !cp_c_dxdx_m(i,j) = cp_sys_in(4+2*n1+n2+k)
        cp_c_dxdx_m(i,j) = cp_sys_in(bg_sys_dim+2*n1+n2+k)
        k = k+1
      enddo
      do j=i+1,nfld
        cp_c_dxdx_m(i,j) = cp_c_dxdx_m(j,i)
      enddo
    enddo
  end subroutine set_cp_c_dxdx_m

end module pert_cosmo_mod

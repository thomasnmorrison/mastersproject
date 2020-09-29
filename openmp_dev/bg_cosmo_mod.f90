! bg_cosmo_mod.f90

! Module for a background cosmology

! to do: add an adaptive time step for the evolve subroutine

#DEFINE BGVAROPTCOSM 1
! BGVAROPTCOSM 1 for a, H, phi_A, dot{phi_A}
! BGVAROPTCOSM 2 for ln(a), H, phi_A, dot{phi_A}

module bg_cosmo_mod
!#include "macros.h"

  use params
  use potential_mod
  use gl_integrator_mod 

implicit none

  integer, parameter :: bg_sys_dim = 2+2*nfld

  real(dl), target :: a_bg   ! background scale factor
  real(dl), target :: hub_bg ! background hubble
  real(dl), target, dimension(nfld) :: f_bg   ! background fields
  real(dl), target, dimension(nfld) :: df_bg  ! background field derivatives

  real(dl) :: a_bg_0         ! scale factor at lattice IC
  real(dl) :: hub_bg_0       ! hubble at lattice IC
  real(dl), dimension(nfld) :: f_bg_0         ! fields at lattice IC
  real(dl), dimension(nfld) :: df_bg_0        ! field derivatives at lattice IC

  real(dl):: bg_sys_temp(bg_sys_dim)            ! temp copy bg variables
  real(dl), allocatable :: bg_sys_g(:,:,:)       ! gl g for bg system
  type(ode_sys), dimension(bg_sys_dim) :: bg_ode ! ode system of background
  !type(ode_sys), dimension(bg_sys_dim) :: bg_ode_cosm

contains

  ! Subroutine to initialize background cosmology equations and memory
  ! allocation
  subroutine init_bg(da, dh, df, ddf)
    procedure(d_template) :: da, dh, df, ddf    

    integer :: i    

    if(.NOT. allocated(bg_sys_g)) then
      call gl_init_g(bg_sys_g, bg_sys_dim)
    endif

    ! Set pointers to derivative functions
    bg_ode(1)%f => da  !diff_bg_a_conf
    bg_ode(2)%f => dh  !diff_bg_hub_conf
    do i=1,nfld
      bg_ode(2+i)%f => df  !diff_bg_f_conf
      bg_ode(2+nfld+i)%f => ddf !diff_bg_df_conf
    enddo

    ! Set pointers to system variables
    bg_ode(1)%y => a_bg
    bg_ode(2)%y => hub_bg
    do i=1,nfld
      bg_ode(2+i)%y => f_bg(i)
      bg_ode(2+nfld+i)%y => df_bg(i)
    enddo

    ! Set keys
    do i=2+1,bg_sys_dim
      if(.NOT. allocated(bg_ode(i)%key)) then
        allocate(bg_ode(i)%key(1))
      endif
    enddo
    do i=1,nfld
      bg_ode(2+i)%key = (/i/)
      bg_ode(2+nfld+i)%key = (/i/)
    enddo
    
    do i=1,bg_sys_dim
      bg_ode(i)%ode_param => null()
    enddo
  end subroutine init_bg

  ! Subroutine to set cp_sys_temp from cp_ode
  subroutine set_bg_sys_temp()
    integer :: i
    do i=1, bg_sys_dim
      bg_sys_temp(i) = bg_ode(i)%y
    enddo
  end subroutine set_bg_sys_temp

  ! Subroutine to set conditions in the background cosmology from supplied
  ! conditions
  subroutine set_bg_cosmo_ic(a_in, hub_in, f_in, df_in)
    real(dl), intent(in) :: a_in
    real(dl), intent(in) :: hub_in
    real(dl), intent(in) :: f_in(:)
    real(dl), intent(in) :: df_in(:)

    a_bg = a_in
    hub_bg = hub_in
    f_bg = f_in
    df_bg = df_in
  end subroutine set_bg_cosmo_ic

  ! Function to calculate the conformal time derivative of background
  ! inflaton
  function diff_bg_f_conf(y, ode_param, key)
    real(dl) :: diff_bg_f_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: df_i ! index locations in y of bg variables

    df_i = 2+nfld+key(1)

    diff_bg_f_conf = y(df_i)
  end function diff_bg_f_conf

  ! Function to calculate the cosmic time derivative of the inflaton
  function diff_bg_f_cosm(y, ode_param, key)
    real(dl) :: diff_bg_f_cosm
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: df_i

    df_i = 2+nfld+key(1)

    diff_bg_f_cosm = y(df_i)
  end function diff_bg_f_cosm

  ! Function to calculate the 2nd conformal time derivative of the inflaton
  function diff_bg_df_conf(y, ode_param, key)
    real(dl) :: diff_bg_df_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: a_i, h_i, f_i, df_i

    a_i = 1; h_i =2; f_i = 2+key(1); df_i = 2+nfld+key(1)
    diff_bg_df_conf = -2._dl*y(h_i)*y(df_i)-y(a_i)**2*bg_dv(y(3:2+nfld),key(1)) 
  end function diff_bg_df_conf

  ! Function to calculate the cosmic time 2nd  derivative of the inflaton
  function diff_bg_df_cosm(y, ode_param, key)
    real(dl) :: diff_bg_df_cosm
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: hub_i, f_i, df_i

    hub_i = 2; f_i = 3; df_i = 2+nfld+key(1)

    !diff_bg_df_cosm = -3._dl*y(hub_i)*y(df_i) - bg_dv(y(f_i))
    diff_bg_df_cosm = -3._dl*y(hub_i)*y(df_i) - bg_dv(y(3:2+nfld),key(1))
  end function diff_bg_df_cosm

  ! Function to calculate the conformal time derivative of the scale
  ! factor
  function diff_bg_a_conf(y, ode_param, key)
    real(dl) :: diff_bg_a_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: a_i, hub_i

    a_i = 1; hub_i = 2
    diff_bg_a_conf = y(a_i)*y(hub_i)
  end function diff_bg_a_conf

  ! Function to calculate the cosmic time derivative of the scale factor
  function diff_bg_a_cosm(y, ode_param, key)
    real(dl) :: diff_bg_a_cosm
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: a_i, hub_i

    a_i = 1; hub_i = 2
#if (BGVAROPTCOSM == 1)
    diff_bg_a_cosm = y(a_i)*y(hub_i)
#elif (BGVAROPTCOSM == 2)
    diff_bg_a_cosm = y(hub_i)
#endif
  end function diff_bg_a_cosm

  ! Function to calulate the conformal time derivative of the (conformal)
  ! Hubble parameter
  function diff_bg_hub_conf(y, ode_param, key)
    real(dl) :: diff_bg_hub_conf
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i
    integer :: a_i, hub_i, df_i

    a_i = 1; hub_i = 2
    !diff_bg_hub_conf = y(hub_i)**2 - 0.5_dl*y(df_i)**2/y(a_i)**4
    !diff_bg_hub_conf = y(hub_i)**2 - 0.5_dl*y(df_i)**2
    diff_bg_hub_conf = 0._dl
    do i=1,nfld
      df_i = 2+nfld+i
      diff_bg_hub_conf = diff_bg_hub_conf + y(df_i)**2
    enddo
    diff_bg_hub_conf = y(hub_i)**2 - 0.5_dl*diff_bg_hub_conf
  end function diff_bg_hub_conf

  ! Function to calculate the cosmic time derivative of the (cosmic) Hubble
  ! parameter
  function diff_bg_hub_cosm(y, ode_param, key)
    real(dl) :: diff_bg_hub_cosm
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    integer :: i
    integer :: df_i

    !diff_bg_hub_cosm = -0.5_dl*y(df_i)**2
    diff_bg_hub_cosm = 0._dl
    do i=1,nfld
      df_i = 2+nfld+i
      diff_bg_hub_cosm = diff_bg_hub_cosm + y(df_i)**2
    enddo
    diff_bg_hub_cosm = -0.5_dl*diff_bg_hub_cosm
  end function diff_bg_hub_cosm

  ! Function to return a zero derivative
  function diff_bg_zero(y, ode_param, key)
    real(dl) :: diff_bg_zero
    real(dl) :: y(:)
    real(dl), pointer, optional :: ode_param(:)
    integer, optional :: key(:)

    diff_bg_zero = 0._dl
  end function diff_bg_zero

  ! Subroutine to step the evolution of the background cosmology
  subroutine step_bg_cosmo(dt)
    real(dl) :: dt  ! cosmic time step

    call gl_integrate(bg_sys_temp, bg_ode, bg_sys_g, dt, bg_sys_dim)
  end subroutine step_bg_cosmo

  ! Subroutine evolve the background cosmology to a specified state
  ! n.b. this subroutine assumes ind corresponds to a monotonic variable
  ! n.b. if ind2 is supplied then the product of the variables corresponding
  !      to ind1 and ind2 is compared to y_end
  ! to do: the ind2 version is going to require another version of the interpolator
  subroutine evolve_bg_cosmo(y_end, ind1, dt, ind2)
    real(dl), intent(in) :: y_end  ! variable y at end point of evolution
    integer, intent(in) :: ind1     ! index of variable y
    real(dl) :: dt
    integer, optional, intent(in) :: ind2

    integer :: i
    real(dl) :: diff, diff1, diff2

    ! check sign of derivative
    do i=1,bg_sys_dim
      bg_sys_temp(i) = bg_ode(i)%y
    enddo
    diff1 = bg_ode(ind1)%f(bg_sys_temp, ode_param=bg_ode(ind1)%ode_param, key=bg_ode(ind1)%key)
    if (present(ind2)) then
      diff2 = bg_ode(ind2)%f(bg_sys_temp, ode_param=bg_ode(ind2)%ode_param, key=bg_ode(ind2)%key)
      ! apply product rule
      diff = bg_ode(ind1)%y*diff2 + bg_ode(ind2)%y*diff1
    else
      diff = diff1
    endif


    if (sign(1._dl,diff)*(bg_ode(ind1)%y-y_end) .LE. 0._dl) then
      do while(sign(1._dl,diff)*(bg_ode(ind1)%y-y_end) .LE. 0._dl)
        call step_bg_cosmo(dt)
      enddo
    else
      do while(sign(1._dl,diff)*(bg_ode(ind1)%y-y_end) .GE. 0._dl)
        call step_bg_cosmo(-dt)
      enddo
    endif
    ! put intoerpolation here
    call gl_newton_root(bg_sys_temp, bg_ode, bg_sys_g, dt, bg_sys_dim, y_end,ind1)
  end subroutine evolve_bg_cosmo


  ! Subroutine to get the homogenious initial conditions on the lattice. Starts
  ! from an initial guess in the distant pass and allows to settle onto the
  ! attractor.
  ! n.b. this subroutine is done in cosmic time
  subroutine get_bg_attractor(y_end, ind, dt, a_in, hub_in, f_in, df_in)
    real(dl) :: y_end  ! inflaton at lattice IC
    integer :: ind
    real(dl) :: dt     ! time step
    real(dl) :: a_in
    real(dl) :: hub_in
    real(dl) :: f_in(:)
    real(dl) :: df_in(:)

    call set_bg_cosmo_ic(a_in, hub_in, f_in, df_in)
    call evolve_bg_cosmo(y_end, ind, dt)

    a_bg_0 = 1.0_dl
    hub_bg_0 = hub_bg
    f_bg_0 = f_bg
    df_bg_0 = df_bg
  end subroutine get_bg_attractor

  ! Subroutine to convert cosmic variables to conformal variables
  subroutine cosm_to_conf()
    hub_bg = a_bg*hub_bg
    !df_bg = a_bg**3*df_bg
    df_bg = a_bg*df_bg
  end subroutine cosm_to_conf

  ! Subroutine to convert conformal variables to cosmic variables
  subroutine conf_to_cosm()
    hub_bg = hub_bg/a_bg
    !df_bg = df_bg/a_bg**3
    df_bg = df_bg/a_bg
  end subroutine conf_to_cosm

  ! Function to calculate the background equation of state from conformal
  ! time variables
  function get_bg_w_conf() result(w)
    real(dl) :: w

    integer :: i
    !w = (0.5_dl*df_bg**2/a_bg**2 - bg_v(f_bg)) &
    !    /(0.5_dl*df_bg**2/a_bg**2 + bg_v(f_bg))
    w = 0._dl
    do i=1,nfld
      w = w + 0.5_dl*df_bg(i)**2
    enddo
    w = (w/a_bg**2 - bg_v(f_bg))/(w/a_bg**2 + bg_v(f_bg))
  end function get_bg_w_conf

  ! Function to calculate the background equation of state from cosmic
  ! time variables
  function get_bg_w_cosm() result(w)
    real(dl) :: w

    integer :: i
    !w = (0.5_dl*df_bg**2 - bg_v(f_bg)) &
    !    /(0.5_dl*df_bg**2 + bg_v(f_bg))
    w = 0._dl
    do i=1,nfld
      w = w + 0.5_dl*df_bg(i)**2
    enddo
    w = (w - bg_v(f_bg))/(w + bg_v(f_bg))
  end function get_bg_w_cosm

  ! Function to calculate conformal hubble
  ! n.b. This function is need to set the initial conditions for the hubble
  !      parameter, bescause the Friedmann equation is a constraint equation.
  function get_bg_hub_conf(a_in, f_in, df_in)
    real(dl) :: get_bg_hub_conf
    real(dl), optional :: a_in
    real(dl), optional :: f_in(:)  ! optional input for inflaton
    real(dl), optional :: df_in(:) ! optional input for inflaton (conf) momentum

    integer :: i

    get_bg_hub_conf = 0._dl
    if (present(a_in) .and. present(f_in) .and. present(df_in)) then
      do i=1,nfld
        get_bg_hub_conf = get_bg_hub_conf + 0.5_dl*df_in(i)**2
      enddo
      get_bg_hub_conf = sqrt((get_bg_hub_conf &
                              + a_in**2*bg_v(f_in))/3._dl)
      !get_bg_hub_conf = sqrt((0.5_dl*df_in**2 &
      !                        + a_in**2*bg_v(f_in))/3._dl)
    else
      do i=1,nfld
        get_bg_hub_conf = get_bg_hub_conf + 0.5_dl*df_bg(i)**2
      enddo
      get_bg_hub_conf = sqrt((get_bg_hub_conf &
                              + a_bg**2*bg_v(f_bg))/3._dl)
      !get_bg_hub_conf = sqrt((0.5_dl*df_bg**2 &
      !                        + a_bg**2*bg_v(f_bg))/3._dl)
    endif
  end function get_bg_hub_conf

  ! Function to calculate cosmic hubble
  ! n.b. This function is need to set the initial conditions for the hubble
  !      parameter, bescause the Friedmann equation is a constraint equation. 
  function get_bg_hub_cosm(f_in, df_in)
    real(dl) :: get_bg_hub_cosm
    real(dl), optional :: f_in(:)  ! optional input for inflaton
    real(dl), optional :: df_in(:) ! optional input for inflaton (cosm) momentum

    integer :: i
    
    get_bg_hub_cosm = 0._dl
    if (present(f_in) .and. present(df_in)) then
      do i=1,nfld
        get_bg_hub_cosm = get_bg_hub_cosm + 0.5_dl*df_in(i)**2
      enddo
      get_bg_hub_cosm = sqrt((get_bg_hub_cosm + bg_v(f_in))/3._dl)
      !get_bg_hub_cosm = sqrt((0.5_dl*df_in**2 + bg_v(f_in))/3._dl)
    else
      do i=1,nfld
        get_bg_hub_cosm = get_bg_hub_cosm + 0.5_dl*df_bg(i)**2
      enddo
      get_bg_hub_cosm = sqrt((get_bg_hub_cosm + bg_v(f_bg))/3._dl)
      !get_bg_hub_cosm = sqrt((0.5_dl*df_bg**2 + bg_v(f_bg))/3._dl)
    endif
  end function get_bg_hub_cosm

  function get_bg_a(y)
    real(dl) :: get_bg_a
    real(dl) :: y(:)
#if (BGVAROPTCOSM == 1)
    get_bg_a = y(1)
#elif (BGVAROPTCOSM == 2)
    get_bg_a = exp(y(1))
#endif
  end function get_bg_a

  function get_bg_hub(y)
    real(dl) :: get_bg_hub
    real(dl) :: y(:)
    get_bg_hub = y(2)
  end function get_bg_hub

  function get_bg_f(y, ind)
    real(dl) :: get_bg_f
    real(dl) :: y(:)
    integer :: ind
    get_bg_f = y(2+ind)
  end function get_bg_f

  function get_bg_df(y, ind)
    real(dl) :: get_bg_df
    real(dl) :: y(:)
    integer :: ind
    get_bg_df = y(2+nfld+ind)
  end function get_bg_df

end module bg_cosmo_mod

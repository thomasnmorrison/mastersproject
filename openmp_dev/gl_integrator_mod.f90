! gl_integrator_mod.f90

! Module for a Gauss-Legendre integrator
! Reference J. C. Butcher, Implicit Runge-Kutta Processes for naming
! convensions.

! to do: write parameters for Butcher's tables of higher order integrators
! to do: test 8th order integrator

#define CONVTEST 0

module gl_integrator_mod

  use params

implicit none

  abstract interface
    !import
    ! Derivative function template
    function d_template(y, ode_param, key) result(f)
      import
      !real(dl), pointer :: y(:)
      real(dl) :: y(:)
      !integer, intent(in) :: dim
      real(dl), optional, pointer :: ode_param(:)
      !integer, optional, pointer :: key(:)
      integer, optional :: key(:)
      real(dl) :: f
    end function d_template
  end interface

  ! Derived data type with the derivative data for an ODE 
  type :: ode_sys
    procedure(d_template), pointer, nopass :: f
    real(dl), pointer :: y
    !integer, pointer :: key(:)
    integer, allocatable :: key(:)
    real(dl), pointer :: ode_param(:)
  end type ode_sys

  ! Butcher's table and parameters fo order 2
  integer, parameter :: nu_gl2 = 1
  real(dl), parameter :: a11_gl2 = 0.5_dl
  real(dl), parameter, dimension(1) :: a_gl2_flat = (/a11_gl2/)
  real(dl), parameter, dimension(1,1) :: a_gl2 = reshape(a_gl2_flat,(/1,1/))
  real(dl), parameter :: b1_gl2 = 1._dl
  real(dl), parameter, dimension(nu_gl2) :: b_gl2 = (/b1_gl2/)

  ! Butcher's table order 4
  integer, parameter :: nu_gl4 = 2
  real(dl), parameter :: a11_gl4 = 0.25_dl, a12_gl4 = -0.0386751345948128822_dl
  real(dl), parameter :: a21_gl4 = 0.53867513459481288_dl, a22_gl4 = 0.25_dl
  real(dl), parameter, dimension(2,2) :: a_gl4 = reshape((/a11_gl4, a21_gl4, a12_gl4, a22_gl4/),(/2,2/))
  real(dl), parameter :: b1_gl4 = 0.5_dl, b2_gl4 = 0.5_dl
  real(dl), parameter, dimension(nu_gl4) :: b_gl4 = (/b1_gl4, b2_gl4/)

  ! Butcher's table order 6
  integer, parameter :: nu_gl6 = 3
  real(dl), parameter :: a11_gl6 = 0.138888888888888_dl, &
                         a12_gl6 = -0.0359766675249389_dl, &
                         a13_gl6 = 0.00978944401530832_dl
  real(dl), parameter :: a21_gl6 = 0.300263194980964_dl, &
                         a22_gl6 = 0.222222222222222_dl, &
                         a23_gl6 = -0.022485417203086_dl
  real(dl), parameter :: a31_gl6 = 0.267988333762469_dl, &
                         a32_gl6 = 0.480421111969383_dl, &
                         a33_gl6 = 0.138888888888888_dl
  real(dl), parameter, dimension(3,3) :: a_gl6 = reshape((/a11_gl6, a21_gl6, a31_gl6, a12_gl6, a22_gl6, a32_gl6, a13_gl6, a23_gl6, a33_gl6/),(/3,3/))
  real(dl), parameter :: b1_gl6 = 0.277777777777777_dl, &
                         b2_gl6 = 0.444444444444444_dl, &
                         b3_gl6 = 0.277777777777777_dl
  real(dl), parameter, dimension(nu_gl6) :: b_gl6 = (/b1_gl6, b2_gl6, b3_gl6/)

  ! Butcher's table order 8
  integer, parameter :: nu_gl8 = 4
  real(dl), parameter :: a11_gl8 = 0.08696371128436346434_dl,&
                         a12_gl8 = -0.02660418008499879331_dl,&
                         a13_gl8 = 0.01262746268940472451_dl,&
                         a14_gl8 = -0.00355514968579568317_dl
  real(dl), parameter :: a21_gl8 = 0.18811811749986807163_dl,&
                         a22_gl8 = 0.16303628871563653565_dl,&
                         a23_gl8 = -0.02788042860247089521_dl,&
                         a24_gl8 = 0.00673550059453815551_dl
  real(dl), parameter :: a31_gl8 = 0.16719192197418877317_dl,&
                         a32_gl8 = 0.35395300603374396651_dl,&
                         a33_gl8 = 0.16303628871563653565_dl,&
                         a34_gl8 = -0.01419069493114114295_dl
  real(dl), parameter :: a41_gl8 = 0.17748257225452261185_dl,&
                         a42_gl8 = 0.31344511474186834679_dl,&
                         a43_gl8 = 0.35267675751627186461_dl,&
                         a44_gl8 = 0.08696371128436346434_dl
  real(dl), parameter, dimension(4,4) :: a_gl8 = reshape((/a11_gl8, a21_gl8, a31_gl8, a41_gl8, &
                                                           a12_gl8, a22_gl8, a32_gl8, a42_gl8, &
                                                           a13_gl8, a23_gl8, a33_gl8, a43_gl8, &
                                                           a14_gl8, a24_gl8, a34_gl8, a44_gl8/),(/4,4/))
  real(dl), parameter :: b1_gl8 = 0.17392742256872692868_dl,&
                         b2_gl8 = 0.32607257743127307130_dl,&
                         b3_gl8 = 0.32607257743127307130_dl,&
                         b4_gl8 = 0.17392742256872692868_dl
  real(dl), parameter, dimension(nu_gl8) :: b_gl8 = (/b1_gl8, b2_gl8, b3_gl8, b4_gl8/)

  ! interfacing variables
  integer ::  nu_gl
  real(dl), allocatable :: b_gl(:)
  !real(dl), allocatable :: c_gl(:,:)
  real(dl), allocatable :: a_gl(:,:)

  ! Convergence parameters
  real(dl), parameter :: eps_0 = 1._dl/dble(2**40)!1._dl/dble(2**32)   ! default convergence criteria
  integer, parameter :: c_max_0 = 2**4   ! default maximum number of iterations solving g

contains

  ! Subroutine to perform integration step with 2nd order GL integrator
  subroutine gl2_integrate(y, z, ode, g, h, dim, eps, c_max)
    real(dl), pointer :: y(:)    ! system configuration
    real(dl), pointer :: z(:)    ! intermediate system configuration copy
    type(ode_sys) :: ode(:)      ! ode for system
    real(dl), pointer :: g(:,:,:)
    real(dl), intent(in) :: h    ! step size
    real(dl), intent(in) :: eps  ! convergence criteria epsilon
    integer, intent(in) :: dim   ! dimensionality of system
    integer, intent(in) :: c_max

    integer :: i,k

    call gl2_solve_g(y, z, ode, g, h, eps, dim, c_max)
    
    z(:) = 0._dl
    do i=1, nu_gl2!a_gl2
      z(:) = z(:) + b_gl2(i)*g(2,:,i)
    enddo

    y = y + h*z
  end subroutine gl2_integrate

  ! Subroutine to perform one integration step of an ODE system using a GL
  ! scheme as initialized by gl_init_scheme
  subroutine gl_integrate(z, ode, g, h, dim, eps_in, c_max_in)
    !real(dl), pointer :: y(:)    ! system configuration
    !real(dl), pointer :: z(:)    ! intermediate system configuration copy
    real(dl) :: z(:)    ! intermediate system configuration copy
    type(ode_sys) :: ode(:)      ! ode for system
    !real(dl), pointer :: g(:,:,:)
    real(dl) :: g(:,:,:)
    real(dl), intent(in) :: h    ! step size
    real(dl), optional :: eps_in  ! convergence criteria epsilon
    integer,  intent(in) :: dim   ! dimensionality of system
    integer, optional :: c_max_in ! loop cut-off

    integer :: i,k
    real(dl) :: eps
    integer :: c_max

    if (.not. present(eps_in)) then
      eps = eps_0
    else
      eps = eps_in
    endif
    if (.not. present(c_max_in)) then
      c_max = c_max_0
    else
      c_max = c_max_in
    endif

    call gl_solve_g(z, ode, g, h, eps, dim, c_max)
    ! print for testing
    !print*, 'g = ', g(2,:,:)    

    z(:) = 0._dl
    do i=1, nu_gl
      z(:) = z(:) + b_gl(i)*g(2,:,i)
    enddo
    !y = y + h*z
    do i=1, dim
      ode(i)%y = ode(i)%y + h*z(i)
    enddo
  end subroutine gl_integrate  

  ! Subroutine to solve for the g vector for the input ode system using a 
  ! 2nd order scheme
  subroutine gl2_solve_g(y, z, ode, g, h, eps, dim, c_max)
    real(dl), pointer :: y(:)
    real(dl), pointer :: z(:)
    type(ode_sys) :: ode(:)
    real(dl) :: g(:,:,:)    ! (old/new,dim,i)
    real(dl), intent(in) :: h
    real(dl), intent(in) :: eps
    integer, intent(in) :: dim, c_max
 
    integer :: k, c

    do k=1, dim
      c = 0
      call gl2_recurse_g(y, z, ode, g, h, eps, dim, k, c_max, c)
    enddo
  end subroutine gl2_solve_g

  ! Subroutine to solve for the g vector for the input ode system using the 
  ! scheme initialized in gl_init_scheme
  subroutine gl_solve_g(z, ode, g, h, eps, dim, c_max)
    !real(dl), pointer :: y(:)
    !real(dl), pointer :: z(:)
    real(dl) :: z(:)
    type(ode_sys) :: ode(:)
    real(dl) :: g(:,:,:)    ! (old/new,dim,i)
    real(dl), intent(in) :: h
    real(dl), intent(in) :: eps
    integer, intent(in) :: dim, c_max

    integer :: k, c

    !do k=1, dim
    !  c = 0
    !  call gl_recurse_g_old(z, ode, g, h, eps, dim, k, c_max, c)
    !  ! print for testing
    !  print*, 'k, g = ', k, g(2,:,:)
    !enddo
    c = 0
    g(:,:,:) = 0._dl
    call  gl_recurse_g(z, ode, g, h, eps, dim, c_max, c)
    !print*, g(:,:,:)
    !print*, 'c = ', c
    !print*, 'c = ', c
    !print*, 'g(2,1,:) = ', g(2,1,:)
    !print*, 'g(2,2,:) = ', g(2,2,:)
  end subroutine gl_solve_g

  ! Subroutine to iterate an estimate of g recursively to convergence for a
  ! 2nd order scheme
  recursive subroutine gl2_recurse_g(y, z, ode, g, h, eps, dim, k, c_max, c)
    real(dl), pointer :: y(:), z(:)
    type(ode_sys) :: ode(:)
    real(dl) :: g(:,:,:)
    real(dl), intent(in) :: h
    real(dl), intent(in) :: eps
    integer, intent(in) :: dim
    integer, intent(in) :: k, c_max
    integer :: c
    integer :: i,j

    c = c + 1            ! count recursions
    g(1,k,:) = g(2,k,:)  ! update values from previous recursion

    ! Iterate guess of g
    do i=1, nu_gl2
      z = y
      do j=1, i-1
        z = z + h*a_gl2(i,j)*g(2,k,j)
      enddo
      do j=i, nu_gl2
        z = z + h*a_gl2(i,j)*g(1,k,j)
      enddo
      g(2,k,i) = ode(k)%f(z, ode(k)%ode_param, ode(k)%key)
    enddo
    
    ! Exception for convergence failure
    if (c .GE. c_max) then
      print*, 'CONVERGENCE FAILURE'
      return
    endif

    ! Check convergence and recurse
    !if (max(abs(g(2,k,:) - g(1,k,:))) .GT. eps) then
    if (maxval(abs(g(2,k,:) - g(1,k,:))) .GT. eps) then
      call gl2_recurse_g(y, z, ode, g, h, eps, dim, k, c, c_max)
    endif
  end subroutine gl2_recurse_g

  ! Subroutine to iterate an estimate of g recursively to convergence for the
  ! scheme initialized by gl_init_scheme
  recursive subroutine gl_recurse_g_old(z, ode, g, h, eps, dim, k, c_max, c)
    real(dl), pointer :: z(:)        ! temporary system configuration
    type(ode_sys) :: ode(:)          ! ode system
    real(dl) :: g(:,:,:)             !
    real(dl), intent(in) :: h        ! time step
    real(dl), intent(in) :: eps      ! congergence criteria
    integer, intent(in) :: dim       ! dimension of system
    integer, intent(in) :: k         ! component being solved
    integer, intent(in) :: c_max    ! recursion cut-off
    integer :: c
    integer :: i,j

    !real(dl) :: test_var  ! only for testing

    !c = c + 1            ! count recursions
    !g(1,k,:) = g(2,k,:)  ! update values from previous recursion

    ! Iterate guess of g
    do i=1, nu_gl
      do j=1, dim
        z(j) = ode(j)%y
      enddo
      do j=1, i-1
        ! bug here z is the intermediate point at which the derivatives are
        ! evaluated, it in principle depends on all components of g. The bug
        ! here is to use only the k^th component.
        z = z + h*a_gl(i,j)*g(2,k,j)
      enddo
      do j=i, nu_gl
        z = z + h*a_gl(i,j)*g(1,k,j)
      enddo
      if (associated(ode(k)%ode_param)) then
        g(2,k,i) = ode(k)%f(z, ode(k)%ode_param, ode(k)%key)
      else
        g(2,k,i) = ode(k)%f(z, key=ode(k)%key)
      endif
    enddo

    ! Exception for convergence failure
    !if (c .GE. c_max) then
    !  print*, 'CONVERGENCE FAILURE'
    !  print*, 'ERROR = ', maxval(abs(g(2,k,:) - g(1,k,:)))
    !  print*, 'eps = ', eps
    !  print*, 'c = ', c
    !  return
    !endif

    ! Check convergence and recurse
    if (maxval(abs(g(2,k,:) - g(1,k,:))) .GT. eps) then
      !print*, 'check 10d'
      call gl_recurse_g_old(z, ode, g, h, eps, dim, k, c, c_max)
    endif
    if (maxval(abs(g(2,k,:) - g(1,k,:))) .LT. eps) then
      !print*, 'CONVERGED'
      return
    endif
    if (c .GE. c_max) then
      print*, 'CONVERGENCE FAILURE'
      print*, 'ERROR = ', maxval(abs(g(2,k,:) - g(1,k,:)))
      print*, 'eps = ', eps
      print*, 'c = ', c
      return
    endif
  end subroutine gl_recurse_g_old

  ! subroutine to recurse g
  ! n.b. the if statements at the end will slow this down. Once it is determined
  !      how many recursions are needed it will be faster to just run
  !      an upper bound on that every time.
  recursive subroutine gl_recurse_g(z, ode, g, h, eps, dim, c_max, c)
    !real(dl), pointer :: z(:)        ! temporary system configuration
    real(dl):: z(:)
    type(ode_sys) :: ode(:)          ! ode system
    real(dl) :: g(:,:,:)             !
    real(dl), intent(in) :: h        ! time step
    real(dl), intent(in) :: eps      ! congergence criteria
    integer, intent(in) :: dim       ! dimension of system
    integer, intent(in) :: c_max    ! recursion cut-off
    integer :: c
    integer :: i,j,k

    ! count recursions
    c = c + 1
    !print*, 'c = ', c
    if (c .GE. c_max) then
      print*, 'CONVERGENCE FAILURE'
      print*, 'ERROR = ', maxval(abs(g(2,:,:) - g(1,:,:)))
      !print*, 'g(2,:,:) = ', g(2,:,:)
      !print*, 'g(1,:,:) = ', g(1,:,:)
      print*, 'eps = ', eps
      print*, 'c = ', c
      return
    endif
    ! update iteration of g
#if (CONVTEST == 1)
    print*, 'ERROR = ', maxval(abs(g(2,:,:) - g(1,:,:)))
#endif
    g(1,:,:) = g(2,:,:)
    ! calculate next iteration of g
    !do i=1,nu_gl
    do i=1,nu_gl
      !print*, 'i = ', i
      ! calculate intermidiate point
      ! i=2 calculation of z is incorrect
      do k=1, dim
        z(k) = ode(k)%y
        do j=1, nu_gl
          !print*, 'k,z,g=', k, z(k), g(1,k,j)
          z(k) = z(k) + h*a_gl(i,j)*g(1,k,j)
        enddo
      enddo
      !print*, 'i = ', i
      !print*, 'z(:) = ', z(:)
      ! n.b. these two loops over k are performed in sequence, so cannot be
      ! combined into one loop
      do k=1, dim
        if (associated(ode(k)%ode_param)) then
          g(2,k,i) = ode(k)%f(z, ode(k)%ode_param, ode(k)%key)
        else
          g(2,k,i) = ode(k)%f(z, key=ode(k)%key)
        endif
      enddo
      !print*, 'i = ', i
      !print*, 'z(:) = ', z(:)
    enddo

    !print*, 'c = ', c
    !print*, 'g1 = ', g(1,:,:)
    !print*, 'g(2,1,:) = ', g(2,1,:)
    !print*, 'g(2,2,:) = ', g(2,2,:)

    ! Check convergence and recurse
    if (maxval(abs(g(2,:,:) - g(1,:,:))) .GT. eps) then
      call gl_recurse_g(z, ode, g, h, eps, dim, c_max, c)
    endif
    if (maxval(abs(g(2,:,:) - g(1,:,:))) .LT. eps) then
      !print*, 'CONVERGED'
#if (CONVTEST == 1)
      print*, 'ERROR, c = ', maxval(abs(g(2,:,:) - g(1,:,:))), c
#endif
      !print*
      return
    endif
    !if (c .GE. c_max) then
    !  print*, 'CONVERGENCE FAILURE'
    !  print*, 'ERROR = ', maxval(abs(g(2,k,:) - g(1,k,:)))
    !  print*, 'eps = ', eps
    !  print*, 'c = ', c
    !  return
    !endif
  end subroutine gl_recurse_g

  ! Subroutine to impliment Newton's root finding method using the gl integrator
  ! for function evaluations
  ! n.b. in the first instance I'm making this non recursive, just using a loop
  ! n.b. this subroutine does not find the time of the root, only the
  !      configuration of the system at the root
  subroutine gl_newton_root(z, ode, g, h, dim, rt, ind, eps_in, c_max_in)
    !real(dl), pointer :: z(:)    ! intermediate system configuration copy
    real(dl) :: z(:)
    type(ode_sys) :: ode(:)      ! ode for system
    real(dl) :: g(:,:,:)
    real(dl), intent(in) :: h    ! step size
    real(dl), optional :: eps_in  ! convergence criteria epsilon
    integer,  intent(in) :: dim   ! dimensionality of system
    real(dl), intent(in) :: rt    ! value at which to converge z(ind)
    integer, intent(in) :: ind    ! index of root variable
    integer, optional :: c_max_in ! loop cut-off

    real(dl) :: eps
    integer :: c, c_max
    real(dl) :: h_tmp  ! updating time step
    integer :: k

    if (.not. present(eps_in)) then; eps = eps_0; else; eps = eps_in; endif
    if (.not. present(c_max_in)) then; c_max = c_max_0; else; c_max = c_max_in; endif
    
    ! iterate Newton's method
    do c=1, c_max
      do k=1, dim
        z(k) = ode(k)%y
      enddo
      h_tmp = -(ode(ind)%y - rt)/ode(ind)%f(z,ode(ind)%ode_param,ode(ind)%key)
      call gl_integrate(z, ode, g, h_tmp, dim)
    enddo
    
    ! check convergence
    if (abs(ode(ind)%y - rt) .GT. eps) then
      print*, 'CONVERGENCE FAILURE'
      print*, 'ERROR = ', abs(ode(ind)%y - rt)
      print*, 'eps = ', eps
      print*, 'c_max = ', c_max
    endif
  end subroutine gl_newton_root

  ! Subroutine to allocate memory to a new/old GL g
  subroutine gl2_init_g(g, dim_in)
    real(dl), allocatable :: g(:,:,:)
    integer :: dim_in

    allocate(g(2,dim_in,1))
  end subroutine gl2_init_g

  ! Subroutine to allocate memory to new/old GL g
  ! This subroutine is to simplify changing the integrator order
  subroutine gl_init_g(g, dim_in, g_0)
    real(dl), allocatable :: g(:,:,:)
    integer :: dim_in
    real(dl), optional :: g_0  ! initial value of g

    if(.NOT. allocated(g)) then
     allocate(g(2,dim_in,nu_gl))
    endif
    if (present(g_0)) then
      g = g_0
    endif
  end subroutine gl_init_g

  ! Subroutine to initialize integration scheme constants
  ! n.b. call this subroutine first
  subroutine gl_init_scheme(order)
    integer, intent(in) :: order

    select case(order)
      case(2)
        nu_gl = nu_gl2
        allocate(b_gl(nu_gl))
        b_gl = b_gl2
        allocate(a_gl(nu_gl,nu_gl))
        a_gl = a_gl2
      case(4)
        nu_gl = nu_gl4
        allocate(b_gl(nu_gl))
        b_gl = b_gl4
        allocate(a_gl(nu_gl,nu_gl))
        a_gl = a_gl4
        print*, 'nu_gl = ', nu_gl
        print*, 'b_gl(1) = ', b_gl(1)
        print*, 'b_gl(2) = ', b_gl(2)
        print*, 'a_gl(1,1) = ', a_gl(1,1)
        print*, 'a_gl(1,2) = ', a_gl(1,2)
        print*, 'a_gl(2,1) = ', a_gl(2,1)
        print*, 'a_gl(2,2) = ', a_gl(2,2)
      case(6)
        nu_gl = nu_gl6
        allocate(b_gl(nu_gl))
        b_gl = b_gl6
        allocate(a_gl(nu_gl,nu_gl))
        a_gl = a_gl6
      case(8)
        nu_gl = nu_gl8
        allocate(b_gl(nu_gl))
        b_gl = b_gl8
        allocate(a_gl(nu_gl,nu_gl))
        a_gl = a_gl8
      case default
        print*, "COULD NOT INITIALIZE SCHEME"
    end select
  end subroutine gl_init_scheme

  ! Subroutine to deallocate integration scheme constants
  subroutine gl_del_scheme()
    deallocate(b_gl)
    deallocate(a_gl)
  end subroutine gl_del_scheme

end module

! potential_mod_new.f90

! Module for functions and parameters relating to the potential.

! The intent is to rewrite potential_mod in a cleaner way, removing redundancy
! and allowing the potential to be selected by preprocessor statements and
! macros.

! to do: write the functions for a constant potential


#define POTOPT 0
! POTOPT 0 for constant potential
! POTOPT 1 for constant chi mass
! POTOPT 2 for localized instability
! POTOPT 3 for instability turn on
!#define INFPWR 1

module potential_mod
#include "macros.h"

  use params

implicit none

#ifdef ONEFLD
  integer, parameter :: nfld = 1
#endif
#ifdef TWOFLD
  integer, parameter :: nfld = 2
#endif

  real(dl), parameter :: mpl = 1.e5 ! Planck mass to mass unit ratio

  real(dl), parameter :: v_0 = 25._dl ! value of constant potential
  real(dl), parameter :: m2 = 1._dl ! phi squared mass, or prefactor in power law
  integer, parameter :: inf_p = 2   ! power law of phi potential

  real(dl), parameter :: phi_p = 8._dl
  real(dl), parameter :: phi_w = 0.05_dl

  real(dl), parameter :: m2_inf = 1._dl
  real(dl), parameter :: lambda_chi = 1._dl
  real(dl), parameter :: m2_p = 1._dl
  real(dl), parameter :: g2 = 2._dl*(m2_inf-m2_p)/phi_p**2

contains

  ! Function to evaluate background potential
  ! n.b. convension on prefactor is 1/inf_p
!  elemental function bg_v(f1)
!    real(dl) :: bg_v
!    real(dl), intent(in) :: f1
!#if (POTOPT == 1)
!    bg_v = m2*f1**inf_p/dble(inf_p)
!#endif
!  end function bg_v

  function bg_v(f)
    real(dl) :: bg_v
    real(dl), intent(in) :: f(:)
#ifdef ONEFLD
    bg_v = pot_v(f(1))
#endif
#ifdef TWOFLD
    bg_v = pot_v(f(1),f(2))
#endif
  end function bg_v

  ! Function to evaluate derivative of background potential
!  elemental function bg_dv(f1)
!    real(dl) :: bg_dv
!    real(dl), intent(in) :: f1
!#if (POTOPT == 1)
!    bg_dv = m2*f1**(inf_p-1)
!#endif
!  end function bg_dv

  function bg_dv(f, ind)
    real(dl) :: bg_dv
    real(dl), intent(in) :: f(:)
    integer, intent(in) :: ind
#ifdef ONEFLD
    bg_dv = pot_dv(f(1))
#endif
#ifdef TWOFLD
    bg_dv = pot_dv(f(1),f(2),ind)
#endif
  end function bg_dv

  ! Function to evaluate 2nd derivative of background potential
!  elemental function bg_ddv(f1)
!    real(dl) :: bg_ddv
!    real(dl), intent(in) :: f1
!#if (POTOPT == 1)
!    bg_ddv = (inf_p-1)*m2*f1**(inf_p-2)
!#endif
!  end function bg_ddv

  function bg_ddv(f, ind1, ind2)
    real(dl) :: bg_ddv
    real(dl), intent(in) :: f(:)
    integer, intent(in) :: ind1, ind2
#ifdef ONEFLD
    bg_ddv = pot_ddv(f(1))
#endif
#ifdef TWOFLD
    bg_ddv = pot_ddv(f(1),f(2),ind1,ind2)
#endif
  end function bg_ddv

! One field potential functions
#ifdef ONEFLD
  elemental function pot_v(f1)
    real(dl) :: pot_v
    real(dl), intent(in) :: f1
#if (POTOPT == 0)
    pot_v = v_0
#elif (POTOPT == 1)
    pot_v = m2*f1**(inf_p)/dble(inf_p)
#endif
  end function pot_v

  elemental function pot_dv(f1)
    real(dl) :: pot_dv
    real(dl), intent(in) :: f1
#if (POTOPT == 0)
    pot_dv = 0._dl
#elif (POTOPT == 1)
    pot_dv = m2*f1**(inf_p-1)
#endif
  end function pot_dv

  elemental function pot_ddv(f1)
    real(dl) :: pot_ddv
    real(dl), intent(in) :: f1
#if (POTOPT == 0)
    pot_ddv = 0._dl
#elif (POTOPT == 1)
    pot_ddv = (inf_p-1)*m2*f1**(inf_p-2)
#endif
  end function pot_ddv
#endif

! Two field potential functions
#ifdef TWOFLD
  ! Function for potential
  elemental function pot_v(f1,f2)
    real(dl) :: pot_v
    real(dl), intent(in) :: f1,f2
#if (POTOPT == 0)
    pot_v = v_0
#elif (POTOPT == 1)
    pot_v = m2*f1**inf_p/dble(inf_p) + 0.5_dl*m2_inf*f2**2
#elif (POTOPT == 2)
    pot_v = m2*f1**inf_p/dble(inf_p) + 0.5_dl*m2_inf*f2**2 &
            - 0.5_dl*(sign(0.5_dl,-f1+phi_p-phi_w)+sign(0.5_dl,f1-phi_p-phi_w))&
            *(m2_p-m2_inf + g2*(f1-phi_p)**2 &
              -0.25_dl*g2**2/(m2_inf-m2_p)*(f1-phi_p)**4)*f2**2 &
            + 0.25_dl*lambda_chi*f2**4
#endif
  end function pot_v

  ! Function for first derivatives of the potential
  elemental function pot_dv(f1,f2,ind)
    real(dl) :: pot_dv
    real(dl), intent(in) :: f1,f2
    integer, intent(in) :: ind
#if (POTOPT == 0)
    pot_dv = 0._dl
#elif (POTOPT == 1)
    if (ind==1) then
      pot_dv = m2*f1**(inf_p-1)
    else if (ind==2) then
      pot_dv = m2_inf*f2
    endif
#elif (POTOPT == 2)
    if (ind==1) then
      pot_dv = m2*f1**(inf_p-1)&
               - (sign(0.5_dl,-f1+phi_p-phi_w)+sign(0.5_dl,f1-phi_p-phi_w))&
               *(g2*(f1-phi_p) - 0.5_dl*g2**2/(m2_inf-m2_p)*(f1-phi_p)**3)*f2**2
    else if (ind==2) then
      pot_dv = m2_inf*f2 &
               - (sign(0.5_dl,-f1+phi_p-phi_w)+sign(0.5_dl,f1-phi_p-phi_w))&
               *(m2_p-m2_inf + g2*(f1-phi_p)**2 &
               -0.25_dl*g2**2/(m2_inf-m2_p)*(f1-phi_p)**4)*f2 &
               + lambda_chi*f2**3
    endif
#endif
  end function pot_dv

  ! Function for second derivatinves of the potential
  elemental function pot_ddv(f1,f2,ind1,ind2)
    real(dl) :: pot_ddv
    real(dl), intent(in) :: f1,f2
    integer, intent(in) :: ind1,ind2
#if (POTOPT == 0)
    pot_ddv = 0._dl
#elif (POTOPT == 1)
    if (ind1==1 .and. ind2==1) then
      pot_ddv = (inf_p-1)*m2*f1**(inf_p-2)
    else if (ind1==2 .and. ind2==2) then
      pot_ddv = m2_inf
    else
      pot_ddv = 0._dl
    endif
#elif (POTOPT == 2)
    if (ind1==1 .and. ind2==1) then
      pot_ddv = (inf_p-1)*m2*f1**(inf_p-2) &
                - (sign(0.5_dl,-f1+phi_p-phi_w)+sign(0.5_dl,f1-phi_p-phi_w))&
                *(g2 - 1.5_dl*g2**2/(m2_inf-m2_p)*(f1-phi_p)**2)*f2**2
    else if (ind1==2 .and. ind2==2) then
      pot_ddv = m2_inf&
                - (sign(0.5_dl,-f1+phi_p-phi_w)+sign(0.5_dl,f1-phi_p-phi_w))&
                *(m2_p-m2_inf + g2*(f1-phi_p)**2 &
                -0.25_dl*g2**2/(m2_inf-m2_p)*(f1-phi_p)**4) &
                + 3*lambda_chi*f2**2
    else
      pot_ddv = - (sign(1._dl,-f1+phi_p-phi_w)+sign(1._dl,f1-phi_p-phi_w))&
                *(g2*(f1-phi_p) - 0.5_dl*g2**2/(m2_inf-m2_p)*(f1-phi_p)**3)*f2
    endif
#endif
  end function pot_ddv
#endif

end module potential_mod

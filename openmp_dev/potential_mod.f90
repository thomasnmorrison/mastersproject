!!!!! potential_mod.f90 !!!!!
!!!!! Module to include the definiton of the potential and it's derivatives.
! The intention is to put functions asociated with the potential into this
! module, then include this module in the hamiltonian module.
! This may require moving some of the parameter definitions from the hamiltonian
! module over to here.

! to do: remove definition of potential and modeldv from hamiltonian
! to do: add the ability to define multiple phi_p for multiple interactions
! to do: longitudinal instability
! to do: the optional arguments for nfld=1 vs nfld=2 is too cumbersum, replace with _1fld, _2fld versions and use preprocessor statement 
! to do: write a subroutine to calculate the trace of the mass matrix

module potential_mod
#include "macros.h"
  use fftw3
  use params

  implicit none

#ifdef ONEFLD
	integer, parameter :: nfld=1
#endif
#ifdef TWOFLD2
	integer, parameter :: nfld=2
#endif
	integer, parameter :: potential_option = 6 ! parameter to choose the form of Delta_V
	! no phi-chi interaction, chi massless: potential_option = 0
	! trapped plus transverse instability: potential_option = 1
	! transverse instability blip: potential_option = 2
	! absolute value transverse blip: potential_option = 3
	! trapped-asymptotically constant: potential_option = 4
	! constant chi mass: potential_option = 5
	! trapped-asymptotically constant (polynomial): potential_option = 6
	integer, parameter :: infl_option = 2 ! parameter to choose form of inflationary potential
	! phi^4: infl_option = 1
	! m^2 phi^2: infl_option = 2
	! constant potential: infl_option = 3

	! general parameters
	real(dl), parameter :: mpl=1.e5												! machine units conversion factor, moved from hampiltonian_conformal.f90
	real(dl), parameter :: m2 = 1.0_dl										! inflaton mass
	real(dl), parameter :: lambda_chi = 10.0_dl!10.0_dl						! chi self interaction
	real(dl), parameter :: phi_p = 8.75_dl!*dsqrt(4.0*twopi)	! interaction potential characteristic phi value

	! potential_option = 1 parameters
	!real(dl), parameter :: g2 = 1.e5!(1.25e4)!/(4.*twopi)!											! phi-chi coupling strength
	real(dl), parameter :: beta2 = 0.!1.e-3*g2							! instability amplitude

	! potential_option = 2 parameters
	real(dl), parameter :: a2 =	1.2e4											! max/min of deltaV
	real(dl), parameter :: b_1 = 0.3_dl										!	width of deltaV (std of Gaussian envelope)
	real(dl), parameter :: c_1 = a2/b_1*exp(0.5_dl)				!	normalization factor
	real(dl), parameter :: c_2 = 1.0_dl/(2.0_dl*b_1**2)		! derived parameter for deltaV

	! potential_option = 4 parameters
	real(dl), parameter :: m2_inf = 1._dl*m2										! asymptotic squared mass of chi
	real(dl), parameter :: m2_p =	-100._dl*m2										! minimum squared mass of chi (at phi=phi_p) 
	! Derived parameters
	real(dl), parameter :: c_3 = 0._dl!2.*sqrt(3.*g2/(m2_inf-m2_p))
	! Also set g2 under "potential_option = 1 parameters"

	! potential_option = 5 parameters
	! set m2_inf under "potential_option = 4 parameters"

	! potential_option = 6 parameters
	real(dl), parameter :: c_4 = 0.05_dl															! width of delta V, with g2 as a derived parameter
	!real(dl), parameter :: c_4 = sqrt(2._dl*(m2_inf-m2_p)/g2)		! width of delta V as a derived parameter
	real(dl), parameter :: g2 = 2._dl*(m2_inf-m2_p)/c_4**2				! coupling strength as a derived parameter
	! Also set g2, m2_inf, m2_p, and phi_p

contains

	! Function to calculate the background potential (no phi chi interaction term)
	elemental function background_V(f1,f2)
    real(dl) :: background_V
    real(dl), intent(in) :: f1,f2
		
		if (infl_option==1) then
			background_V = 0.25*f1**4 + 0.25_dl*lambda_chi*f2**4
		elseif (infl_option==2) then
			background_V = 0.5_dl*m2*f1**2 + 0.25_dl*lambda_chi*f2**4
		endif
  end function background_V

	! function to calculate potential of the ghost fields
	elemental function ghost_V(f1,f2)
		real(dl) :: ghost_V
		real(dl), intent(in) :: f1,f2

		ghost_V = 0.5_dl*m2*f1**2 + 0.5_dl*m2_inf*f2**2
	end function ghost_V

	! function to calculate derivatives of ghost fields potential
	elemental function ghost_dV(f1,f2,ind)
		real(dl) :: ghost_dV
		real(dl), intent(in) :: f1,f2
		integer, intent(in) :: ind

		if (ind==1) then
			ghost_dV = m2*f1
		elseif (ind==2) then
			ghost_dV = m2_inf*f2
		endif
	end function ghost_dV

	! Function to calculate the background potential for 1 field models
	elemental function background_V_1fld(f1)
		real(dl) :: background_V_1fld
		real(dl), intent(in) :: f1

		if (infl_option==1) then
			background_V_1fld = 0.25_dl*f1**4
		elseif (infl_option==2) then
			background_V_1fld = 0.5_dl*m2*f1**2
		elseif (infl_option==3) then
			background_V_1fld = 0.5_dl*14.2_dl**2
		endif
	end function background_V_1fld

	! Function to calculate derivatives of the background potential
	elemental function background_dV(f1,f2,ind)
    real(dl) :: background_dV
    real(dl), intent(in) :: f1,f2
		integer, intent(in) :: ind
		
		if (infl_option==1) then
			if (ind==1) then
				background_dV = f1**3
			elseif (ind==2) then
				background_dV = lambda_chi*f2**3
			endif
		elseif (infl_option==2) then
			if (ind==1) then
				background_dV = m2*f1
			elseif (ind==2) then
				background_dV = lambda_chi*f2**3
			endif
		endif
  end function background_dV

	! Function to calculate derivative of background potential for 1 field models
	elemental function background_dV_1fld(f1)
		real(dl) :: background_dV_1fld
		real(dl), intent(in) :: f1

		if (infl_option==1) then
			background_dV_1fld = f1**3
		elseif (infl_option==2) then
			background_dV_1fld = m2*f1
		elseif (infl_option==3) then
			background_dV_1fld = 0._dl
		endif
	end function background_dV_1fld

	! Function to calculate the interaction potential
	elemental function Delta_V(f1,f2)
    real(dl) :: Delta_V
    real(dl), intent(in) :: f1,f2
		
		if (potential_option==0) then
			Delta_V = 0.0_dl
		elseif (potential_option==1) then
			Delta_V = (0.5_dl*g2*(f1-phi_p)**2 - 0.5_dl*beta2)*f2**2
		elseif (potential_option==2) then
			Delta_V = - c_1*(f1-phi_p)*exp(-c_2*(f1-phi_p)**2)*f2**2
		elseif (potential_option==3) then
			Delta_V = abs(c_1*(f1-phi_p))*exp(-c_2*(f1-phi_p)**2)*f2**2
		elseif (potential_option==4) then
			!Delta_V = 0.5_dl * (m2_inf - 6.*(m2_inf-m2_p) / (5.+dcosh(2.*dsqrt(3.*g2)*(f1-phi_p)/dsqrt(m2_inf-m2_p))) )*f2**2
			Delta_V = 0.5_dl * (m2_inf - 6.*(m2_inf-m2_p) / (5.+dcosh(c_3*(f1-phi_p))) )*f2**2
		elseif (potential_option==5) then
			Delta_V = 0.5 * m2_inf * f2**2
		elseif (potential_option==6) then
			Delta_V = 0.5_dl * (m2_inf + (-sign(0.5_dl,-f1+phi_p-c_4)-sign(0.5_dl,f1-phi_p-c_4))*( m2_p-m2_inf + g2*(f1-phi_p)**2 - g2**2/(4._dl*(m2_inf-m2_p))*(f1-phi_p)**4)) * f2**2
		endif
  end function Delta_V

	! Function to calculate the derivatives of the interaction potential
	elemental function Delta_dV(f1,f2,ind)
    real(dl) :: Delta_dV
    real(dl), intent(in) :: f1,f2
		integer, intent(in) :: ind

		if (potential_option==0) then
			Delta_dV = 0.0_dl
		elseif (potential_option==1) then
			if (ind==1) then
				Delta_dV = g2*(f1-phi_p)*f2**2
			endif
			if (ind==2) then
				Delta_dV = (g2*(f1-phi_p)**2 - beta2)*f2
			endif
		! Double check these derivatives
		elseif (potential_option==2) then
			if (ind==1) then
				Delta_dV = - c_1*(1.0_dl - 2.0_dl*c_2*(f1-phi_p)**2)*exp(-c_2*(f1-phi_p)**2)*f2**2
			endif
			if (ind==2) then
				Delta_dV =  - 2.0_dl*c_1*(f1-phi_p)*exp(-c_2*(f1-phi_p)**2)*f2
			endif
		elseif (potential_option==3) then
			if (ind==1) then
				if (f1>=phi_p) then
					Delta_dV = c_1*(1.0_dl - 2.0_dl*c_2*(f1-phi_p)**2)*exp(-c_2*(f1-phi_p)**2)*f2**2
			 	elseif (f1<phi_p) then
					Delta_dV = - c_1*(1.0_dl - 2.0_dl*c_2*(f1-phi_p)**2)*exp(-c_2*(f1-phi_p)**2)*f2**2
			 	endif
			endif
			if (ind==2) then
				Delta_dV = 2.0_dl*abs(c_1*(f1-phi_p)*exp(-c_2*(f1-phi_p)**2))*f2
			endif
		elseif (potential_option==4) then
			if (ind==1) then
				!Delta_dV = (6.*dsqrt(3.*g2*(m2_inf-m2_p)) * dsinh(2.*dsqrt(3.*g2)*(f1-phi_p)/dsqrt(m2_inf-m2_p)) / (5.+dcosh(2.*dsqrt(3.*g2)*(f1-phi_p)/dsqrt(m2_inf-m2_p)))**2 ) * f2**2
				Delta_dV = (6.*dsqrt(3.*g2*(m2_inf-m2_p)) * dsinh(c_3*(f1-phi_p)) / (5.+dcosh(c_3*(f1-phi_p)))**2 ) * f2**2
			endif
			if (ind==2) then
				Delta_dV = (m2_inf - 6.*(m2_inf-m2_p)/(5.+dcosh(c_3*(f1-phi_p))))*f2
			endif
		elseif (potential_option==5) then
			if (ind==1) then
				Delta_dV = 0.0
			endif
			if (ind==2) then
				Delta_dV = m2_inf*f2
			endif
		elseif (potential_option==6) then
			if (ind==1) then
				Delta_dV = 0.5_dl * ((-sign(0.5_dl,-f1+phi_p-c_4)-sign(0.5_dl,f1-phi_p-c_4))*( 2._dl*g2*(f1-phi_p) - g2**2/(m2_inf-m2_p)*(f1-phi_p)**3)) * f2**2
			elseif (ind==2) then
				Delta_dV = (m2_inf + (-sign(0.5_dl,-f1+phi_p-c_4)-sign(0.5_dl,f1-phi_p-c_4))*( m2_p-m2_inf + g2*(f1-phi_p)**2 - g2**2/(4._dl*(m2_inf-m2_p))*(f1-phi_p)**4)) * f2
			endif
		endif
  end function Delta_dV

	! Function to sum the background and interaction potentials and be called from outside this module
	elemental function potential(f1,f2)
		real(dl) :: potential
    real(dl), intent(in) :: f1,f2
	
		potential = background_V(f1,f2)+Delta_V(f1,f2)
	end function potential

	! Dev function to wrap both 1 field and 2 field potential functions
	elemental function potential_test(f1,f2)
		real(dl) :: potential_test
		real(dl), intent(in) :: f1
		real(dl), optional, intent(in) :: f2

		if (present(f2)) then
			potential_test = background_V(f1,f2)+Delta_V(f1,f2) ! 2 field case
		else
			potential_test = background_V_1fld(f1)	! 1 field case
		endif
	end function potential_test

	! Function to sum derivatived of the background and interaction potentials and be called from outside this module
	elemental function modeldv(f1,f2,ind)
    real(dl) :: modeldv
    real(dl), intent(in) :: f1,f2
		integer, intent(in) :: ind

		modeldv = background_dV(f1,f2,ind) + Delta_dV(f1,f2,ind)
	end function modeldv

	! Dev function to wrap 1 field and 2 field potential derivatives
	! n.b. rearranged order of arguments
	elemental function modeldv_test(ind,f1,f2)
		real(dl) :: modeldv_test
		real(dl), intent(in) :: f1
		integer, intent(in) :: ind
		real(dl), optional, intent(in) :: f2
		
		if (present(f2)) then
			modeldv_test = background_dV(f1,f2,ind) + Delta_dV(f1,f2,ind)
		else
			modeldv_test = background_dV_1fld(f1)
		endif
	end function modeldv_test

	! Function to calculate the trace of the mass matrix
	function trace_m2(f1,f2)
		real(dl) :: trace_m2
		real(dl), intent(in) :: f1(:,:,:),f2(:,:,:)

		real(dl), dimension(nfld) :: m2_infl
		real(dl), dimension(nfld) :: m2_delta
		real(dl) :: f1_hom, f2_hom
		integer :: i

		f1_hom = sum(f1(IRANGE)); f2_hom = sum(f2(IRANGE))

		if (infl_option==1) then
		elseif (infl_option==2) then
			m2_infl(1) = m2
			m2_infl(2) = 0.0
		if (potential_option==1) then
		elseif (potential_option==2) then
		elseif (potential_option==3) then
		elseif (potential_option==4) then
			m2_delta(1) = (36.*g2)*f2_hom**2*((-2.*sinh(c_3*(f1_hom-phi_p))**2)/(5 + cosh(c_3*(f1_hom-phi_p))) + cosh(c_3*(f1_hom-phi_p)))/(5 + cosh(c_3*(f1_hom-phi_p)))**2
			m2_delta(2) = (m2_inf - 6.*(m2_inf-m2_p) / (5.+dcosh(c_3*(f1_hom-phi_p))) )
		elseif (potential_option==5) then
			m2_delta(1) = 0.0
			m2_delta(2) = m2_inf
		endif; endif
		
		trace_m2 = sum(m2_infl) + sum(m2_delta) 
	end function trace_m2

end module potential_mod









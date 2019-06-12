! correlator_int_mod.f90

! Moduled to generate initial power spectra by integration of mode function from Minkowski initial conditions

! to do: write a subroutine to integrate mode functions in conformal time
! to do: write a subroutine to get a_scl from a background solution
! to do: fix nos so that it covers the corners of the lattice - done, but covers modes not on lattice
! to do: double check that set_cor_cosmic(j) is correct

module correlator_int_mod

#include "macros.h"
	use params
	use analysis
	use potential_mod

integer, parameter :: os = 16 									! factor by which k is oversampled relative to lattice
integer, parameter :: nos = max(nx,ny,nz)*os		! number times dkos to reach max k on lattice, as is this has an implicit top hat filter at he Nyquist frequency
real(dl), parameter :: dkos = dk/dble(os)				! oversampled mode spacing

integer, parameter :: nstep_mode = 2**8							! number of intgration steps

real(dl), dimension(2*nfld,2*nfld,nos) :: cor_rad						! power on some oversampled radial profile
real(dl), dimension(nfld) :: m2_diag															! masses of fields
real(dl) :: t_bg																						! time variable
real(dl) :: dt_bg																						! time step size
real(dl), dimension(nfld) :: mode_amp												! mode function amplitude
real(dl), dimension(nfld) :: dmode_amp											! time derivative of mode function amplitude
real(dl) :: a_bg																						! scale factor of background
real(dl) :: hub_bg																					! Hubble parameter of background

contains

	! subroutine to initialize the m2_diag vector
	! n.b. assume constant masses, assume fields are eigen vectors of m2 matrix
	subroutine init_m2_diag()
		integer :: n			! field index
		do n=1, nfld
			if (n==1) then
				m2_diag(n) = m2
			elseif (n==2) then
				m2_diag(n) = m2_inf
			else
				print*, "init_m2_diag not build for nfld>2"
			endif
		enddo
	end subroutine init_m2_diag

	! subroutine to update background scale factor given the clock variable
!	subroutine get_a_bg(clock, clock_opt)
!		real(dl) :: clock			! clock varialble, ie. cosmic time, conformal time, inflaton
!		integer :: clock_opt	! option to choose between type of clock (maybe make as a preprosseor)
!	end subroutine get_a_bg

	! subroutine to step background, step cosmic time, step a_bg, step hub_bg
	! n.b. written for constant hubble
	subroutine background_evolve_cosmic()
		t_bg = t_bg + dt_bg
		a_bg = exp(hub_bg*t_bg)		! can include number of e-folds in here or not
		hub_bg = hub_bg						! assuming constant hubble background
	end subroutine background_evolve_cosmic

	! subroutine to initialize cor_rad by integrating from Minkowski solution in cosmic time
	subroutine init_cor_rad_cosmic(n_ef, hub_in)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl), intent(in) :: hub_in		! input initial hubble parameter
		real(dl) :: k2_mode								! comoving wavenumber at start of lattice sim
		integer :: i											! track integration step
		integer :: j											! index k_mode

		real(dl), parameter :: norm = nvol/(mpl**2*dx**3)		! normalization constant
		
		cor_rad(:,:,:) = 0._dl
		call init_m2_diag()
		! loop this over k_mode
		do j=1, nos
			k2_mode = (j*dkos)**2
			call mode_init_cosmic(n_ef, hub_in, k2_mode)
			! loop this over integration steps
			do i=1,nstep_mode
				call mode_evolve_cosmic()
			enddo
			call set_cor_cosmic(n_ef,j)
		enddo
		
		cor_rad(:,:,:) = norm*cor_rad(:,:,:)
	end subroutine init_cor_rad_cosmic

	! subroutine to initialize field variables for cosmic time integration
	subroutine mode_init_cosmic(n_ef, hub_in, k2)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl), intent(in) :: hub_in		! input initial hubble parameter
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)
		
		t_bg = 0._dl
		dt_bg = n_ef/hub_bg/dble(nstep_mode)
		a_bg = 1.0_dl
		hub_bg = hub_in
		mode_amp(:) = 1._dl/sqrt(2._dl*sqrt(m2_diag(:)+k2*exp(2._dl*n_ef)))
		dmode_amp(:) = 0._dl
	end subroutine mode_init_cosmic

	! subroutine to perform one integration step on the mode function
	subroutine mode_evolve_cosmic()
		mode_amp(:) = mode_amp(:) + 0.5_dl*dmode_amp(:)*dt_bg
		dmode_amp(:) = dmode_amp(:) + 0.5_dl*(0.25_dl/mode_amp(:)**3 - (k2*exp(2._dl*n_ef)/a_bg**2 + m2_diag(:) -2.25_dl*hub_bg**2)*mode_amp(:))*dt_bg
		call background_evolve_cosmic()
		dmode_amp(:) = dmode_amp(:) + 0.5_dl*(0.25_dl/mode_amp(:)**3 - (k2*exp(2._dl*n_ef)/a_bg**2 + m2_diag(:) -2.25_dl*hub_bg**2)*mode_amp(:))*dt_bg
		mode_amp(:) = mode_amp(:) + 0.5_dl*dmode_amp(:)*dt_bg
	end subroutine mode_evolve_cosmic

	! subroutine to set power from mode functions in cosmic time
	subroutine set_cor_cosmic(n_ef, l)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		integer :: l			! mode indexing
		integer :: m, n		! field index

		do n=1,nfld
			cor_rad(n, n, l) = mode_amp(n)**2/exp(-n_ef)**3/a_bg**3
			cor_rad(n+1, n+1, l) = exp(-n_ef)**3*a_bg**3*(2.25_dl*hub_bg**2*mode_amp(n)**2 + dmode_amp(n)**2 - 3._dl*hub_bg*mode_amp(n)*dmode_amp(n) + 0.25_dl/mode_amp(n)**2)
			cor_rad(n+1, n, l) = -1.5_dl*hub*mode_amp(n)**2 + dmode_amp(n)*mode_amp(n)
			cor_rad(n, n+1, l) = cor_rad(n+1, n, l)
		enddo
	end subroutine set_cor_cosmic

end module correlator_int_mod

! correlator_int_mod.f90

! Moduled to generate initial power spectra by integration of mode function from Minkowski initial conditions
! n.b. The calculation of zeta modes assumes contribution is purely from the phi field. Should generalize to the multifield case.

! to do: write a subroutine to integrate mode functions in conformal time
! to do: write a subroutine to get a_scl from a background solution
! to do: fix nos so that it covers the corners of the lattice - done, but covers modes not on lattice
! to do: 
! to do: check the normalization on the zeta mode
! to do:
! to do: 
! to do: define background dphi, rho, and P

module correlator_int_mod

#include "macros.h"
	use params
	use analysis
	use potential_mod

integer, parameter :: os_mode = 2**4 									! factor by which k is oversampled relative to lattice
integer, parameter :: nos_mode = max(nx,ny,nz)*os_mode		! number times dkos to reach max k on lattice, as is this has an implicit top hat filter at he Nyquist frequency
real(dl), parameter :: dkos_mode = dk/dble(os_mode)					! oversampled mode spacing

integer, parameter :: nstep_mode = 2**14									! number of intgration steps

real(dl), dimension(2*nfld,2*nfld,nos_mode) :: cor_rad			! power on some oversampled radial profile
real(dl), dimension(nos_mode)	:: z_amp_rad									! zeta mode amplitudes on oversampled radial profile
real(dl), dimension(nos_mode) :: z_shift_rad								! zeta mode phase shifts on oversampled radial profile

real(dl), dimension(nfld) :: m2_diag												! masses of fields
real(dl) :: t_bg																						! time variable
real(dl) :: dt_bg 																					! time step size
real(dl) :: a_bg																						! scale factor of background
real(dl) :: hub_bg																					! Hubble parameter of background
real(dl) :: dphi_bg																					! time derivative of background inflaton field

real(dl), dimension(nfld) :: mode_amp												! mode function amplitude for fields
real(dl), dimension(nfld) :: dmode_amp											! time derivative of mode function amplitude for fields
real(dl) :: zmode_amp																				! mode function amplitude for zeta
real(dl) :: dzmode_amp																			! time derivative of mode function amplitude for zeta
real(dl) :: zphase_shift																		! phase difference: zeta phase - inflaton phase
real(dl) :: dzphase_shift																		! time derivative of phase difference


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
		!real(dl), parameter :: norm = 1._dl/(mpl**2*dx**3*nvol)		! normalization constant
		
		cor_rad(:,:,:) = 0._dl
		call init_m2_diag()
		! loop this over k_mode
		do j=1, nos_mode
			k2_mode = (j*dkos_mode)**2
			call mode_init_cosmic(n_ef, hub_in, k2_mode)
			!print*, 'j, k2, Det cor = ', j, k2_mode, mode_amp(2 )
			! loop this over integration steps
			!print*, 'A, dA = ', mode_amp(2), dmode_amp(2)
			!do i=1,nstep_mode
				call mode_evolve_cosmic(n_ef,k2_mode)
				!if (modulo(i,2**6)==0) then
				!	print*, 't, A, dA = ', t_bg, mode_amp(2), dmode_amp(2)
				!endif
			!enddo
			!print*, 'A, dA = ', mode_amp(2), dmode_amp(2)
			call set_cor_cosmic(n_ef,j)
			!print*, 'a = ', exp(-n_ef)*a_bg
			!print*, 'k = ', sqrt(k2_mode)
			!print*, 'cor = ', cor_rad(3,3,j), cor_rad(4,4,j)
			!print*, (exp(-n_ef)**3*a_bg**3)
			!print*, 'Det cor = ', cor_rad(3,3,j),cor_rad(4,4,j),cor_rad(3,4,j),cor_rad(4,3,j), cor_rad(3,3,j)*cor_rad(4,4,j) - cor_rad(3,4,j)*cor_rad(4,3,j)
		enddo
		cor_rad(:,:,:) = norm*cor_rad(:,:,:)	
	end subroutine init_cor_rad_cosmic

	! subroutine to initialize cor_rad by integrating from Minkowski solution in cosmic time. Also sets initial zeta
	subroutine init_cor_rad_cosmic_z(n_ef, dphi_in, hub_in)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl), intent(in) :: dphi_in		! input initial time derivative of inflaton
		real(dl), intent(in) :: hub_in		! input initial hubble parameter
		real(dl) :: k2_mode								! comoving wavenumber at start of lattice sim
		integer :: i											! track integration step
		integer :: j											! index k_mode

		real(dl), parameter :: norm = nvol/(mpl**2*dx**3)		! normalization constant for the fields
		real(dl), parameter :: norm_z = 1._dl								! normalization constant for zeta
		
		cor_rad(:,:,:) = 0._dl
		call init_m2_diag()
		! loop this over k_mode
		do j=1, nos_mode
			k2_mode = (j*dkos_mode)**2
			call mode_init_cosmic(n_ef, hub_in, k2_mode)
			call mode_init_cosmic_z(n_ef, dphi_in, k2_mode)
			call mode_evolve_cosmic_z(n_ef,k2_mode)
			call set_cor_cosmic(n_ef,j)
			call set_zeta_cosmic(j)
		enddo
		cor_rad(:,:,:) = norm*cor_rad(:,:,:)
		z_amp_rad(:) = norm_z*z_amp_rad(:)
		z_shift_rad(:) = norm_z*z_shift_rad(:)
	end subroutine init_cor_rad_cosmic_z

	! subroutine to initialize field variables for cosmic time integration
	! n.b. mode_amp = a**(3./2.)*phi_k, not (a/a_0)**(3./2.)*phi_k like worked out pencil and paper
	subroutine mode_init_cosmic(n_ef, hub_in, k2)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl), intent(in) :: hub_in		! input initial hubble parameter
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)
		
		t_bg = 0._dl
		dt_bg = n_ef/(hub_in*dble(nstep_mode))
		a_bg = 1.0_dl
		hub_bg = hub_in
		mode_amp(:) = 1._dl/sqrt(2._dl*sqrt(m2_diag(:)+k2*exp(2._dl*n_ef)))
		dmode_amp(:) = 0._dl
	end subroutine mode_init_cosmic

	! subroutine to initialize the zeta variables for cosmic time integration
	! n.b. must be called after mode_init_cosmic
	! to do: set zmode_amp
	! to do: after testing merge this subroutine with mode_init_cosmic
	subroutine mode_init_cosmic_z(n_ef, dphi_in, k2)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl), intent(in) :: dphi_in		! input initial time derivative of inflaton
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)

		dphi_bg = dphi_in
		zmode_amp = -2._dl*k2*exp(2._dl*n_ef)*mode_amp(1)**3/(3._dl*dphi_bg)!-2._dl*exp(2._dl*n_ef)*k2*mode_amp(1)**3/(3._dl*dphi_bg)
		dzmode_amp = 0._dl
		zphase_shift = 0.25_dl*twopi!-0.25_dl*twopi
		dzphase_shift = 0._dl

		!print*, "k*a_0 = ", sqrt(k2)*exp(n_ef)
		!print*, "mode_amp(1) initial = ", mode_amp(1)
		!print*, "zmode_amp initial = ", zmode_amp
		!print*, "zphase_shift initial = ", zphase_shift

	end subroutine mode_init_cosmic_z

	! subroutine to perform one integration step on the mode function
	subroutine mode_evolve_cosmic(n_ef, k2)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)
		integer :: i

		do i=1,nstep_mode
			mode_amp(:) = mode_amp(:) + 0.5_dl*dmode_amp(:)*dt_bg
			dmode_amp(:) = dmode_amp(:) + 0.5_dl*(0.25_dl/mode_amp(:)**3 - (k2*exp(2._dl*n_ef)*exp(-2._dl*hub_bg*t_bg) + m2_diag(:) - 2.25_dl*hub_bg**2)*mode_amp(:))*dt_bg
			!if (modulo(i,2**4)==1) then
			!	print*, 't, A, dA = ', t_bg, mode_amp(2), dmode_amp(2)
			!endif
			call background_evolve_cosmic()
			dmode_amp(:) = dmode_amp(:) + 0.5_dl*(0.25_dl/mode_amp(:)**3 - (k2*exp(2._dl*n_ef)*exp(-2._dl*hub_bg*t_bg) + m2_diag(:) - 2.25_dl*hub_bg**2)*mode_amp(:))*dt_bg
			mode_amp(:) = mode_amp(:) + 0.5_dl*dmode_amp(:)*dt_bg
		enddo
	end subroutine mode_evolve_cosmic

	! subroutine to perform one integration step of the mode function and zeta mode
	! n.b. the intention is that this subroutine will supercede mode_evolve_cosmic once tested
	! to do:
	subroutine mode_evolve_cosmic_z(n_ef, k2)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)
		integer :: i

		do i=1,nstep_mode
			!print*, "zmode_amp = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			mode_amp(:) = mode_amp(:) + 0.5_dl*dmode_amp(:)*dt_bg
			!print*, "zmode_amp1 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			!print*, dt_bg
			dmode_amp(:) = dmode_amp(:) + 0.5_dl*(0.25_dl/mode_amp(:)**3 - (k2*exp(2._dl*n_ef)*exp(-2._dl*hub_bg*t_bg) + m2_diag(:) - 2.25_dl*hub_bg**2)*mode_amp(:))*dt_bg
			!print*, "zmode_amp2 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			zmode_amp = zmode_amp + 0.5_dl*dt_bg*dzmode_amp
			!print*, dt_bg
			!print*, "zmode_amp3 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			call dzmode_amp_step_cosmic(0.5_dl*dt_bg, k2)
			!print*, "zmode_amp4 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			zphase_shift = zphase_shift + 0.5_dl*dt_bg*dzphase_shift
			!print*, "zmode_amp5 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			call dzphase_shift_step_cosmic(0.5_dl*dt_bg, n_ef, k2)
			!print*, "zmode_amp6 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			call background_evolve_cosmic()
			!print*, "zmode_amp7 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			call dzphase_shift_step_cosmic(0.5_dl*dt_bg, n_ef, k2)
			!print*, "zmode_amp8 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			zphase_shift = zphase_shift + 0.5_dl*dt_bg*dzphase_shift
			!print*, "zmode_amp9 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			call dzmode_amp_step_cosmic(0.5_dl*dt_bg, k2)
			!print*, "zmode_amp10 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			zmode_amp = zmode_amp + 0.5_dl*dt_bg*dzmode_amp
			!print*, "zmode_amp11 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			dmode_amp(:) = dmode_amp(:) + 0.5_dl*(0.25_dl/mode_amp(:)**3 - (k2*exp(2._dl*n_ef)*exp(-2._dl*hub_bg*t_bg) + m2_diag(:) - 2.25_dl*hub_bg**2)*mode_amp(:))*dt_bg
			!print*, "zmode_amp12 = ", zmode_amp
			!print*, "dzmode_amp = ", dzmode_amp
			mode_amp(:) = mode_amp(:) + 0.5_dl*dmode_amp(:)*dt_bg
		enddo
	end subroutine mode_evolve_cosmic_z

	! subroutine to perform an integration step on dzphase_shift
	! n.b. Assumes rho + P = dphi**2
	subroutine dzphase_shift_step_cosmic(dt_in, n_ef, k2)
		real(dl) :: dt_in
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)
		dzphase_shift = dzphase_shift + dt_in*(-0.5_dl*exp(-3._dl*n_ef)/mode_amp(1)**2 + sin(zphase_shift)*k2*dphi_bg*mode_amp(1)/(3._dl*dphi_bg**2*a_bg**3.5*zmode_amp))
	end subroutine dzphase_shift_step_cosmic

	! subroutine to perform an integration step on dzmode_amp
	! n.b. Assumes rho + P = dphi**2
	subroutine dzmode_amp_step_cosmic(dt_in, k2)
		real(dl) :: dt_in
		real(dl) :: k2										! comoving wavenumber squared (with a=1 at start of sim, end of mode integration)
		dzmode_amp = dzmode_amp - dt_in*(cos(zphase_shift) * k2 * dphi_bg*mode_amp(1)) / (3._dl*dphi_bg**2) / (a_bg**3.5)
	end subroutine dzmode_amp_step_cosmic

	! subroutine to set power from mode functions in cosmic time
	subroutine set_cor_cosmic(n_ef, l)
		real(dl), intent(in) :: n_ef			! number of e-folds before simulation start when modes are matched onto Minkowski space solution
		integer :: l			! mode indexing
		integer :: m, n		! field index

		do n=1,nfld
			cor_rad(2*n-1, 2*n-1, l) = mode_amp(n)**2/(exp(-n_ef)**3*a_bg**3)
			cor_rad(2*n, 2*n, l) = exp(-n_ef)**3*a_bg**3*(2.25_dl*hub_bg**2*mode_amp(n)**2 + dmode_amp(n)**2 - 3._dl*hub_bg*mode_amp(n)*dmode_amp(n) + 0.25_dl/mode_amp(n)**2)
			cor_rad(2*n, 2*n-1, l) = -1.5_dl*hub_bg*mode_amp(n)**2 + dmode_amp(n)*mode_amp(n)
			cor_rad(2*n-1, 2*n, l) = cor_rad(2*n, 2*n-1, l)
		enddo
	end subroutine set_cor_cosmic

	! subroutine to store computed zeta modes, should be called along with set_cor_cosmic
	! to do:
	subroutine set_zeta_cosmic(l)
		integer :: l											! mode indexing
		
		z_amp_rad(l) = zmode_amp!/mode_amp(1)		! the division is so zeta can be set linear in phi_k instead of the grv directly
		z_shift_rad(l) = zphase_shift
		!print*, "l = ", l
		!print*, "mode_amp(1) = ", mode_amp(1)
		!print*, "z_amp_rad = ", z_amp_rad(l)
		!print*, "z_shift_rad = ", z_shift_rad(l)
	end subroutine set_zeta_cosmic

end module correlator_int_mod

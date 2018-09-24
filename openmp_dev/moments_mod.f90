!!!!! moments_mod.f90 !!!!!
!!!!! module that allows for the calulation of statistical moments

! to do: include required modules
! to do: update Makefile

module moments_mod

! Inculde other modules
#include "macros.h"
	use params

contains

	! function to calculate the mean of f1 over the lattice
	function get_mean_3d(f1, nsize)
		real(dl) :: get_mean_3d
		real(dl), intent(in) :: f1(:,:,:)
		integer, dimension(1:3) :: nsize

		integer :: n1, n2, n3
		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		get_mean_3d = sum(f1)/dble(n1)/dble(n2)/dble(n3)

	end function get_mean_3d

	! function to compute the n_mom moment of f1 around cm 
	function get_moment_3d(f1, cm, n_moment, nsize)
		real(dl) :: get_moment_3d
		real(dl), intent(in) :: f1(:,:,:)
		real(dl), intent(in) :: cm
		integer, intent(in) :: n_moment
		integer, dimension(1:3) :: nsize

		integer :: n1, n2, n3
		integer :: i, j, k
		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
	
		get_moment_3d = sum((f1(:,:,:)-cm)**n_moment)/dble(n1)/dble(n2)/dble(n3)

	end function get_moment_3d



end module moments_mod

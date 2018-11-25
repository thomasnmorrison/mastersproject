!!!!! sample_sites_mod.f90 !!!!!
!!!!! functionality for sampling at different lattice sites !!!!!

! to do: select sample sites without replacement

module sample_sites_mod
#include "macros.h"
	use hamiltonian

	implicit none

	integer, parameter :: n_sample = 32!512! !1!
	! Hard coded for three dimensions
	integer, dimension(n_sample,3) :: sample_site
	real(dl), dimension(n_sample) :: sample_zeta
	real(dl), dimension(nfld,n_sample) :: sample_zeta_f !evolve zeta for each field
	real(dl), dimension(2,nfld,n_sample) :: sample_dzeta_f
	real(dl), dimension(2,n_sample) :: sample_dzeta
	real(dl), dimension(nfld, n_sample) :: sample_field
	real(dl), dimension(nfld, n_sample) :: sample_fieldp

contains

! to do: make alternate sample site read in subroutine
	subroutine init_sites()
		integer :: i
		integer :: j
		real(dl) :: site
		integer :: k
		if (run_hom==.FALSE.) then
			k = n_sample
		elseif (run_hom==.TRUE.) then
			k = 1
		endif
		open(unit=93,file="sample_sites.out")
		!do i=1,n_sample
		do i=1,k
			do j=1,3
				call random_number(site)
				site = nx*site
				sample_site(i,j) = int(site) + 1
			enddo
			write(93,'(30(I5,2x))') sample_site(i,:)
		enddo
	end subroutine init_sites

	subroutine read_sites()
		integer :: i
		integer :: k
		if (run_hom==.FALSE.) then
			k = n_sample
		elseif (run_hom==.TRUE.) then
			k = 1
		endif
		open(unit=93,file="sample_sites.out")
		!open(unit=92,file="sample_sites_test.out")
		!do i=1,n_sample
		do i=1,k
			read(93,*) sample_site(i,:)
			!write(92,'(30(I5,2x))') sample_site(i,:)
		enddo
	end subroutine read_sites

	subroutine init_sites_hom()
		integer :: j
		do j=1,3
			sample_site(1,j) = 1
		enddo
	end subroutine init_sites_hom

	subroutine init_sample_zeta()
		sample_zeta(:) = 0.0_dl
		sample_zeta_f(:,:) = 0.0_dl
		sample_dzeta(:,:) = 0.0_dl
		sample_dzeta_f(:,:,:) = 0.0_dl
	end subroutine init_sample_zeta

	subroutine get_sample_field()
		integer :: i
		integer :: j		
		integer :: k
		if (run_hom==.FALSE.) then
			k = n_sample
		elseif (run_hom==.TRUE.) then
			k = 1
		endif
		do i=1,nfld
			!do j=1,n_sample
			do j=1,k
				sample_field(i,j) = fld(i,sample_site(j,1),sample_site(j,2),sample_site(j,3))
				sample_fieldp(i,j) = fldp(i,sample_site(j,1),sample_site(j,2),sample_site(j,3))
			enddo			
		enddo
	end subroutine get_sample_field

	subroutine write_lat_sample(time)
		real(dl) :: time
		integer :: i
		integer :: k
		if (run_hom==.FALSE.) then
			k = n_sample
		elseif (run_hom==.TRUE.) then
			k = 1
		endif
		call get_sample_field()
		!do i=1,n_sample
		do i=1,k
			write(94,'(30(ES22.15,2x))') time, sample_field(:,i), sample_fieldp(:,i), sample_zeta(i), sample_zeta_f(:,i)
		enddo
		write(94,*) 	
	end subroutine write_lat_sample

end module sample_sites_mod

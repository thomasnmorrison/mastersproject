!!!! zeta.f90 !!!!
!!!! Compute the zeta parameter
!!!! 

! to do: make better use of arrays so that set_fld_temp can be generalized to nfld

! to do: test zeta_lat calculation
! to do: output equation of state

module zeta_mod

!	use, intrinsic :: iso_c_binding
!  include 'fftw3.f03'
!	use fftw3
!	use params
#include "macros.h"
	use hamiltonian
	use sample_sites_mod
	use moments_mod

	real(dl) :: zeta = 0.0_dl
	real(dl), dimension(2) :: dzeta = (/0.0_dl,0.0_dl/)
	real(dl), dimension(IRANGE) :: zeta_lat
	real(dl), dimension(2,IRANGE) :: dzeta_lat !dzeta/dtau
	!real(dl), dimension(2,IRANGE) :: dzeta_lat_temp	! 1st component is numerator of dzeta, 2nd component is denominator of dzeta
	
	! additional pointers needed for calculating dzeta
	real(C_DOUBLE), pointer :: ggrad1(:,:,:)
	real(C_DOUBLE), pointer :: ggrad2(:,:,:)
	real(C_DOUBLE), pointer :: ggrad3(:,:,:)
	real(C_DOUBLE), pointer :: ggrad4(:,:,:)
	real(C_DOUBLE), pointer :: zf1(:,:,:)
	real(C_DOUBLE), pointer :: zf2(:,:,:)
	real(C_DOUBLE), pointer :: zf3(:,:,:)!new
	real(C_DOUBLE), pointer :: zf4(:,:,:)!new
	real(C_DOUBLE), pointer :: zfp1(:,:,:)
	real(C_DOUBLE), pointer :: zfp2(:,:,:)
	real(C_DOUBLE), pointer :: zlap1(:,:,:)
	real(C_DOUBLE), pointer :: zlap2(:,:,:)
	complex(C_DOUBLE_COMPLEX), pointer :: Fk3(:,:,:)

contains
	!!!! Step zeta integration !!!!
!	subroutine zeta_step(zeta, dzeta, dt)
!		real(dl), dimension(2) :: dzeta
!		real(dl) :: dt
!		
!		zeta = zeta + 0.5_dl*sum(dzeta)*dt	!!
!	end subroutine zeta_step


	subroutine zeta_init()
		zeta = 0.0_dl
		dzeta(:) = 0.0_dl
		zeta_lat(IRANGE) = 0.0_dl
		dzeta_lat(:,IRANGE) = 0.0_dl
	end subroutine zeta_init

	!!!! Step zeta integration using trapazoid rule
	subroutine zeta_step(dt, nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		call get_dzeta(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		zeta = zeta + 0.5_dl*sum(dzeta)*dt
		sample_zeta(:) = sample_zeta(:) + 0.5_dl*(sample_dzeta(1,:)+sample_dzeta(2,:))*dt
		!sample_zeta_f(:,:) = sample_zeta_f(:,:) + 0.5_dl*(sample_dzeta_f(1,:,:)+sample_dzeta_f(2,:,:))*dt
		zeta_lat(IRANGE) = zeta_lat(IRANGE) + 0.5_dl*(dzeta_lat(1,IRANGE)+dzeta_lat(2,IRANGE))*dt
	end subroutine zeta_step

	subroutine set_fld_temp()
		if (nfld == 2) then
			zf1 = fld(1,IRANGE)
			zf2 = fld(2,IRANGE)
			zfp1 = fldp(1,IRANGE)
			zfp2 = fldp(2,IRANGE)
		elseif (nfld ==1) then
			zf1 = fld(1,IRANGE)
			zf2 = 0.0_dl*fld(1,IRANGE)
			!zf2 = fld(1,IRANGE)!THIS IS WRONG ONLY FOR TESTING
			zfp1 = fldp(1,IRANGE)
			zfp2 = 0.0_dl*fldp(1,IRANGE)
			!zfp2 = fldp(1,IRANGE)!THIS IS WRONG ONLY FOR TESTING
		endif
	end subroutine set_fld_temp

	! subroutint to set temp fields in get_dzeta_test
#ifdef GHOST
	subroutine set_fld_temp_test(ind)
		integer :: ind
		
		if (ind <= nfld) then
			zf1 = fld(ind,IRANGE)
			zfp1 = fldp(ind,IRANGE)
		elseif (ind > nfld) then
			zf1 = ghst(ind-nfld,IRANGE)
			zfp1 = ghstp(ind-nfld,IRANGE)
	end subroutine set_fld_temp_test
#endif

	!rename ggrad to gg in order to not overload the namespace
!	subroutine get_dzeta(nsize, dk, dzeta, yscl, ysclp, f1, f2, fp1, fp2, gg1, gg2, lap1,lap2, gg3, gg4, Fk1, planf, planb)
!	subroutine get_dzeta(nsize, dk, yscl, ysclp, f1, f2, fp1, fp2, gg1, gg2, lap1,lap2, gg3, gg4, Fk1, planf, planb)
!	subroutine get_dzeta(nsize, dk, yscl, ysclp, f1, f2, fp1, fp2, gg1, gg2, lap1,lap2, gg3, gg4, Fk1, planf, planb)
!		integer, dimension(1:3), intent(in) :: nsize
!		real*8 :: dk	
!		real(dl), dimension(2) :: dzeta	!making this an array would allow for using a higher order integration scheme
!		real(dl) :: yscl
!		real(C_DOUBLE), pointer :: f1(:,:,:), f2(:,:,:), fp1(:,:,:), fp2(:,:,:), gg1(:,:,:), gg2(:,:,:), lap1(:,:,:), lap2(:,:,:), gg3(:,:,:), gg4(:,:,:)
!    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:)
!		type(C_PTR) :: planf, planb
!
!		real(dl) :: hub
!		integer :: n1, n2, n3
!
!		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
!		hub = !whatever in terms of yscl and ysclp
!
!		dzeta(1) = dzeta(2)	!store old value
!		
!		lap1 = fp1
!		lap2 = f1
!		call graddotgrad(nsize,lap1,lap2,gg1,Fk1,Fk2,Fk3, planf, planb)
!		lap1 = fp2
!		lap2 = f2
!		call graddotgrad(nsize,lap1,lap2,gg2,Fk1,Fk2,Fk3, planf, planb)
!		lap1 = f1
!		lap2 = f2
!		call laplacian_spectral(n1,n2,n3,lap1,Fk1,dk,planf,planb)
!		call laplacian_spectral(n1,n2,n3,lap2,Fk1,dk,planf,planb)
!		call gradient_squared_3d_spectral([n1,n2,n3],f1,Fk1,Fk2,gg3,dk,planf,planb)
!		call gradient_squared_3d_spectral([n1,n2,n3],f2,Fk1,Fk2,gg4,dk,planf,planb)
!
!		f1 = (yscl*(gg1 + gg2 + fp1*lap1 + fp2*lap2) + hub*yscl**4*(gg3 + gg4)) / (3.0_dl*fp1**2 + 3.0_dl*fp2**2 + yscl**4*(gg3 + gg4))	
!
!		dzeta(2) = sum(f1(IRANGE))/dble(n1)/dble(n2)/dble(n3)
!
!		dzeta(2) = dzeta(2) +
!		dzeta(2) = yscl**2*dzeta(2) /3.0_dl
!
!	end subroutine get_dzeta

	!!!! Updates new value of dzeta into dzeta(2), stores old value of dzeta into dzeta(1)
	subroutine get_dzeta(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		real(dl) :: hub
		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		!hub = -ysclp/(6.0_dl*yscl**2)/dble(n1)/dble(n2)/dble(n3) !!!!Double check that this should not be dividing number of lattice sites
		hub = -ysclp/(6.0_dl*yscl**2)

		dzeta(1) = dzeta(2)	!store old value
		dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
		sample_dzeta(1,:) = sample_dzeta(2,:)
		!sample_dzeta_f(1,:,:) = sample_dzeta_f(2,:,:)

		call set_fld_temp()
		zlap1 = zfp1
		zlap2 = zf1
		call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
		!if (nfld==2) then
			zlap1 = zfp2
			zlap2 = zf2
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad2,Fk1,Fk2,Fk3, planf, planb)
		!endif
		zlap1 = zf1
		call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)
		!if (nfld==2) then
			zlap2 = zf2
			call laplacian_spectral(n1,n2,n3,zlap2,Fk1,dk,planf,planb)
		!endif
		zf3 = zf1
		call graddotgrad(nsize,dk,zf1,zf3,ggrad3,Fk1,Fk2,Fk3, planf, planb)
		!if (nfld==2) then
			zf4 = zf2
			call graddotgrad(nsize,dk,zf2,zf4,ggrad4,Fk1,Fk2,Fk3, planf, planb)
		!endif		
		!!!! Use graddotgrad instead of gradient_squared_3d_spectral, to be on the same side should use a field copy for one input
		!call gradient_squared_3d_spectral([n1,n2,n3],zf1,Fk1,Fk2,ggrad3,dk,planf,planb)
		!call gradient_squared_3d_spectral([n1,n2,n3],zf2,Fk1,Fk2,ggrad4,dk,planf,planb)

		!zf1 = (yscl*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2) + hub*yscl**4*(ggrad3 + ggrad4)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))!!!Included extra term
		!if (nfld==2) then
			zf1 = (yscl**2*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))!!!Without extra term and squaring yscl to account for conformal time
		!endif
		!if (nfld==1) then
		!	zf1 = (yscl**2*(ggrad1 + zfp1*zlap1)) / (3.0_dl*zfp1**2 + yscl**4*(ggrad3))!!!Without extra term and squaring yscl to account for conformal time
		!endif
		!zf1 = (yscl*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2)

		
		do i=1,n_sample
			!sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,1),sample_site(i,1))/dble(n1)/dble(n2)/dble(n3)!!!!!!!!!!!!SHOULD NOT DIVIDE BY LATTICE SITES!!!!!!!!!!
			!sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,1),sample_site(i,1))!looks like I'm choosing the wrong site
			sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,2),sample_site(i,3))
		enddo
		dzeta(2) = sum(zf1(IRANGE))/dble(n1)/dble(n2)/dble(n3)
		dzeta_lat(2,IRANGE) = zf1(IRANGE)

	!testing
		do i=1,nfld
			if (i==1) then
				zf2 = (yscl**2*(ggrad1 + zfp1*zlap1)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))
			elseif (i==2) then
				zf2 = (yscl**2*(ggrad2 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))
			endif
			do j=1,n_sample
				!sample_dzeta_f(2,i,j) = zf2(sample_site(j,1),sample_site(j,2),sample_site(j,3))
			enddo	
		enddo
		!print*, ggrad1(6,6,6)
		!print*, zfp1(6,6,6)
		!print*, zlap1(6,6,6)
		!print*, ggrad3(6,6,6)
		!print*, zf1(6,6,6)
	end subroutine get_dzeta

	subroutine get_dzeta_test(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		real(dl) :: hub
		integer :: n1, n2, n3
		integer :: i,j

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		hub = -ysclp/(6.0_dl*yscl**2)

		dzeta(1) = dzeta(2)	!store old value
		dzeta_lat(1,IRANGE) = dzeta_lat(2,IRANGE)
		sample_dzeta(1,:) = sample_dzeta(2,:)
		!sample_dzeta_f(1,:,:) = sample_dzeta_f(2,:,:)

!#ifdef GHOST
	!do n=1,2*nfld
!#elif
	!do n=1,nfld
!#endif
	!enddo
		call set_fld_temp()
		zlap1 = zfp1
		zlap2 = zf1
		call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
		!if (nfld==2) then
			zlap1 = zfp2
			zlap2 = zf2
			call graddotgrad(nsize,dk,zlap1,zlap2,ggrad2,Fk1,Fk2,Fk3, planf, planb)
		!endif
		zlap1 = zf1
		call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)
		!if (nfld==2) then
			zlap2 = zf2
			call laplacian_spectral(n1,n2,n3,zlap2,Fk1,dk,planf,planb)
		!endif
		zf3 = zf1
		call graddotgrad(nsize,dk,zf1,zf3,ggrad3,Fk1,Fk2,Fk3, planf, planb)
		!if (nfld==2) then
			zf4 = zf2
			call graddotgrad(nsize,dk,zf2,zf4,ggrad4,Fk1,Fk2,Fk3, planf, planb)
		!endif		

		zf1 = (yscl**2*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))!!!Without extra term and squaring yscl to account for conformal time

		do i=1,n_sample
			sample_dzeta(2,i) = zf1(sample_site(i,1),sample_site(i,2),sample_site(i,3))
		enddo
		dzeta(2) = sum(zf1(IRANGE))/dble(n1)/dble(n2)/dble(n3)
		dzeta_lat(2,IRANGE) = zf1(IRANGE)

		!do i=1,nfld
		!	if (i==1) then
		!		zf2 = (yscl**2*(ggrad1 + zfp1*zlap1)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))
		!	elseif (i==2) then
		!		zf2 = (yscl**2*(ggrad2 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))
		!	endif
		!	do j=1,n_sample
				!sample_dzeta_f(2,i,j) = zf2(sample_site(j,1),sample_site(j,2),sample_site(j,3))
		!	enddo	
		!enddo
	end subroutine get_dzeta_test


! calculates grad(f1).grad(f2) and outputs in gg (f1 and f2 are also overwritten)
! check if I can do this by updating a partial sum
	subroutine graddotgrad(nsize, dk, f1, f2, gg, Fk1, Fk2, Fk3, planf, planb)
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk
		real(C_DOUBLE), pointer :: f1(:,:,:), f2(:,:,:), gg(:,:,:)!, g1(:,:,:), g2(:,:,:), g3(:,:,:)
		complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:), Fk3(:,:,:)!, Fk(:,:,:)
		type(C_PTR) :: planf, planb

		integer :: n1, n2, n3
		integer :: nn1, nn2, nn3
    integer :: i,j,k,ii,jj,kk
		!real(C_DOUBLE) :: kcur!new
    complex(C_DOUBLE_COMPLEX), parameter :: i_imag = (0.,1.)!new
		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1

		call fftw_execute_dft_r2c(planf, f1, Fk1)
		call fftw_execute_dft_r2c(planf, f2, Fk2)

		!do k=1,n3; if (k>nn3) then; kk=nn3+1-k;

		!!!!!!!! putting in i_imag instead of factoring it out, to test if the c2r fft handles things properly if i_imag is factored out
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble((i-1)*dk)*Fk1(i,j,k)
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble((i-1)*dk)*Fk2(i,j,k)
					Fk3(i,j,k) = i_imag*dble((i-1)*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		!gg = -f1*f2
		gg = f1*f2 !using i_imag

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(jj*dk)*Fk1(i,j,k)
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(jj*dk)*Fk2(i,j,k)
					Fk3(i,j,k) = i_imag*dble(jj*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)
		!gg = gg-f1*f2
		gg = gg+f1*f2 !using i_imag
		
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(kk*dk)*Fk1(i,j,k)
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					!Fk3(i,j,k) = dble(kk*dk)*Fk2(i,j,k)
					Fk3(i,j,k) = i_imag*dble(kk*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)
		!gg = gg-f1*f2
		gg = gg+f1*f2 !using i_imag
		gg = gg /dble(n1)**2 /dble(n2)**2 /dble(n3)**2

!		!!!!Alrternative calculation with fewer lines, more variables!!!!
!		call fftw_execute_dft_r2c(planf, f1, Fk)
!		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
!			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
!				do i=1,nn1
					! this is now a vector so I will need separate DFTs
!					Fk1(i,j,k) = dble((i-1)*dk)*Fk(i,j,k)
!					Fk2(i,j,k) = dble(jj*dk)*Fk(i,j,k)
!					Fk3(i,j,k) = dble(kk*dk)*Fk(i,j,k)					
!				enddo
!			enddo
!		enddo
!		call fftw_execute_dft_c2r(planb, Fk1, g1)
!		call fftw_execute_dft_c2r(planb, Fk2, g2)
!		call fftw_execute_dft_c2r(planb, Fk3, g3)
!		
!		call fftw_execute_dft_r2c(planf, f2, Fk)
!		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
!			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
!				do i=1,nn1
!					! this is now a vector so I will need separate DFTs
!					Fk1(i,j,k) = dble((i-1)*dk)*Fk(i,j,k)
!					Fk2(i,j,k) = dble(jj*dk)*Fk(i,j,k)
!					Fk2(i,j,k) = dble(kk*dk)*Fk(i,j,k)					
!				enddo
!			enddo
!		enddo
!		call fftw_execute_dft_c2r(planb, Fk1, f1)
!		call fftw_execute_dft_c2r(planb, Fk2, f2)
!		call fftw_execute_dft_c2r(planb, Fk3, gg)

!		f1(:,:,:) = -f1(:,:,:)*g1(:,:,:) - f2(:,:,:)*g2(:,:,:) - gg(:,:,:)*g3(:,:,:)
!		f1 = f1 /dble(n1) /dble(n2) /dble(n3)

	end subroutine graddotgrad

! function to calculate averaged characteristic function
!	function char_func(f1)
!		real(dl) :: char_func
!		real(dl) :: f1(:,:,:)
!		real(dl) :: char_lat(:,:,:)
!		real(dl) :: f1_mean
!		integer, dimension(1:3) :: n_size = (/nx,ny,nz/)
		
!		f1_mean = get_mean_3d(f1, n_size)
!		char_lat(IRANGE) = dexp(-log(yscl)*(f1(IRANGE)-f1_mean))
!		char_func = sum(char_lat(IRANGE))/dble(nx)/dble(ny)/dble(nz)
!!!!! complete this function
!	end function char_func

!!! Add characteristic function to output
	subroutine write_zeta(time)
		real(dl) :: time
		real(dl) :: zeta_moment1, zeta_moment2, zeta_moment3, zeta_moment4!, zeta_char
		integer, dimension(1:3) :: n_size = (/nx,ny,nz/)

		!n_size(1) = nx; n_size(2) = ny; n_size(3) = nz
		
		zeta_moment1 = get_mean_3d(zeta_lat, n_size)
		zeta_moment2 = get_moment_3d(zeta_lat, zeta_moment1, 2, n_size)
		zeta_moment3 = get_moment_3d(zeta_lat, zeta_moment1, 3, n_size)
		zeta_moment4 = get_moment_3d(zeta_lat, zeta_moment1, 4, n_size)
		!zeta_char = char_func(zeta_lat)
		write(96,'(30(ES22.15,2x))') time, dzeta(1), dzeta(2), zeta, zeta_moment1, zeta_moment2, zeta_moment3, zeta_moment4!, zeta_char
		!write(96,*)

	end subroutine write_zeta

end module zeta_mod







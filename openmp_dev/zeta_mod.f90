!!!! zeta.f90 !!!!
!!!! Compute the zeta parameter
!!!! 

! to do: figure out the rho + P averaging, not required during inflation since there are no oscillations
! to do: make better use of arrays so that set_fld_temp can be generalized to nfld

module zeta_mod

!	use, intrinsic :: iso_c_binding
!  include 'fftw3.f03'
!	use fftw3
!	use params
#include "macros.h"
	use hamiltonian

	real(dl) :: zeta = 0.0_dl
	real(dl), dimension(2) :: dzeta = (/0.0_dl,0.0_dl/)
	
	! additional pointers needed for calculating dzeta
	real(C_DOUBLE), pointer :: ggrad1(:,:,:)
	real(C_DOUBLE), pointer :: ggrad2(:,:,:)
	real(C_DOUBLE), pointer :: ggrad3(:,:,:)
	real(C_DOUBLE), pointer :: ggrad4(:,:,:)
	real(C_DOUBLE), pointer :: zf1(:,:,:)
	real(C_DOUBLE), pointer :: zf2(:,:,:)
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

	!!!! Step zeta integration using trapazoid rule
	subroutine zeta_step(dt, nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		real(dl) :: dt
		integer, dimension(1:3), intent(in) :: nsize
		real*8 :: dk	
    complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:),Fk2(:,:,:),Fk3(:,:,:)
		type(C_PTR) :: planf, planb

		call get_dzeta(nsize, dk, Fk1, Fk2, Fk3, planf, planb)
		zeta = zeta + 0.5_dl*sum(dzeta)*dt
	end subroutine zeta_step

	subroutine set_fld_temp()
		zf1 = fld(1,IRANGE)
		zf2 = fld(2,IRANGE)
		zfp1 = fldp(1,IRANGE)
		zfp2 = fldp(2,IRANGE)
	end subroutine set_fld_temp

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

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		hub = -ysclp/(6.0_dl*yscl**2)/dble(n1)/dble(n2)/dble(n3)

		dzeta(1) = dzeta(2)	!store old value
		call set_fld_temp()
		zlap1 = zfp1
		zlap2 = zf1
		call graddotgrad(nsize,dk,zlap1,zlap2,ggrad1,Fk1,Fk2,Fk3, planf, planb)
		zlap1 = zfp2
		zlap2 = zf2
		call graddotgrad(nsize,dk,zlap1,zlap2,ggrad2,Fk1,Fk2,Fk3, planf, planb)
		zlap1 = zf1
		zlap2 = zf2
		call laplacian_spectral(n1,n2,n3,zlap1,Fk1,dk,planf,planb)
		call laplacian_spectral(n1,n2,n3,zlap2,Fk1,dk,planf,planb)
		call gradient_squared_3d_spectral([n1,n2,n3],zf1,Fk1,Fk2,ggrad3,dk,planf,planb)
		call gradient_squared_3d_spectral([n1,n2,n3],zf2,Fk1,Fk2,ggrad4,dk,planf,planb)

		zf1 = (yscl*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2) + hub*yscl**4*(ggrad3 + ggrad4)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2 + yscl**4*(ggrad3 + ggrad4))	

		!zf1 = (yscl*(ggrad1 + ggrad2 + zfp1*zlap1 + zfp2*zlap2)) / (3.0_dl*zfp1**2 + 3.0_dl*zfp2**2)
		
		dzeta(2) = sum(zf1(IRANGE))/dble(n1)/dble(n2)/dble(n3)

	end subroutine get_dzeta

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

		n1 = nsize(1); n2 = nsize(2); n3 = nsize(3)
		nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1

		call fftw_execute_dft_r2c(planf, f1, Fk1)
		call fftw_execute_dft_r2c(planf, f2, Fk2)

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = dble((i-1)*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = dble((i-1)*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)

		gg = -f1*f2

		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = dble(jj*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = dble(jj*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)
		gg = gg-f1*f2

		
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = dble(kk*dk)*Fk1(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f1)
		do k=1,n3; if(k>nn3) then; kk = k-n3-1; else; kk=k-1; endif
			do j=1,n2; if(j>nn2) then; jj = j-n2-1; else; jj=j-1; endif
				do i=1,nn1
					Fk3(i,j,k) = dble(kk*dk)*Fk2(i,j,k)
				enddo
			enddo
		enddo
		call fftw_execute_dft_c2r(planb, Fk3, f2)
		gg = gg-f1*f2
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

	subroutine write_zeta(time)
		real(dl) :: time

		write(96,'(30(ES22.15,2x))') time, dzeta(1), dzeta(2), zeta
		!write(96,*)

	end subroutine write_zeta

end module zeta_mod







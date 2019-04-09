!!!! correlator.f90 !!!!
! Module provides subroutines that will generate an initial conditions power spectrum
! the power spectrum generated can be picked up by init_fields_sr to generate initial
! conditions for the fields.

! to do: write a subroutine that will generate power spectrum with correct horizon crossing

module correlator_mod

#include "macros.h"
	use params
	use analysis
	use potential_mod

integer, parameter :: os = 16 
integer, parameter :: nos = max(nx,ny,nz)*os**2
real(dl), parameter :: dxos = dx/os							! oversampled grid spacing
real(dl), parameter :: dkos = dk/(2*os)					! oversampled mode spacing

!complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld,nos) :: cor
real(dl), dimension(2*nfld,2*nfld,nos) :: cor
real(dl), dimension(2*nfld,2*nfld,nos) :: cor_rad					! radial profile
!real(dl), dimension(2*nfld,2*nfld,nnx,ny,nz) :: cor_lat		! all modes on lattice
complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld,nnx,ny,nz) :: cor_lat		! all modes on lattice
real(dl), dimension(nfld) :: m2_diag

contains

	! initialize the m2_diag vector
	! assume constant masses
	! assume fields are eigen vectors of m2 matrix
	! to do: put in cases for different forms of potential
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

	! Generate power spectrum using a step function to go from Minkowski space fluctuations
	! subhorizon to on the slow roll attractor superhorizon with Hankle function solutions. 
	! Assumes a_0=1
	! to do: take as input the initial fields
	subroutine slow_roll_cor(hub)
		real(dl) :: hub			! input initial Hubble parameter
		integer :: k,n,m		! 
		real(dl) :: kk			! comoving wavenumber (double check the units)
		
		call init_m2_diag()
		! initalize to zero power
		!cor(:,:,:) = (0._dl,0._dl)
		cor(:,:,:) = 0._dl
		! Loop over sampled modes and calculate power
		do k=1,nos; kk = (dble(k)-0.5_dl)*dkos
			! select super/sub horizon
			if (kk <= hub) then
				do n=1,nfld
					! field-field
					cor(2*(n-1),2*(n-1),k) = 0.5_dl*(1._dl + (hub/kk)**2)/(kk*mpl**2)
					! field-momentum
					cor(2*(n-1),2*n,k) = -0.5_dl*(kk**2+m2_diag(n))/(3._dl*hub)*(1._dl + (hub/kk)**2)/(kk*mpl**2)
					cor(2*n,2*(n-1),k) = cor(2*(n-1),2*n,k)
					! momentum-momentum
					cor(2*n,2*n,k) = 0.5_dl*((kk**2+m2_diag(n))/(3._dl*hub))**2*(1._dl + (hub/kk)**2)/(kk*mpl**2)
				enddo
			else
				do n=1,nfld	
					! field-field
					cor(2*(n-1),2*(n-1),k) = 0.5_dl*(kk**2+m2_diag(n))**-0.5/(mpl**2)
					! momentum-momentum
					cor(2*n,2*n,k) = 2.0_dl*(kk**2+m2_diag(n))**0.5/(mpl**2)
				enddo
			end if
		enddo

		!testing
		cor(:,:,:) = 0._dl
		do k=1,nos;
			do n=1,2*nfld; do m=1,2*nfld
				if (n==m) then
					cor(n,m,k) = 1._dl
				endif
			enddo
			print*, cor(n,:,k)
			enddo
		enddo
	
	end subroutine slow_roll_cor

	! Generate power spectrum using a step function to go from Minkowski space fluctuations
	! subhorizon to on the slow roll attractor superhorizon with Hankle function solutions. 
	! Assumes a_0=1
	! to do: take as input the initial fields
	! to do: check that factor of dk is correct for rad2 <= hub**2
	! to do: check that only using modes i=1,nnx is sufficient
	! to do: use an if k /=0 statement to not divide by 0
	! to do: make cor_lat Hermetian on field indicies
	! to do: density of states term
	subroutine slow_roll_cor_lat(hub)
		real(dl) :: hub			! input initial Hubble parameter
		integer :: i,j,k		!	lattice indicies
		integer :: ii,jj,kk	!
		integer :: m,n			! field indicies 
		real(dl) :: rad2		! comoving wavenumber (double check the units)
		
		call init_m2_diag()
		! initalize to zero power
		cor_lat(:,:,:,:,:) = 0._dl*(1._dl,0._dl)
		! Loop over all modes on the lattice and set power
		do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
		do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			do i=1,nnx
				rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
				rad2 = rad2*dk**2
				if (rad2 == 0._dl) then
					do n=1,nfld
						! field-field
						cor_lat(2*n-1,2*n-1,LATIND) = 1._dl*(1._dl,0._dl)
						! field-momentum
						cor_lat(2*n-1,2*n,LATIND) = 0._dl*(1._dl,0._dl)
						cor_lat(2*n,2*n-1,LATIND) = cor_lat(2*(n-1),2*n,LATIND)
						! momentum-momentum
						cor_lat(2*n,2*n,LATIND) = 1._dl*(1._dl,0._dl)
					enddo
				elseif (rad2 <= hub**2) then
					do n=1,nfld
						! field-field
						!cor_lat(2*n-1,2*n-1,LATIND) = 0.5_dl*(1._dl + hub**2/rad2)/(sqrt(rad2)*mpl**2)
						cor_lat(2*n-1,2*n-1,LATIND) = (1._dl,0._dl)/sqrt(rad2)**3	! TESTING
						! field-momentum
						!cor_lat(2*n-1,2*n,LATIND) = -0.5_dl*(rad2+m2_diag(n))/(3._dl*hub)*(1._dl + hub**2/rad2)/(sqrt(rad2)*mpl**2)
						!cor_lat(2*n,2*n-1,LATIND) = cor_lat(2*n-1,2*n,LATIND)
						cor_lat(2*n-1,2*n,LATIND) = 0.5_dl*(1._dl,0._dl)/sqrt(rad2)**3	! TESTING
						cor_lat(2*n,2*n-1,LATIND) = 0.5_dl*(1._dl,0._dl)/sqrt(rad2)**3	! TESTING
						! momentum-momentum
						!cor_lat(2*n,2*n,LATIND) = 0.5_dl*((rad2+m2_diag(n))/(3._dl*hub))**2*(1._dl + hub**2/rad2)/(sqrt(rad2)*mpl**2)
						cor_lat(2*n,2*n,LATIND) = (1._dl,0._dl)/sqrt(rad2)**3	! TESTING
					enddo
				else
					do n=1,nfld	
						! field-field
						!cor_lat(2*n-1,2*n-1,LATIND) = 0.5_dl*(rad2+m2_diag(n))**-0.5/(mpl**2)
						cor_lat(2*n-1,2*n-1,LATIND) = (1._dl,0._dl)/sqrt(rad2)	! TESTING
						! momentum-momentum
						!cor_lat(2*n,2*n,LATIND) = 2.0_dl*(rad2+m2_diag(n))**0.5/(mpl**2)
						cor_lat(2*n,2*n,LATIND) = sqrt(rad2)*(1._dl,0._dl)	! TESTING
					enddo
				end if				
			enddo
		enddo
		enddo
	
	end subroutine slow_roll_cor_lat

	! Subroutine to initialize a correlation matrix for the Minkowski space vacuum fluctuations
	! n.b. the correlation matrix that is initialized here will have the same value as the initial
	! spectrum, ie. all normalizations must match between the two.
	subroutine mink_cor_lat()
		integer :: i,j,k		!	lattice indicies
		integer :: ii,jj,kk	! fourier modes
		integer :: m,n			! field indicies 
		real(dl) :: rad2		! comoving wavenumber (double check the units)

		cor_lat(:,:,:,:,:) = 0._dl*(1._dl,0._dl)	! initialize zero power in all modes
		call init_m2_diag()												! initialize mass matrix

		! Loop over all modes on the lattice and set power
		do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
		do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			do i=1,nnx
				rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
				rad2 = rad2*dk**2
				do n=1,nfld
					cor_lat(2*n-1,2*n-1,LATIND) = (1._dl,0._dl)/(2._dl*sqrt(rad2 + m2_diag(n))*mpl**2)	! field-field
					cor_lat(2*n,2*n,LATIND) = (1._dl,0._dl)*sqrt(rad2 + m2_diag(n))/(2._dl*mpl**2)			! momentum-momentum
				enddo
			enddo
		enddo
		enddo

		cor_lat = nvol**2*cor_lat	! nvol from FFT convension (redundant with an nvol divison when initializing fields)
	end subroutine mink_cor_lat

	! Subroutine to intitialize a correlation matrix for the Minkowski space vacuum fluctuations squeezed in the pi
	! direction or phi direction.
	subroutine mink_cor_lat_scale(fld_scl, dfld_scl)
		real(dl) :: fld_scl		! factor by which field fluctuation power is scaled
		real(dl) :: dfld_scl	! factor by which momentum fluctuation power is scaled
		integer :: i,j,k			!	lattice indicies
		integer :: ii,jj,kk		! fourier modes
		integer :: m,n				! field indicies 
		real(dl) :: rad2			! comoving wavenumber (double check the units)

		cor_lat(:,:,:,:,:) = 0._dl*(1._dl,0._dl)	! initialize zero power in all modes
		call init_m2_diag()												! initialize mass matrix

		! Loop over all modes on the lattice and set power
		do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
		do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			do i=1,nnx
				rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
				rad2 = rad2*dk**2
				do n=1,nfld
					cor_lat(2*n-1,2*n-1,LATIND) = fld_scl*(1._dl,0._dl)/(2._dl*sqrt(rad2 + m2_diag(n))*mpl**2)	! field-field
					cor_lat(2*n,2*n,LATIND) = dfld_scl*(1._dl,0._dl)*sqrt(rad2 + m2_diag(n))/(2._dl*mpl**2)			! momentum-momentum
				enddo
			enddo
		enddo
		enddo

		cor_lat = nvol**2*cor_lat	! nvol from FFT convension (redundant with an nvol divison when initializing fields)
	end subroutine mink_cor_lat_scale

	! Subroutine to perform a linear transformation on the spectrum of Minkowski space vacuum fluctuations
	! built for one field model
	! to do: declare transformation matrix
	! to do: construct transformation matrix from scl, ar, and rot
	! to do: transform cor_lat using transformation matrix
	! to do check row/column on tran_mat
	subroutine mink_cor_lat_transform(scl, ar, rot)
		real(dl) :: scl		! Overall scal pactor for P_tensor
		real(dl) :: ar		! Transformation to the aspect ratio of P_phiphi to P_pipi
		real(dl) :: rot		! Rotation angle of P_tensor, generates phi-pi correlation
		integer :: i,j,k			!	lattice indicies
		integer :: ii,jj,kk		! fourier modes
		integer :: m,n				! field indicies 
		real(dl) :: rad2			! comoving wavenumber (double check the units)
		complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: tran_mat	! transformation matrix

		! Construct transfromation matrix
		tran_mat(:,:) = 0._dl*(1._dl,0._dl)
		do n=1,nfld
			tran_mat(2*n-1,2*n-1) = scl*cos(rot)*ar*(1._dl,0._dl)
			tran_mat(2*n-1,2*n) = -scl*sin(rot)/ar*(1._dl,0._dl)
			tran_mat(2*n,2*n-1) = scl*sin(rot)*ar*(1._dl,0._dl)
			tran_mat(2*n,2*n) = scl*cos(rot)/ar*(1._dl,0._dl)
		enddo

	end subroutine mink_cor_lat_transform

	! Subroutine to generate power spectrum by scaling Minkowski space spectrum with Det(P) = const
	! n.b. fld_scl*dfld_scl 
	subroutine mink_cor_lat_scale2(fld_scl, dfld_scl)
		real(dl) :: fld_scl		! scaling factor for <fld**2>
		real(dl) :: dfld_scl	! scaling factor for <dfld**2>
		integer :: i,j,k			!	lattice indicies
		integer :: ii,jj,kk		! fourier modes
		integer :: m,n				! field indicies 
		real(dl) :: rad2			! comoving wavenumber (double check the units)
		
		cor_lat(:,:,:,:,:) = 0._dl*(1._dl,0._dl)	! initialize zero power in all modes
		call init_m2_diag()												! initialize mass matrix

		! Loop over all modes on the lattice and set power
		do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
		do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			do i=1,nnx
				rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
				rad2 = rad2*dk**2
				do n=1,nfld
					cor_lat(2*n-1,2*n-1,LATIND) = fld_scl*(1._dl,0._dl)/(2._dl*sqrt(rad2 + m2_diag(n))*mpl**2)	! field-field
					cor_lat(2*n,2*n,LATIND) = dfld_scl*(1._dl,0._dl)*sqrt(rad2 + m2_diag(n))/(2._dl*mpl**2)			! momentum-momentum
					cor_lat(2*n-1,2*n,LATIND) = -(1._dl,0._dl)*sqrt(fld_scl*dfld_scl-1._dl)/(2._dl*mpl**2)				! field-momentum
					cor_lat(2*n,2*n-1,LATIND) = cor_lat(2*n-1,2*n,LATIND)
				enddo
			enddo
		enddo
		enddo

		cor_lat = nvol**2*cor_lat	! nvol from FFT convension (redundant with an nvol divison when initializing fields)
	end subroutine mink_cor_lat_scale2

	! Subroutine to generate power spectrum by scaling Minkowski space spectrum with Det(P) = const
	! n.b. don't use dfld_scl = 0
	subroutine mink_cor_lat_scale3(dfld_scl, cross_fac)
		real(dl) :: dfld_scl	! scaling factor for <dfld**2>
		real(dl) :: cross_fac	! factor to determin the cross spectrum of <fld*dfld>
		integer :: i,j,k			!	lattice indicies
		integer :: ii,jj,kk		! fourier modes
		integer :: m,n				! field indicies 
		real(dl) :: rad2			! comoving wavenumber (double check the units)
		
		cor_lat(:,:,:,:,:) = 0._dl*(1._dl,0._dl)	! initialize zero power in all modes
		call init_m2_diag()												! initialize mass matrix

		! Loop over all modes on the lattice and set power
		do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
		do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif
			do i=1,nnx
				rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
				rad2 = rad2*dk**2
				do n=1,nfld
					cor_lat(2*n-1,2*n-1,LATIND) = (4._dl*cross_fac**2+1._dl)/dfld_scl*(1._dl,0._dl)/(2._dl*sqrt(rad2 + m2_diag(n))*mpl**2)	! field-field
					cor_lat(2*n,2*n,LATIND) = dfld_scl*(1._dl,0._dl)*sqrt(rad2 + m2_diag(n))/(2._dl*mpl**2)			! momentum-momentum
					cor_lat(2*n-1,2*n,LATIND) = -(1._dl,0._dl)*cross_fac/(2._dl*mpl**2)				! field-momentum
					cor_lat(2*n,2*n-1,LATIND) = cor_lat(2*n-1,2*n,LATIND)
				enddo
			enddo
		enddo
		enddo

		cor_lat = nvol**2*cor_lat	! nvol from FFT convension (redundant with an nvol divison when initializing fields)
	end subroutine mink_cor_lat_scale3

	! Output initial power spectrum
	! to do: rewrite this subroutine for cor_lat
	subroutine write_cor()
		integer :: n
		integer :: k
		real(dl) :: kk
		
		do k=1,nos; kk = (k-0.5_dl)*dkos
			do n=1,2*nfld	
				write(91,'(30(ES22.15,2x))') kk, cor(n,:,k)
			enddo 
		enddo

		close(91)
	end subroutine write_cor

	! Output initial power spectrum
	! to do: calculate spectrum as in analysis_spec.f90
	! to do: the spectrum calculation needs modification due to the squaring
	subroutine write_cor_lat()
		real(dl), dimension(2*nfld,2*nfld,1:ns) :: S

    integer :: i,j,k,ii,jj,kk
    real(dl) :: p, c(2), W(2*nfld,2*nfld,ns)
    integer :: l
		integer :: m,n
		integer, parameter :: bin_rule = 0

    W = 0.
    S = 0.

		! This spectrum calulation is taken from analysis_spec.f90
		do m=1,2*nfld; do n=1,2*nfld
    do k=1,nz; if (k<=nnz) then; kk=k-1; else; kk=nz+1-k; endif
    do j=1,ny; if (j<=nny) then; jj=j-1; else; jj=ny+1-j; endif
    do i=1,nx; if (i<=nnx) then; ii=i-1; else; ii=nx+1-i; endif	! include negative momentum modes on x
			if (bin_rule==0) then
    		p = sqrt(dble(ii**2 + jj**2 + kk**2)); l=floor(p)
    		c = (1.0 - (/l-p,l+1-p/)**2)**2
				S(m,n,l+1:l+2) = S(m,n,l+1:l+2) + c*real(cor_lat(m,n,ii+1,j,k))
    		W(m,n,l+1:l+2) = W(m,n,l+1:l+2) + c
			! if using the floor routine
			else
				p = sqrt(dble(ii**2 + jj**2 + kk**2)); l=floor(p)
				S(m,n,l+1) = real(cor_lat(m,n,ii+1,j,k))
			endif
    enddo
    enddo
    enddo
		enddo; enddo  
    where (W /= 0.0) S = S/W!/dble(nx)/dble(ny)/dble(nz)

		! make output to match the formatting of spectrum.out
		if (nfld==2) then
		do i=1,ns
			write(91,'(30(ES22.15,2x))') (i-1)*dk, S(1,1,i), S(2,2,i), S(2,1,i), 0._dl, S(3,3,i), S(4,4,i), S(4,3,i), 0._dl, 0._dl, 0._dl, 0._dl, 0._dl, 0._dl, S(3,1,i), S(4,1,i), S(3,2,i), S(4,2,i)
		enddo
		endif

		close(91)
	end subroutine write_cor_lat

end module correlator_mod

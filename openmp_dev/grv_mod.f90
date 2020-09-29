! grv_mod.f90

module grv_mod

  use params

implicit none

	!integer, parameter :: dl = kind(1.d0)
	!real(dl), parameter :: twopi = 6.2831853071795864769252867665590

contains

	! Subroutine to initialize random number generator
	subroutine init_rng(seed_in)
		integer, intent(in) :: seed_in				! rng seed
		integer :: i
		integer :: nseed											! size of random rumber generated
		integer, allocatable :: seed(:)				

		call random_seed(SIZE=nseed)
		print*, 'nseed = ', nseed
    allocate(seed(nseed))
		seed = seed_in*(/ (i-1, i=1,nseed) /)
    call random_seed(PUT=seed)
    deallocate(seed)
	end subroutine init_rng

	! Function to return a complex Gaussian random variable with 0 mean and 1 variance
	function get_grv_complex(n) result(grv)
		integer :: n
		complex(dl), dimension(n) :: grv
		complex, parameter :: w = (0._dl, twopi)				! i*2*pi, used for random phase
		real(dl) :: a(n), p(n)													! amplitude, phase

		call random_number(a); call random_number(p)
		grv(:) = sqrt(-1._dl*log(a)) * exp(w*p)
	end function get_grv_complex

	! Function to return a real Gaussian random variable with 0 mean and 1 variance
	function get_grv_real(n) result(grv)
		integer, intent(in) :: n
		real(dl), dimension(n) :: grv
		real(dl) :: a(n), p(n)													! amplitude, phase

		call random_number(a); call random_number(p)		! a and p are dimension n
		grv(:) = (1._dl,0._dl)*sign(sqrt(-1._dl*log(a)),p-0.5_dl)
	end function get_grv_real

	! subroutine to test the variance of get_grv_real and get_grv_complex
	subroutine grv_test(n,n_max)
		integer :: n, n_max		
		integer :: i
		real(dl), dimension(n) :: grv_r
		complex(dl), dimension(n) :: grv_c
		real(dl), dimension(n) :: sum2_r
		real(dl), dimension(n) :: sum2_c
		real(dl), dimension(n) :: sum_r
		complex(dl), dimension(n) :: sum_c
	
		sum2_r = 0._dl; sum2_c=0._dl; sum_r=0._dl; sum_c=(0._dl,0._dl)
		do i=1, n_max
			grv_r = get_grv_real(n)
			grv_c = get_grv_complex(n)
			sum2_r = sum2_r + grv_r**2
			sum2_c = sum2_c + grv_c*conjg(grv_c)
			sum_r = sum_r + grv_r
			sum_c = sum_c + grv_c
		enddo
		sum2_r = sum2_r/dble(n_max); sum_r = sum_r/dble(n_max)
		sum2_c = sum2_c/dble(n_max); sum_c = sum_c/dble(n_max)
		print*, "real grv mean = ", sum_r
		print*, "complex grv mean = ", sum_c
		print*, "real grv var = ", sum2_r
		print*, "complex grv var = ", sum2_c
	end subroutine grv_test


end module grv_mod

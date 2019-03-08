module params
  integer, parameter :: dl = kind(1.d0)
  real(dl), parameter :: twopi = 6.2831853071795864769252867665590

	logical, parameter :: run_hom = .FALSE.!.TRUE.!  set to true for a homogeneous field calculation
	!integer, parameter :: nfld=2
	!integer, parameter :: os = 16	! Oversampling factor
  integer, parameter :: nx=64, ny=64, nz=64
  integer, parameter :: nnx=nx/2+1, nny=ny/2+1, nnz=nz/2+1
  integer, parameter :: nn = min(nnx,nny,nnz)

  real(dl), parameter :: nvol = dble(nx)*dble(ny)*dble(nz)
  real(dl), parameter :: len = 50._dl!10.!0.053_dl!0.18_dl !0.25_dl!10.!
  real(dl), parameter :: dx = len/dble(nx)
  real(dl), parameter :: dk = twopi / len

end module params

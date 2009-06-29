  !
  !----------------------------------------------------------------
  module parameters
  !
  integer, parameter :: DP = selected_real_kind(14,200)
  end module parameters
  !----------------------------------------------------------------
  program testpade
  !
  use parameters, only : DP
  implicit none
  integer, parameter :: N = 11
  !
  integer :: i, Nw
  real(DP) :: zr, zi, ur, ui
  complex(DP) :: a(N), z(N), u(N), w, padapp
  !
  do i = 1, N
    read(50,*) zr,zi, ur, ui
    z(i) = dcmplx(zr,zi)
    u(i) = dcmplx(ur,ui)
  enddo
  !
  call pade_coeff ( N, z, u, a)
  !
  Nw = 100 
  do i = 1, Nw
    w = dcmplx ( float(i-1)/float(Nw-1) * 50.d0, 0.d0 )
    call pade_eval ( N, z, a, w, padapp)
    write(100,*) aimag(w), real(w), real (padapp), aimag(padapp)
  enddo
  !
  end program testpade
  !

  !
  !----------------------------------------------------------------
  module gspace
  !----------------------------------------------------------------
  !
  use parameters, only : DP, nat
  !
  integer :: ngm, ngl, ngm0, ngms
  ! number of G-vectors corresponding to ecut
  ! number of G-shells
  ! number of G-vectors corresponding to ecut0
  ! number of G-vectors corresponding to ecuts
  integer :: nr1, nr2, nr3, nr
  integer :: nr1s, nr2s, nr3s, nrs
  ! nr?s is for the Sigma cutoff
  ! in espresso this is reserved for the box cutoff,
  ! so we will need to change name at some point
  !
  integer, allocatable :: igtongl(:), nl(:)
  real(DP), allocatable :: g(:,:), gl(:)
  real(DP) :: at (3,3), bg (3,3)
  integer, allocatable :: nls(:)
  ! fft correspondence for the Sigma cutoff
  !
  real(DP) :: tau (3,nat)
  !
  integer, allocatable :: gmap(:,:)
  ! map of G-vectors for folding k+q into first BZ
  real(DP) :: g0vec(3,27)
  ! the folding vectors

  !
  end module gspace
  !
  !----------------------------------------------------------------
  module constants
  !----------------------------------------------------------------
  !
  use parameters, only : DP, alat
  !
  real(DP), parameter :: pi = 3.14159265358979d0
  real(DP), parameter :: twopi = 2.d0 * pi
  real(DP), parameter :: two = 2.d0 
  real(DP), parameter :: tpiba  = 2.d0 * pi / alat
  real(DP), parameter :: tpiba2 = tpiba**2
  real(DP), parameter :: four = 4.d0
  real(DP), parameter :: fpi = four * pi
  real(DP), parameter :: zero = 0.d0
  real(DP), parameter :: one = 1.d0
  real(DP), parameter :: ryd2ev = 13.6058
  real(DP), parameter :: e2 = 2.d0 ! the square of the electron charge 
  complex(DP), parameter :: czero = (0.d0, 0.d0)
  complex(DP), parameter :: cone = (1.d0, 0.d0)
  complex(DP), parameter :: ci = (0.d0, 1.d0)
  !
  end module constants
  !
  !----------------------------------------------------------------
  module kspace
  !----------------------------------------------------------------
  !
  use parameters, only : DP
  !
  real(kind=DP), allocatable :: xk (:,:), wk(:)
  real(kind=DP), allocatable :: xq (:,:), wq(:), eval_occ(:,:)
  ! k-point grid for the calculation of the screened Coulomb interaction
  ! corresponding weights
  !
  end module kspace
  !


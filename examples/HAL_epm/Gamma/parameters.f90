  !
  !----------------------------------------------------------------
  module parameters
  !----------------------------------------------------------------
  !
  ! TYPE DEFINITIONS
  !
  integer, parameter :: dbl = selected_real_kind(14,200)
  integer, parameter :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: i4b = selected_int_kind(9)
  !
  ! STRUCTURE AND POTENTIAL FOR SILICON
  !
!  logical :: tepm = .false.
  logical :: tepm = .true.
  ! switch between empirical potential and local pseudopotential calc
  integer, parameter :: nat = 2, nbnd_occ = 4, nbnd = 8, nbnd_sig = 8
  real(dbl), parameter :: alat = 10.26
  ! lattice parameter of silicon (5.43 A)
  ! without the factor 4 we have the volume of 8 atoms...
  real(dbl), parameter :: omega = 1080.49/4.d0 
  ! unit cell volume of silicon, bohr^3 
  real(dbl), parameter :: v3 = -0.21, v8 = 0.04, v11 = 0.08
  ! diamond structure - only the symmetric form factor is nonzero, Ry
  ! from: Cohen & Bergstresser, PRB 141, 789 (1966)
  !
  ! PLANE WAVES CUTOFF
  !
!  real(dbl), parameter :: ecutwfc = 20.0  ! balde/tosatti: 2 Ry
  real(dbl), parameter :: ecutwfc = 5.0  ! balde/tosatti: 2 Ry
  ! energy cutoff for wave functions in k-space ( in Rydbergs )
  ! 5 Ry and 40 Ry give essentially the same eps^-1(0,0,q)
  ! If I use 2 Ry it bombs out since we shift the G-sphere
  ! for the refolding and most points are lost =O
!  real(dbl), parameter :: ecut0 = 2.0 ! NOTE: ecut0 must be SMALLER than ecut 
  real(dbl), parameter :: ecut0 = 1.1 ! NOTE: ecut0 must be SMALLER than ecut 
  ! energy cutoff for trial diagonalization of Hamiltonian (input to CG - Rydbergs )
!  real(dbl), parameter :: ecuts = 2.5
  real(dbl), parameter :: ecuts = 1.1
  ! energy cutoff for G, W, and Sigma, Ry
  !
  ! BRILLOUIN ZONE SAMPLING
  !
  integer, parameter :: nk0 = 1
  integer, parameter :: nq1 = 6, nq2 = 6, nq3 = 6
  integer, parameter :: q1 = 1, q2 = 1, q3 = 1    
  integer, parameter :: nq = nq1 * nq2 * nq3
  integer, parameter :: nksq = nq, nks = 2 * nksq  
  !
  ! ENERGY SAMPLING
  !
!  real(DP), parameter :: wsigmamin = -100.d0, wsigmamax = 100.d0, deltaw = 0.5, wcoulmax = 100.d0
  real(DP), parameter :: wsigmamin = -100.d0, wsigmamax = 100.d0, deltaw = 0.5, wcoulmax = 100.d0
  ! frequency range for the self-energy (wsigmamin<0, sigmamax>0) - eV
  ! and for the Coulomb  (wcoulmax>0)
  real(DP), parameter :: wgreen_safe = 50.d0
  ! the largest frequency of the Green's function must be at least this value
  ! (otherwise we miss the empty states and the Sigma comes out wrong)
  !
  ! CONVERGENCE PARAMETERS AND ENERGY OFFSET 
  !
  real(dbl), parameter :: eps = 1.d-10
  ! threshold for conjugate gradient convergence
  integer, parameter :: maxter = 200
  ! max number of iterations in cong grad diagonalization
  integer, parameter :: nmax_iter = 30
  ! max n. of iterations in solve_linter
  integer, parameter :: nmix_ph = 6
  ! number of previous iterations used in modified Broyden mixing
  real(DP), parameter :: tr2_ph = 1.d-10
  ! convergence threshold for dvscf
! real(DP), parameter :: tr_cgsolve = 1.0d-8 
  real(DP), parameter :: tr_cgsolve = 1.0d-10
  ! threshold for the iterative solution of the linear system
  ! with thresh>1d-5 the potential does not converge
! real(DP), parameter :: eta = 0.03
  real(DP), parameter :: eta = 0.04
  ! smearing for Green's function, and GW exp^ideltaw sum - Ry
 real(DP), parameter :: alpha_pv = 25.234/13.606 ! this is 2(emax-emin), epm
!  real(DP), parameter :: alpha_pv = 24.804/13.606 ! this is 2(emax-emin), 20 Ry local pp
  ! parameter for the projection over the valence manifold
 real(DP), parameter :: eshift = 10.462/13.6058 ! epm
!  real(DP), parameter :: eshift = 6.12048/13.6058 ! 20 Ry local pp
  ! Energy shift to align the top of the valence band to zero
  !
  ! I/O PARAMETERS
  !
  integer, parameter :: DIRECT_IO_FACTOR = 8
  integer, parameter :: iunwfc0  = 77  ! q-grid
  integer, parameter :: iunwfc   = 78  ! k and k+q grids
  integer, parameter :: iubar    = 80
  integer, parameter :: iudwf    = 79
  integer, parameter :: iudwfp   = 81
  integer, parameter :: iudwfm   = 82
  integer, parameter :: iuncoul  = 86  ! screened Coulomb interaction
  integer, parameter :: iungreen = 87  ! Green's function
  integer, parameter :: iunsigma = 88  ! Self-energy
  integer, parameter :: stdout   = 6
  !
  end module parameters
  !

  !
  !----------------------------------------------------------------
  function allowed (nr)
  !----------------------------------------------------------------
  !
  ! FG - adapted from Modules/fft_scalar.f90 - removed the IBM thing
  !
  ! find if the fft dimension is a good one
  ! a "bad one" is implemented but with awful performances 
  !
  implicit none
  integer :: nr
  logical :: allowed
  integer :: pwr (5), mr, i, fac, p, maxpwr
  integer :: factors( 5 ) = (/ 2, 3, 5, 7, 11 /)
  !
  ! find the factors of the fft dimension
  !
  mr  = nr
  pwr = 0
  factors_loop: do i = 1, 5
     fac = factors (i)
     maxpwr = NINT ( LOG( REAL (mr) ) / LOG( REAL (fac) ) ) + 1
     do p = 1, maxpwr
        if ( mr == 1 ) EXIT factors_loop
        if ( MOD (mr, fac) == 0 ) then
           mr = mr / fac
           pwr (i) = pwr (i) + 1
        endif
     enddo
  end do factors_loop
  !
  IF ( nr /= ( mr * 2**pwr (1) * 3**pwr (2) * 5**pwr (3) * 7**pwr (4) * 11**pwr (5) ) ) &
     CALL error (' allowed ', ' what ?!? ', 1 )
  if ( mr /= 1 ) then
     ! fft dimension contains factors > 11 : no good in any case
     allowed = .false.
  else
     ! fftw and all other cases: no factors 7 and 11
     allowed = ( ( pwr(4) == 0 ) .and. ( pwr(5) == 0 ) )
  endif
  !
  end function allowed
  !

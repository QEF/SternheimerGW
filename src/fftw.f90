  !
  !-----------------------------------------------------------------------
  subroutine cfft3(f,nr1,nr2,nr3,sign)
  !-----------------------------------------------------------------------
  !
  !  driver routine for 3d fft using fftw
  !  FG - from http://www.fftw.org/fftw2_doc/fftw_5.html
  !  and from Paolo Giannozzi's fftw wrapper in the old CPV
  !  
  !  The plans 1 and 2 are for the wavefunction grid
  !  Plans 3 and 4 are for the selfenergy grid
  !
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  !
  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  !   
  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
  !   
  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
  !   
  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
  !   
  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)
  !
  integer :: nr1, nr2, nr3, sign
  complex(DP) f(nr1*nr2*nr3)
  real(DP) :: fac
!@  integer :: ip, plan(4)
  integer :: ip
  integer*8  plan(4) ! see: http://www.fftw.org/fftw2_doc/fftw_5.html#SEC63
                     ! without this it does not work on Opteron with PGI
  save plan
  data plan /0,0,0,0/
  !
  !  plans 1,2 for the wavefunction grid
  !
  if       (sign.eq.1) then
     ip=1
  else if (sign.eq.-1) then
     ip=2
  else 
     call error('cfft3','sign unexpected',sign)
  end if
  !
  if (plan(ip).eq.0) call fftw3d_f77_create_plan  &
     &   (plan(ip),nr1,nr2,nr3,sign,FFTW_ESTIMATE+FFTW_IN_PLACE)
  !
  call fftwnd_f77_one(plan(ip), f, 0)
  !
  if (sign.eq.-1) then
     fac = 1.d0/float(nr1*nr2*nr3)
     call DSCAL(2*nr1*nr2*nr3, fac, f, 1)
  end if
  !
  return
  end
  !
  !-----------------------------------------------------------------------
  subroutine cfft3s(f,nr1,nr2,nr3,sign)
  !-----------------------------------------------------------------------
  !
  !  driver routine for 3d fft using fftw
  !  FG - from http://www.fftw.org/fftw2_doc/fftw_5.html
  !  and from Paolo Giannozzi's fftw wrapper in the old CPV
  !  
  !  The plans 1 and 2 are for the wavefunction grid
  !  Plans 3 and 4 are for the selfenergy grid
  !
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  !
  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  !   
  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
  !   
  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
  !   
  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
  !   
  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)
  !
  integer :: nr1, nr2, nr3, sign
  complex(DP) f(nr1*nr2*nr3)
  real(DP) :: fac
!@  integer :: ip, plan(4)
  integer :: ip
  integer*8  plan(4) ! see: http://www.fftw.org/fftw2_doc/fftw_5.html#SEC63
                     ! without this it does not work on Opteron with PGI
  save plan
  data plan /0,0,0,0/
  !
  ! plans 3,4 for the selfenergy grid
  !
  if       (sign.eq.1) then
     ip=3
  else if (sign.eq.-1) then
     ip=4
  else 
     call error('cfft3','sign unexpected',sign)
  end if
  !
  if (plan(ip).eq.0) call fftw3d_f77_create_plan  &
     &   (plan(ip),nr1,nr2,nr3,sign,FFTW_ESTIMATE+FFTW_IN_PLACE)
  !
  call fftwnd_f77_one(plan(ip), f, 0)
  !
  if (sign.eq.-1) then
     fac = 1.d0/float(nr1*nr2*nr3)
     call DSCAL(2*nr1*nr2*nr3, fac, f, 1)
  end if
  !
  return
  end
  !

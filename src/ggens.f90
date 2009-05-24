  !
  !----------------------------------------------------------------
  subroutine ggens ( gcutms )
  !----------------------------------------------------------------
  !
  ! determine G-vectors within the cutoff from the
  ! array already created in ggen
  ! FG - 
  !
  use parameters
  use gspace
  implicit none
  !
  real(kind=DP) :: gcutms
  !
  ! local variables
  !
  integer :: n1, n2, n3, i, j, k, ipol, ng, igl
  !
  do ng = 1, ngm
    if ( gl( igtongl (ng) ) .le. gcutms ) ngms = ng
  enddo
  !
  allocate ( nls(ngms) )
  !
  !     Now set nl with the correct fft correspondence
  !
  do ng = 1, ngms
     ! n1 is going to be i+1, folded to positive when <= 0
     n1 = nint (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, ng) * at (3, 1) ) + 1
     if (n1.lt.1) n1 = n1 + nr1s
     ! n2 is going to be j+1, folded to positive when <= 0
     n2 = nint (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, ng) * at (3, 2) ) + 1
     if (n2.lt.1) n2 = n2 + nr2s
     ! n3 is going to be k+1, folded to positive when <= 0
     n3 = nint (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, ng) * at (3, 3) ) + 1
     if (n3.lt.1) n3 = n3 + nr3s
     !
     if (n1.le.nr1s.and.n2.le.nr2s.and.n3.le.nr3s) then
       nls (ng) = n1 + (n2 - 1) * nr1s + (n3 - 1) * nr1s * nr2s
     else
        call error('ggens','Mesh too small?',ng)
     endif
  enddo
  !
  write(6,'(4x,"ngms  = ",i5)') ngms
  !
  ! total number of real-space grid points
  !
  nrs = nr1s * nr2s * nr3s
  write(6,'(4x,"nr1s  = ",i10)') nr1s
  write(6,'(4x,"nr2s  = ",i10)') nr2s
  write(6,'(4x,"nr3s  = ",i10)') nr3s
  write(6,'(4x,"nrs   = ",i10)') nrs
  !
  end subroutine ggens
  !

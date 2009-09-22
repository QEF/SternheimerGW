  !
  !-----------------------------------------------------------------------
  subroutine kinetic_psinc ( g2kin, t_psinc )
  !-----------------------------------------------------------------------
  !
  ! ******** THIS IS JUST CRAZY... NOT TESTED AND IN ANY CASE USELESS *****
  ! ******** THE KIN EN IN PSINC IS NOT STRICTLY LOCALIZED, SO THIS IS REALLY A
  ! WASTE OF (COMPUTER) TIME
  !
  ! construct k-dependent kinetic energy in psinc basis starting from
  ! the kinetic energy in G-space
  !
  !-----------------------------------------------------------------------
  !
  use gspace
  use constants
  use parameters
  implicit none
  !
  real(dbl) :: g2kin(ngm)
  complex(kind=DP) :: t_psinc(nr,nr)
  !
  integer :: ir1, ir2, ig, i, j, k
  complex(kind=DP) :: d(nr), td (ngm), prod
  real(dbl) :: len, r1(3), r2(3)
  !
  do ir2 = 1, nr
   !
   ! define D_2 in real space
   !
   d = czero 
   d ( ir2 ) = cone
   !
   ! FFT to G-space
   !
   call cfft3 ( d, nr1, nr2, nr3, -1)
   !
   ! apply kinetic energy in G-space: td = T * D_2
   !
   do ig = 1, ngm
     td (ig) = g2kin (ig) * d ( nl(ig) )
   enddo
   !
   ! find real-space grid point
   !
   call findr (ir2, i, j, k)
   r2 = float(i)/float(nr1) * at(:,1) &
      + float(j)/float(nr2) * at(:,2) &
      + float(k)/float(nr3) * at(:,3)
   !
   do ir1 = ir2, nr
     !
     ! define D_1 in real space 
     !
     d = czero
     d ( ir1 ) = cone
     !
     ! FFT to G-space
     !
     call cfft3 ( d, nr1, nr2, nr3, -1)
     !
     ! scalar product in G-space: < D_1 | T | D_2 >
     !
     prod = czero
     do ig = 1, ngm
       prod = prod + conjg ( d ( nl(ig) ) ) * td (ig) 
     enddo
     t_psinc ( ir1, ir2) = prod
     t_psinc ( ir2, ir1) = prod
     !
     call findr (ir1, i, j, k)
     r1 = float(i)/float(nr1) * at(:,1) &
        + float(j)/float(nr2) * at(:,2) &
        + float(k)/float(nr3) * at(:,3)
     !
     call dist (r1, r2, len)
     !
     write (6,'(2f15.8)') len * alat, abs(t_psinc(ir1,ir2))
     !
   enddo
  enddo
  !
  return
  end subroutine kinetic_psinc
  !
  !-----------------------------------------------------------------------
  subroutine findr (ir, i, j, k)
  !-----------------------------------------------------------------------
  !
  use gspace
  implicit none
  integer :: ir, i, j, k
  !
  i = mod ( ir-1, nr1) + 1
  j = mod ( (ir - (i-1))/nr1, nr2) + 1
  k = (ir - (i-1) - (j-1)*nr2 )/(nr1 * nr2) + 1
  !
  return
  end subroutine findr
  !
  !-----------------------------------------------------------------------
  subroutine dist (a, b, d)
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  real(DP) :: a(3), b(3), d
  !
  d = (a(1)-b(1))**2.d0 + (a(2)-b(2))**2.d0 + (a(3)-b(3))**2.d0
  d = sqrt(d)
  !
  return
  end subroutine dist
  !

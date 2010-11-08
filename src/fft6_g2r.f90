  !
  !----------------------------------------------------------------
  subroutine fft6_g2r ( f_g, f_r)
  !----------------------------------------------------------------
  ! 
  !----------------------------------------------------------------
  !
  use parameters, only : dbl, omega
  use constants, only : czero
  use gspace, only : ngms, nrs, nr1s, nr2s, nr3s, nls
  !
  implicit none
  integer :: ig, igp, ir, irp
  complex(dbl) :: f_r (nrs,nrs), f_g(ngms,ngms), aux(nrs)
  !
  do ig = 1, ngms
    aux = czero
    do igp = 1, ngms
      aux(nls(igp)) = f_g(ig,igp)
    enddo
    call cfft3s ( aux, nr1s, nr2s, nr3s,  1)
    f_r(ig,1:nrs) = aux / omega
  enddo
  !
  ! the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
  ! following teh convention set in the paper
  ! [because the standard transform is sum_G f(G) exp(iGr) ]
  !
  do irp = 1, nrs
    aux = czero
    do ig = 1, ngms
      aux(nls(ig)) = conjg ( f_r(ig,irp) )
    enddo
    call cfft3s ( aux, nr1s, nr2s, nr3s,  1)
    f_r(1:nrs,irp) = conjg ( aux )
  enddo
  !
  end subroutine fft6_g2r
  !

SUBROUTINE fft6_g2r (f_g, f_r)
! 6-D fourier transform used for W(G,G') and G(G,G') 
! which agrees with convention in FG's paper. Although in matrix elements
! apparently there is some problem with the exchange part of the energy requiring
! the hack using -k instead of k when looking at quasi-particle corrections.

  USE io_files,         ONLY : prefix, iunigk
  USE gwsigma,          ONLY : ngmsig, nrsig, nr1sig, nr2sig, nr3sig, nlsig
  USE qpoint,           ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE kinds,            ONLY : DP 
  USE cell_base,        ONLY : omega, alat
  USE control_gw,       ONLY : lgamma
  USE wvfct,            ONLY : npw,npwx, igk

  IMPLICIT NONE

  INTEGER                          :: ios
  INTEGER                          :: ig, igp, ir, irp
  COMPLEX(DP)                      :: f_r (nrsig,nrsig)
  COMPLEX(DP)                      :: f_g(ngmsig,ngmsig)
  COMPLEX(DP)                      :: aux(nrsig)
  COMPLEX(DP)                      :: czero

 !HL CONSTANTS should be defined GLOBALLY! Since they are CONSTANT and they DON'T CHANGE!  
 !ALLOCATE( f_g(ngmsig, ngmsig) )
 !ALLOCATE( f_r(nrsig, nrsig) )

  czero = (0.0d0, 0.0d0)
  f_r(:,:) = czero

  ! The greens_function and coulomb should both be stored with the correct folding from the igkq
  ! to the gamma centered mesh. 

  do ig = 1, ngmsig
    aux(:) = czero
    do igp = 1, ngmsig
      aux(nlsig(igp)) = f_g(ig,igp)
    enddo

  ! cft3s
  ! +-1 parallel 3d fft for rho and for the potential. +-2 parallel 3d fft for wave fxns. 
  ! +G -> R, and minus (-) R -> G \int_R f(R)exp(-iG*R)/Omega
  ! call cft3 (evc_r, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
  ! I think

    call cft3 (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, +1)

     f_r(ig,1:nrsig) = aux / omega

  !  do irp = 1, nrsig
  !     f_r(ig, irp) = aux(irp) / omega
  !     write(6,*)f_r(ig,irp) 
  !  enddo
  enddo

  ! the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
  ! following teh convention set in the paper
  ! [because the standard transform is sum_G f(G) exp(iGr) ]

  do irp = 1, nrsig
    aux = czero
    do ig = 1, ngmsig
      aux(nlsig(ig)) = conjg( f_r(ig,irp) )
    enddo
  ! HL call cft3s ( aux, nr1sig, nr2sig, nr3sig,  1)
    call cft3 (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, +1)
    f_r(1:nrsig,irp) = conjg ( aux )
  enddo

END SUBROUTINE fft6_g2r

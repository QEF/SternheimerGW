SUBROUTINE fft6_g2r (f_g, f_r)
! 6-D fourier transform used for W(G,G') and G(G,G') which agress with convention in FG's paper.
  USE gwsigma,          ONLY : ngmsig, nrsig, nr1sig, nr2sig, nr3sig, nlsig
  USE kinds,            ONLY : DP 
  USE cell_base,        ONLY : omega, alat

  IMPLICIT NONE

  INTEGER                          :: ig, igp, ir, irp
  COMPLEX(DP)                      :: f_r (nrsig,nrsig)
  COMPLEX(DP)                      :: f_g(ngmsig,ngmsig)
  COMPLEX(DP)                      :: aux(nrsig)
  COMPLEX(DP)                      :: czero

 !HL CONSTANTS should be defined GLOBALLY! Since they are CONSTANT and they DON'T CHANGE!  
 !ALLOCATE( f_g(ngmsig, ngmsig) )
 !ALLOCATE( f_r(nrsig, nrsig) )

  czero = (0.0d0, 0.0d0)

!Initializing f_r causes things to hang... I think since the size of this array is only determined 
!Upon execution of ggensig hardwiring the array size to nrsig is causing lots and lots of problems.
! f_r(:,:) = (0.0d0, 0.0d0) 


  f_r(:,:) = czero

  do ig = 1, ngmsig
    !write(6,*)ig
    aux(:) = czero

    do igp = 1, ngmsig
      aux(nlsig(igp)) = f_g(ig,igp)
    enddo

!   cft3s
!   +-1 parallel 3d fft for rho and for the potential. +-2 parallel 3d fft for wave fxns. 
!   +G -> R, and minus (-) R -> G \int_R f(R)exp(-iG*R)/Omega
!   call cft3s (evc_r, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)

!   I think

    call cft3s (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, +1)

!   f_r(ig,1:nrsig) = aux / omega
    do irp = 1, nrsig
       f_r(ig, irp) = aux(irp) / omega
      !write(6,*)f_r(ig,irp) 
    enddo
  enddo

  !
  ! the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
  ! following teh convention set in the paper
  ! [because the standard transform is sum_G f(G) exp(iGr) ]
  !

  do irp = 1, nrsig
    aux = czero

    do ig = 1, ngmsig
      aux(nlsig(ig)) = conjg ( f_r(ig,irp) )
    enddo

! HL call cft3s ( aux, nr1sig, nr2sig, nr3sig,  1)
! Do I need to write a new plan for this FFT?  is there any difference between dimension in nr1sig, and 
! "nr1sigx" I don't think so but why does quantum espresso insist on repeated dimension indices in all 
! these FFT calls. 

    call cft3s (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, +1)
    f_r(1:nrsig,irp) = conjg ( aux )
    
  enddo
  !
END SUBROUTINE fft6_g2r

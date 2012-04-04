SUBROUTINE sigma_r2g(sigma, sigma_g)

! 6-D fourier transform used for taking Sigma back in to G space after the product has been formed in 
! the cell periodic fxns in real space. 

  USE io_files,          ONLY : prefix, iunigk
  USE gwsigma,           ONLY : ngmsig, nrsig, nr1sig, nr2sig, nr3sig, nlsig
  USE kinds,             ONLY : DP
  USE control_gw,        ONLY : lgamma
  USE cell_base,         ONLY : omega, alat
  USE wvfct,             ONLY : npw, npwx, igk
  USE freq_gw,           ONLY : nwsigma
  USE qpoint,            ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gsmooth,           ONLY : nrxxs, nrx1s, nrx2s, nrx3s, nr1s, nr2s, nr3s, nls
  USE gvect,             ONLY : nrxx, g, nl, nr1, nr2, nr3, nrx1, nrx2, nrx3
  USE fft_scalar,        ONLY : cfft3d

  IMPLICIT NONE

  INTEGER     :: ios 
  INTEGER     :: ig, igp, ir, irp, iw

!  Self-Energy Grid
! COMPLEX(DP) :: sigma (nrsig,nrsig,nwsigma)
! COMPLEX(DP) :: sigma_g(ngmsig,ngmsig,nwsigma) 
  COMPLEX(DP) :: sigma (nrsig,nrsig,1)
  COMPLEX(DP) :: sigma_g(ngmsig,ngmsig,1) 
  COMPLEX(DP) :: aux(nrsig)

  COMPLEX(DP) :: czero

! HL
  sigma_g = (0.0d0,0.0d0)
!Only 1 frequency omega =0 when doing sigma_ex.
!  do iw = 1, nwsigma
      iw = 1
      do ir = 1, nrsig
        aux = (0.0d0, 0.0d0)
        do irp = 1, nrsig
          aux(irp) = sigma(ir,irp,iw)
        enddo

      !Sigma Grid.
      !call cft3s (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)
       call cfft3d (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)

      !Smooth Grid
      !call cft3s (aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)

        do igp = 1, ngmsig
           !Sigma Grid.
           sigma (ir,igp,iw) = aux( nlsig(igp))
           !Smooth Grid. 
           !sigma (ir,igp,iw) = aux(nls(igp))
        enddo
      enddo

      do igp = 1, ngmsig
        aux = czero
        do ir = 1, nrsig
          aux(ir) = conjg ( sigma(ir,igp,iw) )
        enddo

      !Sigma Grid.
       !call cft3s (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)
        call cfft3d (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)

        do ig = 1, ngmsig
           sigma (ig,igp,iw) = conjg ( aux( nlsig( ig )) ) * omega
        enddo
      enddo
!  enddo 

! HL Static Approximation single frequency.
! everything beyond ngmsig is garbage

   do ig = ngmsig + 1, nrsig
    do igp = ngmsig + 1, nrsig
       !do iw = 1, nwsigma
        sigma (ig,igp,1) = (0.0d0, 0.0d0)
       !enddo
    enddo
   enddo
   
    sigma_g = sigma(1:ngmsig,1:ngmsig,:)
  
END SUBROUTINE sigma_r2g


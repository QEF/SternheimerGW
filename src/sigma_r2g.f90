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

  IMPLICIT NONE

  INTEGER     :: ios 
  INTEGER     :: ig, igp, ir, irp, iw
!  COMPLEX(DP) :: sigma (nrsig,nrsig,nwsigma)
!  COMPLEX(DP) :: sigma_g(ngmsig,ngmsig,nwsigma) 

  COMPLEX(DP) :: sigma (nrxxs,nrxxs,1)
  COMPLEX(DP) :: sigma_g(ngmsig,ngmsig,1) 
  COMPLEX(DP) :: aux(nrxxs)

  COMPLEX(DP) :: czero

! HL
  sigma_g = (0.0d0,0.0d0)

  !do iw = 1, nwsigma
   iw = 1
      do ir = 1, nrxxs
        aux = (0.0d0, 0.0d0)
        do irp = 1, nrxxs
          !aux(irp) = sigma(ir,irp,iw,1)
          aux(irp) = sigma(ir,irp,iw)
        enddo

        !call cft3 (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)
        call cft3s (aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)

        do igp = 1, ngmsig
          !sigma (ir,igp,1,1) = aux( nlsig(igp))
          sigma (ir,igp,iw) = aux(nls(igp))
        enddo
      enddo

      do igp = 1, ngmsig
        aux = czero
        do ir = 1, nrxxs
          aux(ir) = conjg ( sigma(ir,igp,iw) )
        enddo

       !call cft3 (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)

        call cft3s (aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)

        do ig = 1, ngmsig
           sigma (ig,igp,iw) = conjg ( aux( nls( ig )) ) * omega
        enddo
      enddo
 ! enddo 

 !HL Static Approximation single frequency.
 !everything beyond ngmsig is garbage

   do ig = ngmsig + 1, nrxxs
    do igp = ngmsig + 1, nrxxs
 !      do iw = 1, nwsigma
        !sigma (ig,igp,iw) = (0.0d0, 0.0d0)
        sigma (ig,igp,1) = (0.0d0, 0.0d0)
 !      enddo
    enddo
   enddo

! #ifdef __PARA
!     if (me.eq.1.and.mypool.eq.1) then
! #endif
   
       sigma_g = sigma(1:ngmsig,1:ngmsig,:)

! #ifdef __PARA
!    endif
! #endif
  
END SUBROUTINE sigma_r2g


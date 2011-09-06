SUBROUTINE sigma_r2g(sigma, sigma_g)

! 6-D fourier transform used for taking Sigma back in to G space after the product has been formed in 
! the cell periodic fxns in real space. 

  USE io_files,         ONLY : prefix, iunigk
  USE gwsigma,           ONLY : ngmsig, nrsig, nr1sig, nr2sig, nr3sig, nlsig
  USE qpoint,           ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE kinds,            ONLY : DP
  USE control_gw,       ONLY : lgamma
  USE cell_base,        ONLY : omega, alat
  USE wvfct,            ONLY : npw,npwx, igk
  USE freq_gw,          ONLY : nwsigma

  IMPLICIT NONE

  INTEGER     :: ios 
  INTEGER     :: ig, igp, ir, irp, iw
  COMPLEX(DP) :: sigma (nrsig,nrsig,nwsigma,1)
  COMPLEX(DP) :: sigma_g(ngmsig,ngmsig,nwsigma) 
  COMPLEX(DP) :: aux(nrsig)
  COMPLEX(DP) :: czero



!Again reading in mapping indices for FFTs. 
  
    do iw = 1, nwsigma
      do ir = 1, nrsig
        aux = (0.0d0, 0.0d0)
        do irp = 1, nrsig
          aux(irp) = sigma(ir,irp,iw,1)
        enddo

!       call cfft3s ( aux, nr1s, nr2s, nr3s, -1)
!       call cft3s (aux, nr1sig, nr2sig, nr3s, ngm1sig, ngm2sig, ngmsig3, -2)

        call cft3s (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)

        do igp = 1, ngmsig
          sigma (ir,igp,iw,1) = aux( nlsig(igp ) )
        enddo
      enddo

      do igp = 1, ngmsig
        aux = czero

        do ir = 1, nrsig
          aux(ir) = conjg ( sigma(ir,igp,iw,1) )
        enddo

!       call cfft3s ( aux, nr1s, nr2s, nr3s, -1)
  
        call cft3s (aux, nr1sig, nr2sig, nr3sig, nr1sig, nr2sig, nr3sig, -1)

        do ig = 1, ngmsig
          sigma (ig,igp,iw,1) = conjg ( aux( nlsig( ig )) ) * omega
        enddo
      enddo
    enddo

    !
    ! everything beyond ngmsig is garbage
    !

    do ig = ngmsig + 1, nrsig
     do igp = ngmsig + 1, nrsig
      do iw = 1, nwsigma
         sigma (ig,igp,iw,1) = (0.0d0, 0.0d0)
      enddo
     enddo
    enddo

! #ifdef __PARA
!     if (me.eq.1.and.mypool.eq.1) then
! #endif
 
   sigma_g = sigma(1:ngmsig,1:ngmsig,:,1)

! #ifdef __PARA
!    endif
! #endif
  
END SUBROUTINE sigma_r2g


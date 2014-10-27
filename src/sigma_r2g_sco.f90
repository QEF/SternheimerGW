SUBROUTINE sigma_r2g_sco(sigma, sigma_g)
! 6-D fourier transform used for taking sigma_ex back in to G space after 
! the product has been formed in real space.

  USE io_files,          ONLY : prefix, iunigk
  USE kinds,             ONLY : DP
  USE control_gw,        ONLY : lgamma
  USE cell_base,         ONLY : omega, alat
  USE wvfct,             ONLY : npw, npwx, igk
  USE freq_gw,           ONLY : nwsigma
  USE qpoint,            ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,           ONLY : ngmsco, nrsco, nr1sco, nr2sco, nr3sco, nlsco
  USE fft_scalar,    ONLY : cfft3d

  IMPLICIT NONE

  INTEGER     :: ios 
  INTEGER     :: ig, igp, ir, irp, iw
  COMPLEX(DP) :: sigma (nrsco,nrsco, nwsigma)
  COMPLEX(DP) :: sigma_g (ngmsco,ngmsco, nwsigma) 
  COMPLEX(DP) :: aux(nrsco)
  COMPLEX(DP) :: czero

! HL
    sigma_g = (0.0d0,0.0d0)
    do iw = 1, nwsigma
      do ir = 1, nrsco
        aux = (0.0d0, 0.0d0)
        do irp = 1, nrsco
          aux(irp) = sigma(ir,irp,iw)
        enddo
        call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
        do igp = 1, ngmsco
           sigma (ir,igp, iw) = aux(nlsco(igp))
        enddo
      enddo
      do igp = 1, ngmsco
        aux = czero
        do ir = 1, nrsco
          aux(ir) = conjg ( sigma(ir,igp,iw) )
        enddo
        call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
        do ig = 1, ngmsco
           sigma (ig,igp,iw) = conjg ( aux( nlsco( ig )) ) * omega
        enddo
      enddo
    enddo

    do ig = ngmsco + 1, nrsco
       do igp = ngmsco + 1, nrsco
          do iw = 1, nwsigma
             sigma (ig, igp, iw) = (0.0d0, 0.0d0)
          enddo 
       enddo
    enddo
    sigma_g = sigma(1:ngmsco,1:ngmsco,:)
  
END SUBROUTINE sigma_r2g_sco


SUBROUTINE sigma_r2g_ex(sigma, sigma_g)
! 6-D fourier transform used for taking sigma_ex back in to G space after 
! the product has been formed in real space.

  USE io_files,          ONLY : prefix, iunigk
  USE kinds,             ONLY : DP
  USE control_gw,        ONLY : lgamma
  USE cell_base,         ONLY : omega, alat
  USE wvfct,             ONLY : npw, npwx, igk
  USE freq_gw,           ONLY : nwsigma
  USE qpoint,            ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,           ONLY : ngmsex, nrsex, nr1sex, nr2sex, nr3sex, nlsex
  USE fft_scalar,        ONLY : cfft3d

  IMPLICIT NONE

  INTEGER     :: ios 
  INTEGER     :: ig, igp, ir, irp, iw
  COMPLEX(DP) :: sigma (nrsex,nrsex)
  COMPLEX(DP) :: sigma_g(ngmsex,ngmsex) 
  COMPLEX(DP) :: aux(nrsex)
  COMPLEX(DP) :: czero

! HL
      sigma_g = (0.0d0,0.0d0)
      czero   = (0.0d0, 0.0d0)
      do ir = 1, nrsex
        aux = (0.0d0, 0.0d0)
        do irp = 1, nrsex
          aux(irp) = sigma(ir,irp)
        enddo

      ! call cft3s (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, -1)
        call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, -1)

        do igp = 1, ngmsex
           sigma (ir,igp) = aux( nlsex(igp))
        enddo
      enddo

      do igp = 1, ngmsex
        aux = czero
        do ir = 1, nrsex
          aux(ir) = conjg ( sigma(ir,igp) )
        enddo

       !call cft3s (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, -1)
        call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, -1)

        do ig = 1, ngmsex
           sigma (ig,igp) = conjg ( aux( nlsex( ig )) ) * omega
        enddo
      enddo

     do ig = ngmsex + 1, nrsex
        do igp = ngmsex + 1, nrsex
           sigma (ig,igp) = (0.0d0, 0.0d0)
        enddo
     enddo
   
     sigma_g = sigma(1:ngmsex,1:ngmsex)
  
END SUBROUTINE sigma_r2g_ex


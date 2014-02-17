subroutine fft6(f_g, f_r, conv)
  USE kinds,         ONLY : DP
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE fft_scalar,    ONLY : cfft3d
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at

IMPLICIT NONE

INTEGER :: conv
COMPLEX(DP)  :: f_g(ngmsco, ngmsco)
COMPLEX(DP)  :: f_r(nrsco, nrsco)
COMPLEX(DP)  :: aux (nrsco)
COMPLEX(DP)  :: ci, czero
INTEGER :: ig, igp, irr, icounter, ir, irp

ci = dcmplx(0.0d0, 1.d0)
czero = dcmplx(0.0d0, 0.0d0)

if(conv.eq.1) then
            do ig = 1, ngmgrn
               aux(:) = czero
               do igp = 1, ngmgrn
                  aux(nlsco(igp)) = f_g(ig,igp)
               enddo
               call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
               do irp = 1, nrsco
                  f_r(ig, irp) = aux(irp) / omega
               enddo
            enddo
            do irp = 1, nrsco
               aux = czero
                    do ig = 1, ngmgrn
                           aux(nlsco(ig)) = conjg(f_r(ig,irp))
                    enddo
               call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
               f_r(1:nrsco,irp) = conjg ( aux )
            enddo
else if (conv.eq.-1) then
    do ir = 1, nrsco
      aux = (0.0d0, 0.0d0)
      do irp = 1, nrsco
         aux(irp) = f_r(ir,irp)
      enddo
      call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
      do igp = 1, ngmsco
         f_r (ir, igp) = aux(nlsco(igp))
      enddo
    enddo
    do igp = 1, ngmsco
      aux = czero
      do ir = 1, nrsco
        aux(ir) = conjg (f_r(ir,igp))
      enddo
      call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
      do ig = 1, ngmsco
         f_r(ig, igp) = conjg ( aux( nlsco( ig )) ) * omega
      enddo
    enddo
else 
       write(6,'("Shit Buzz...")')
end if

end subroutine fft6

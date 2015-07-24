subroutine fft6_g(f_g, f_r, fc, gmapsym, eigv, isym, nig0, conv)
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba2, tpiba, omega, alat, at
  USE fft_base,       ONLY : dffts
  USE gvect,         ONLY : g, ngm, nl
  USE fft_interfaces, ONLY : invfft, fwfft
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE fft_custom,     ONLY : fft_cus, set_custom_grid, ggent, gvec_init

IMPLICIT NONE

TYPE(fft_cus) fc 
COMPLEX(DP)  :: f_g(fc%ngmt, fc%ngmt)
COMPLEX(DP)  :: f_r(fc%dfftt%nnr, fc%dfftt%nnr)
COMPLEX(DP)  :: aux (fc%dfftt%nnr)
COMPLEX(DP)  :: ci, czero
INTEGER :: ig, igp, irr, icounter, ir, irp
INTEGER :: conv

INTEGER     :: isym
INTEGER     :: gmapsym  (ngm, nrot) 
COMPLEX(DP) :: eigv     (ngm, nrot)  
COMPLEX(DP) :: phase
INTEGER     :: nig0
COMPLEX(DP) :: pwg0(fc%dfftt%nnr)

ci = dcmplx(0.0d0, 1.d0)
czero = dcmplx(0.0d0, 0.0d0)
if(conv.eq.1) then
            if(nig0.gt.1) then
               pwg0(:) = dcmplx(0.0d0, 0.0d0)
               phase = eigv(nig0, isym)
               pwg0(fc%nlt(nig0)) = dcmplx(1.0d0, 0.0d0)
               CALL invfft('Custom', pwg0(:), fc%dfftt)
            endif
            do ig = 1, fc%ngmt
               aux(:) = czero
               do igp = 1, fc%ngmt
                  phase = eigv(igp,isym)
                  aux(fc%nlt(gmapsym(igp,invs(isym)))) = f_g(ig,igp)!*conjg(phase)
               enddo
               call invfft('Custom', aux, fc%dfftt)
               do irp = 1, fc%dfftt%nnr
                  f_r(ig, irp) = aux(irp) / omega
               enddo
            enddo
            if(nig0.gt.1) then
              do ig=1, fc%ngmt
                do ir = 1, fc%dfftt%nnr  
                   f_r(ig,ir) = f_r(ig,ir)*pwg0(ir) 
                enddo
              enddo
            endif              
            if(nig0.gt.1) then
               pwg0(:) = dcmplx(0.0d0, 0.0d0)
               !pwg0(fc%nlt(gmapsym(nig0,invs(isym)))) = dcmplx(1.0d0, 0.0d0)
               pwg0(fc%nlt(nig0)) = dcmplx(1.0d0, 0.0d0)
               CALL invfft('Custom', pwg0(:), fc%dfftt)
               pwg0(:) = conjg(pwg0(:))
            endif
            do irp = 1, fc%dfftt%nnr 
               aux = czero
               do ig = 1, fc%ngmt
!                 phase = eigv(ig, isym)
                  aux(fc%nlt(gmapsym(ig,invs(isym)))) = conjg(f_r(ig,irp))!*conjg(phase)
               enddo
               call invfft('Custom', aux, fc%dfftt)
               f_r(1:fc%dfftt%nnr,irp) = conjg ( aux )
            enddo
            if(nig0.gt.1) then
              do irp=1, fc%dfftt%nnr
                do ir = 1, fc%dfftt%nnr  
                   f_r(ir, irp) = f_r(ir,irp)*pwg0(ir) 
                enddo
              enddo
            endif              
else if (conv.eq.-1) then
    do ir = 1, fc%dfftt%nnr
      aux = (0.0d0, 0.0d0)
      do irp = 1, fc%dfftt%nnr
         aux(irp) = f_r(ir,irp)
      enddo
      call fwfft('Custom', aux, fc%dfftt)
      do igp = 1, fc%ngmt
         f_r (ir, igp) = aux(fc%nlt(igp))
      enddo
    enddo
    do igp = 1, fc%ngmt
      aux = czero
      do ir = 1, fc%dfftt%nnr
         aux(ir) = conjg (f_r(ir,igp))
      enddo
      call fwfft ('Custom', aux, fc%dfftt)
      do ig = 1, fc%ngmt
         f_r(ig, igp) = conjg (aux(fc%nlt(ig))) * omega
      enddo
    enddo
    f_g(1:fc%ngmt, 1:fc%ngmt) = f_r(1:fc%ngmt,1:fc%ngmt)
else 
    call errore (' FFT routines',' Wrong switch',1)
end if
end subroutine fft6_g

!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
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
COMPLEX(DP) :: f_g(fc%ngmt, fc%ngmt)
COMPLEX(DP) :: f_r(fc%dfftt%nnr, fc%dfftt%nnr)
COMPLEX(DP) :: aux (fc%dfftt%nnr)
COMPLEX(DP) :: ci, czero
COMPLEX(DP) :: eigv     (ngm, nrot)  
COMPLEX(DP) :: phase
COMPLEX(DP) :: pwg0(fc%dfftt%nnr)
INTEGER   :: nig0
INTEGER   :: ig, igp, irr, icounter, ir, irp
INTEGER   :: conv
INTEGER   :: isym
INTEGER   :: gmapsym  (ngm, nrot) 

ci = dcmplx(0.0d0, 1.d0)
czero = dcmplx(0.0d0, 0.0d0)
if(conv.eq.1) then
            if((nig0.gt.1) .and. (nig0.le.fc%ngmt)) then
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
            if((nig0.gt.1) .and. (nig0.le.fc%ngmt)) then
              do ig=1, fc%ngmt
                do ir = 1, fc%dfftt%nnr  
                   f_r(ig,ir) = f_r(ig,ir)*pwg0(ir) 
                enddo
              enddo
            endif              
            if((nig0.gt.1) .and. (nig0.le.fc%ngmt)) then
               pwg0(:) = dcmplx(0.0d0, 0.0d0)
               pwg0(fc%nlt(nig0)) = dcmplx(1.0d0, 0.0d0)
               CALL invfft('Custom', pwg0(:), fc%dfftt)
               pwg0(:) = conjg(pwg0(:))
            endif
            do irp = 1, fc%dfftt%nnr 
               aux = czero
               do ig = 1, fc%ngmt
                  aux(fc%nlt(gmapsym(ig,invs(isym)))) = conjg(f_r(ig,irp))
               enddo
               call invfft('Custom', aux, fc%dfftt)
               f_r(1:fc%dfftt%nnr,irp) = conjg ( aux )
            enddo
            if((nig0.gt.1) .and. (nig0.le.fc%ngmt)) then
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

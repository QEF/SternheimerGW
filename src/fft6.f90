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
!> Provide routines for 6 dimensional FFT.
!!
!! This module allows to do a Fourier transform of quantities that have to
!! spacial coordinates into reciprocal space and the reverse
!! \f{equation}{
!!   f(r, r') \longrightarrow_{\text{fwfft6}}  f(G, G')
!!            \longrightarrow_{\text{invfft6}} f(r, r')~.
!! \f}
MODULE fft6_module

  IMPLICIT NONE

CONTAINS

subroutine fft6(f_g, f_r, fc, conv)
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba2, tpiba, omega, alat, at
  USE fft_base,       ONLY : dffts
  USE fft_interfaces, ONLY : invfft, fwfft
  USE fft_custom,     ONLY : fft_cus, set_custom_grid, ggent, gvec_init

IMPLICIT NONE

TYPE(fft_cus) fc 
COMPLEX(DP)  :: f_g(fc%ngmt, fc%ngmt)
COMPLEX(DP)  :: f_r(fc%dfftt%nnr, fc%dfftt%nnr)
COMPLEX(DP)  :: aux (fc%dfftt%nnr)
COMPLEX(DP)  :: ci, czero
INTEGER :: ig, igp, irr, icounter, ir, irp
INTEGER :: conv

ci = dcmplx(0.0d0, 1.d0)
czero = dcmplx(0.0d0, 0.0d0)

if(conv.eq.1) then
            do ig = 1, fc%ngmt
               aux(:) = czero
               do igp = 1, fc%ngmt
                  aux(fc%nlt(igp)) = f_g(ig,igp)
               enddo
               call invfft('Custom', aux, fc%dfftt)
               do irp = 1, fc%dfftt%nnr
                  f_r(ig, irp) = aux(irp) / omega
               enddo
            enddo
            do irp = 1, fc%dfftt%nnr 
               aux = czero
                    do ig = 1, fc%ngmt
                           aux(fc%nlt(ig)) = conjg(f_r(ig,irp))
                    enddo
               call invfft('Custom', aux, fc%dfftt)
               f_r(1:fc%dfftt%nnr,irp) = conjg ( aux )
            enddo
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
         f_r(ig, igp) = conjg ( aux( fc%nlt( ig )) ) * omega
      enddo
    enddo
    f_g(1:fc%ngmt, 1:fc%ngmt) = f_r(1:fc%ngmt,1:fc%ngmt)
else 
    call errore (' FFT routines',' Wrong switch',1)
end if
end subroutine fft6

END MODULE fft6_module

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
SUBROUTINE construct_w(scrcoul_g, scrcoul_pade_g, w_ryd)

  USE cell_base,     ONLY : omega
  USE constants,     ONLY : pi
  USE control_gw,    ONLY : eta, godbyneeds, padecont, modielec, do_imag
  USE disp,          ONLY : nq1, nq2, nq3
  USE freq_gw,       ONLY : fiu, nfs
  USE gwsigma,       ONLY : gcutcorr
  USE kinds,         ONLY : DP
  USE mp_global,     ONLY : mp_global_end
  USE timing_module, ONLY : time_construct_w

  implicit none

  complex(DP) :: scrcoul_pade_g (gcutcorr, gcutcorr)
  complex(DP) :: z(nfs), a(nfs)
  complex(DP)  :: scrcoul_g    (gcutcorr, gcutcorr, nfs) 

  real(DP) :: w_ryd
  real(DP) :: rcut

  integer :: ig, igp
  integer  :: iwim

  CALL start_clock(time_construct_w)

   rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
   scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)
   if(.NOT.modielec) then
     do ig = 1, gcutcorr
        do igp = 1, gcutcorr
           do iwim = 1, nfs
               z(iwim) = fiu(iwim)
               a(iwim) = scrcoul_g (ig,igp,iwim)
           enddo
           if (padecont) then
               call pade_eval ( nfs, z, a, cmplx(0.0_dp, w_ryd, kind=dp), scrcoul_pade_g (ig,igp))
           else if (godbyneeds .and. do_imag) then
               scrcoul_pade_g(ig,igp) = -a(2)/(w_ryd**2+(a(1))**2)
           else if (godbyneeds .and. .not. do_imag) then
               scrcoul_pade_g(ig,igp) = a(2)/(cmplx(w_ryd**2, 0.0_dp, kind=dp)-(a(1)-(0.0d0,1.0d0)*eta)**2)
           else
                write(6,'("No screening model chosen!")')
                stop
                call mp_global_end()
           endif
        enddo
     enddo
   else
     CALL errore(__FILE__, "modielec not maintained anymore", 1)
   endif

   CALL stop_clock(time_construct_w)

end SUBROUTINE construct_w

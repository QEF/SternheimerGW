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
SUBROUTINE coulpade(scrcoul_g, xq_ibk, vcut)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, modielec, truncation
  USE freq_gw,       ONLY : fiu, nfs
  USE gwsigma,       ONLY : sigma_c_st, gcutcorr
  USE gvect,         ONLY : g, ngm, nl
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE lr_symm_base,  ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE truncation_module, ONLY : truncate, vcut_type

  IMPLICIT NONE

  !> the truncated Coulomb potential
  TYPE(vcut_type), INTENT(IN) :: vcut

  complex(DP) ::  scrcoul_g   (gcutcorr, gcutcorr, nfs)
  complex(DP) :: z(nfs), u(nfs), a(nfs)
  complex(DP) :: phase
  complex(DP) :: eigv     (ngm, nrot)  

  real(DP) :: qg2, qg, qxy, qz, q_G(3)
  real(DP) :: rcut, spal, zcut, factor
  real(DP) :: xq_ibk(3), xq_ibz(3)

  integer :: gmapsym  (ngm, nrot) 
  integer :: ig, igp, irr, icounter, ir, irp
  integer :: iwim, iw, ikq
  integer :: iqstart, iqstop, iqs, nkr
  integer :: iq, ipol, iqrec, isym

  logical :: pade_catch
  logical :: found_q
  logical :: limq, inv_q, found

!Rotate G_vectors for FFT.
   rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
   if(.not.modielec) THEN
       do iw = 1, nfs
         do ig = 1, gcutcorr
            q_G = tpiba * (g(:,ig) + xq_ibk)
            factor = truncate(truncation, vcut, q_G)
            do igp = 1, gcutcorr
              scrcoul_g(ig, igp, iw) = scrcoul_g(ig, igp, iw) * factor 
            end do
         enddo!ig
       enddo!nfs
  endif
    if(.not.modielec) THEN
        if(godbyneeds) THEN
          do ig = 1, gcutcorr
            do igp = 1, gcutcorr 
!For godby-needs plasmon pole the algebra is done assuming real frequency*i.
!that is: the calculation is done at i*wp but we pass a real number as the freq.
               do iw = 1, nfs
                  z(iw) = dcmplx(aimag(fiu(iw)), 0.0d0)
                  u(iw) = scrcoul_g(ig, igp, iw)
               enddo
               call godby_needs_coeffs(nfs, z, u, a)
               do iw = 1, nfs 
!Just overwrite scrcoul_g with godby-needs coefficients.
                  scrcoul_g (ig, igp, iw) = a(iw)
               enddo
          enddo
         enddo
       else if (padecont) THEN
         do igp = 1, gcutcorr
          do ig = 1, gcutcorr
!Pade input points on the imaginary axis
             do iw = 1, nfs
                z(iw) = fiu(iw)
                u(iw) = scrcoul_g (ig, igp, iw)
             enddo
             call pade_coeff ( nfs, z, u, a)
!Overwrite scrcoul with Pade coefficients to be passed to pade_eval.
             do iw = 1, nfs 
                scrcoul_g (ig, igp, iw) = a(iw)
             enddo
          enddo !enddo on ig
       enddo  !enddo on igp
       else if(.not.padecont.and..not.godbyneeds) THEN
                 WRITE(6,'("No screening model chosen!")')
       endif
    endif
end SUBROUTINE coulpade

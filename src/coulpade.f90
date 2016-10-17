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
SUBROUTINE coulpade(num_g_corr, num_gp_corr, scrcoul_g, xq_ibk, vcut)

  USE cell_base,         ONLY : tpiba
  USE control_gw,        ONLY : godbyneeds, padecont, modielec, truncation
  USE freq_gw,           ONLY : fiu, nfs
  USE gvect,             ONLY : g
  USE kinds,             ONLY : DP
  USE truncation_module, ONLY : truncate, vcut_type

  IMPLICIT NONE

  !> the number of G vectors in the correlation grid
  INTEGER, INTENT(IN) :: num_g_corr

  !> the number of G' vectors in the correlation grid
  INTEGER, INTENT(IN) :: num_gp_corr

  !> the truncated Coulomb potential
  TYPE(vcut_type), INTENT(IN) :: vcut

  complex(DP) ::  scrcoul_g   (num_g_corr, num_gp_corr, nfs)
  complex(DP) :: z(nfs), u(nfs), a(nfs)

  real(DP) :: q_G(3)
  real(DP) :: factor
  real(DP) :: xq_ibk(3)

  integer :: ig, igp
  integer :: iw

!Rotate G_vectors for FFT.
   if(.not.modielec) THEN
       do iw = 1, nfs
         do ig = 1, num_g_corr
            q_G = tpiba * (g(:,ig) + xq_ibk)
            factor = truncate(truncation, vcut, q_G)
            do igp = 1, num_gp_corr
              scrcoul_g(ig, igp, iw) = scrcoul_g(ig, igp, iw) * factor 
            end do
         enddo!ig
       enddo!nfs
  endif
    if(.not.modielec) THEN
        if(godbyneeds) THEN
          do ig = 1, num_g_corr
            do igp = 1, num_gp_corr 
!For godby-needs plasmon pole the algebra is done assuming real frequency*i.
!that is: the calculation is done at i*wp but we pass a real number as the freq.
               do iw = 1, nfs
                  z(iw) = cmplx(aimag(fiu(iw)), 0.0_dp, kind=dp)
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
         do igp = 1, num_g_corr
          do ig = 1, num_gp_corr
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

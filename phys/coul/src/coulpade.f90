!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
!> Evaluate the convolution of \f$\epsilon^{-1}\f$ and Coulomb potential. 
MODULE coulpade_module

  IMPLICIT NONE

CONTAINS

!> Evaluate the convolution of \f$\epsilon^{-1}\f$ and Coulomb potential. 
!!
!! This routine takes the inverse of the dielectric function
!! \f$\epsilon^{-1}(G, G', \omega)\f$ and multiplies it with the bare Coulomb
!! potential to obtain the screened Coulomb interaction.
SUBROUTINE coulpade(xq_ibk, vcut, scrcoul_g)

  USE cell_base,          ONLY : tpiba
  USE control_gw,         ONLY : truncation
  USE gvect,              ONLY : g
  USE kinds,              ONLY : DP
  USE truncation_module,  ONLY : truncate, vcut_type

  IMPLICIT NONE

  !> the current q vector
  REAL(dp), INTENT(IN) :: xq_ibk(3)

  !> the truncated Coulomb potential
  TYPE(vcut_type), INTENT(IN) :: vcut

  !> *on input*: the inverse of the dielectric constant <br>
  !! *on output*: the coefficients used to evaluate the screened Coulomb
  !! interaction at an arbitrary frequency
  COMPLEX(dp), INTENT(INOUT) ::  scrcoul_g(:, :, :)

  !> the number of G vectors in the correlation grid
  INTEGER :: num_g_corr

  !> the vector q + G
  REAL(dp) :: q_G(3)

  !> the strength of the truncated Coulomb potential
  REAL(dp) :: factor

  !> loop variables for G and G'
  INTEGER :: ig, igp

  !> total number of frequencies
  INTEGER :: num_freq

  !> loop variable for the frequency
  INTEGER :: ifreq

  ! initialize helper variable
  num_g_corr = SIZE(scrcoul_g, 1)
  num_freq = SIZE(scrcoul_g, 3)

  ! sanity check for the array size
  IF (SIZE(scrcoul_g, 2) /= num_g_corr) &
    CALL errore(__FILE__, "input array should have same dimension for G and G'", 1)

  !                 -1
  ! evaluate W = eps  V
  !
  ! loop over frequencies
  DO ifreq = 1, num_freq

    ! outer loop over G
    DO ig = 1, num_g_corr

       ! determine V at q + G'
       q_G = tpiba * (g(:,ig) + xq_ibk)
       factor = truncate(truncation, vcut, q_G)

       ! inner loop over G'
       DO igp = 1, num_g_corr
         scrcoul_g(ig, igp, ifreq) = scrcoul_g(ig, igp, ifreq) * factor
       END DO ! igp

    END DO ! ig

  END DO ! ifreq

END SUBROUTINE coulpade

END MODULE coulpade_module

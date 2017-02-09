!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2017
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
!> Provide the routines to evaluate the analytic continuation.
MODULE coulpade_module

  IMPLICIT NONE

CONTAINS

!> Evaluate the analytic continuation of the screened Coulomb interaction. 
!!
!! This routine takes the inverse of the dielectric function
!! \f$\epsilon^{-1}(G, G', \omega)\f$ and multiplies it with the bare Coulomb
!! potential to obtain the screened Coulomb interaction. Then it evaluates the
!! coefficients to perform an analytic continuation to the full complex plane
!! from the given frequency mesh.
!!
!! There are a few different methods implemented to evaluate the analytic
!! continuation
SUBROUTINE coulpade(xq_ibk, freq, vcut, scrcoul_g)

  USE cell_base,          ONLY : tpiba
  USE control_gw,         ONLY : godbyneeds, padecont, paderobust, truncation, tr2_gw
  USE freqbins_module,    ONLY : freqbins_type, freqbins_symm
  USE godby_needs_module, ONLY : godby_needs_coeffs
  USE gvect,              ONLY : g
  USE kinds,              ONLY : DP
  USE pade_module,        ONLY : pade_coeff_robust
  USE truncation_module,  ONLY : truncate, vcut_type

  IMPLICIT NONE

  !> the current q vector
  REAL(dp), INTENT(IN) :: xq_ibk(3)

  !> the frequency grid used for the calculation
  TYPE(freqbins_type), INTENT(IN) :: freq

  !> the truncated Coulomb potential
  TYPE(vcut_type), INTENT(IN) :: vcut

  !> *on input*: the inverse of the dielectric constant <br>
  !! *on output*: the coefficients used to evaluate the screened Coulomb
  !! interaction at an arbitrary frequency
  COMPLEX(dp), INTENT(INOUT) ::  scrcoul_g(:, :, :)

  !> the number of G vectors in the correlation grid
  INTEGER :: num_g_corr

  !> frequency used for Pade coefficient (will be extended if frequency
  !! symmetry is used)
  COMPLEX(dp), ALLOCATABLE :: z(:)

  !> value of the screened Coulomb interaction on input mesh
  COMPLEX(dp), ALLOCATABLE :: u(:)

  !> coefficients of the Pade approximation
  COMPLEX(dp), ALLOCATABLE :: a(:)

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
  num_freq = SIZE(freq%solver)

  ! sanity check for the array size
  num_g_corr = SIZE(scrcoul_g, 1)
  IF (SIZE(scrcoul_g, 2) /= num_g_corr) &
    CALL errore(__FILE__, "input array should have same dimension for G and G'", 1)
  IF (SIZE(scrcoul_g, 3) /= freq%num_freq()) &
    CALL errore(__FILE__, "frequency dimension of Coulomb inconsistent with frequency mesh", 1)

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

  !
  ! analytic continuation to the complex plane
  !
  !! 1. Godby-Needs plasmon-pole model - assumes that the function can be accurately
  !!    represented by a single pole and uses the value of the function at two
  !!    frequencies \f$\omega = 0\f$ and \f$\omega = \omega_{\text{p}}\f$ to determine
  !!    the parameters.
  IF (godbyneeds) THEN
    CALL godby_needs_coeffs(AIMAG(freq%solver(2)), scrcoul_g)

  !! 2. Pade expansion - evaluate Pade coefficients for a continued fraction expansion
  !!    using a given frequency grid; symmetry may be used to extend the frequency grid
  !!    to more points.
  ELSE IF (padecont) THEN

    ! allocate helper arrays
    ALLOCATE(u(freq%num_freq()))
    ALLOCATE(a(freq%num_freq()))

    ! use symmetry to extend the frequency mesh
    CALL freqbins_symm(freq, z, scrcoul_g)

    ! evalute Pade approximation for all G and G'
    DO igp = 1, num_g_corr
     DO ig = 1, num_g_corr

       ! set frequency and value used to determine the Pade coefficients
       u = scrcoul_g(ig, igp, :)

       ! evaluate the coefficients
       CALL pade_coeff(freq%num_freq(), z, u, a)

       ! store the coefficients in the same array
       scrcoul_g(ig, igp, :) = a

     ENDDO ! ig
  ENDDO ! igp

  !! 3. robust Pade expansion - evaluate Pade coefficients using a circular frequency
  !!    mesh in the complex plane
  ELSEIF (paderobust) THEN
    CALL pade_coeff_robust(freq%solver, tr2_gw, scrcoul_g)
  ELSE
    CALL errore(__FILE__, "No screening model chosen!", 1)
  END IF

END SUBROUTINE coulpade

END MODULE coulpade_module

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
!> Evaluate the analytic continuation of the screened Coulomb interaction.
!!
MODULE construct_w_module

  IMPLICIT NONE

CONTAINS

  !> Construct the screened Coulomb interaction for an arbitrary frequency.
  SUBROUTINE construct_w(gmapsym, grid, freq_in, scrcoul_coeff, freq_out, scrcoul)

    USE control_gw,         ONLY : godbyneeds, padecont, paderobust
    USE freqbins_module,    ONLY : freqbins_type
    USE godby_needs_module, ONLY : godby_needs_model
    USE kinds,              ONLY : dp
    USE pade_module,        ONLY : pade_eval_robust
    USE sigma_grid_module,  ONLY : sigma_grid_type
    USE timing_module,      ONLY : time_construct_w

    !> The symmetry map from the irreducible point to the current one
    INTEGER,                  INTENT(IN)  :: gmapsym(:)

    !> the FFT grids on which the screened Coulomb interaction is evaluated
    TYPE(sigma_grid_type),    INTENT(IN)  :: grid

    !> the frequency grid on which W was evaluated
    TYPE(freqbins_type),      INTENT(IN)  :: freq_in

    !> the coefficients of the screened Coulomb potential used for the analytic continuation
    COMPLEX(dp),              INTENT(IN)  :: scrcoul_coeff(:,:,:)

    !> the frequency for which the screened Coulomb potential is evaluated
    COMPLEX(dp),              INTENT(IN)  :: freq_out

    !> The screened Coulomb interaction symmetry transformed and parallelized over images.
    !! The array is appropriately sized to do a FFT on the output.
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: scrcoul(:,:)

    !> Counter on the G and G' vector
    INTEGER ig, igp

    !> corresponding point to G' in global G list
    INTEGER igp_g

    !> allocation error flag
    INTEGER ierr

    !> helper array to extract the current coefficients
    COMPLEX(dp), ALLOCATABLE :: coeff(:)

    !> complex constant of zero
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND = dp)

    CALL start_clock(time_construct_w)

    !
    ! create and initialize output array
    ! allocate space so that we can perform an in-place FFT on the array
    !
    ALLOCATE(scrcoul(grid%corr%dfftt%nnr, grid%corr_par%dfftt%nnr), STAT = ierr)
    IF (ierr /= 0) THEN
      CALL errore(__FILE__, "allocation of screened Coulomb potential failed", 1)
      RETURN
    END IF
    scrcoul = zero

    !
    ! construct screened Coulomb interaction
    !
    !! The screened Coulomb interaction is interpolated with either Pade or
    !! Godby-Needs analytic continuation. We only evaluate W at the irreducible
    !! mesh, but any other point may be obtained by
    !! \f{equation}{
    !!   W_{S q}(G, G') = W_{q}(S^{-1} G, S^{-1} G')~.
    !! \f}
    ALLOCATE(coeff(freq_in%num_freq()))
    DO igp = 1, grid%corr_par%ngmt
      !
      ! get the global corresponding index
      igp_g = grid%corr_par%ig_l2gt(igp)

      DO ig = 1, grid%corr%ngmt

        ! symmetry transformation of the coefficients
        coeff = scrcoul_coeff(gmapsym(ig), gmapsym(igp_g), :)

        IF (padecont) THEN
          !
          ! Pade analytic continuation
          CALL pade_eval(freq_in%num_freq(), freq_in%solver, coeff, freq_out, scrcoul(ig, igp))

        ELSE IF (paderobust) THEN
          !
          ! robust Pade analytic continuation
          CALL pade_eval_robust(coeff, freq_out, scrcoul(ig, igp))

        ELSE IF (godbyneeds) THEN
          !
          ! Godby-Needs Pole model
          scrcoul(ig, igp) = godby_needs_model(freq_out, coeff)

        ELSE
          CALL errore(__FILE__, "No screening model chosen!", 1)

        END IF

      END DO ! ig
    END DO ! igp

    CALL stop_clock(time_construct_w)

  END SUBROUTINE construct_w

END MODULE construct_w_module

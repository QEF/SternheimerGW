!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2017
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
!> This module can be used to plot the screened Coulomb interaction along the
!! real frequency axis.
MODULE plot_coulomb_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC plot_coulomb

CONTAINS

  !> Main driver for plotting routines that calls the appropriate actual
  !! plotting routine
  SUBROUTINE plot_coulomb(freq, coulomb)

    USE analytic_module,    ONLY: godby_needs, pade_approx
    USE control_gw,         ONLY: model_coul
    USE freqbins_module,    ONLY: freqbins_type
    USE kinds,              ONLY: dp

    !> the frequency grid used for the calculation
    TYPE(freqbins_type), INTENT(IN) :: freq

    !> the screened coulomb interaction for this frequency
    COMPLEX(dp),         INTENT(IN) :: coulomb(:,:,:)

    SELECT CASE (model_coul)
    CASE (godby_needs)
      CALL plot_coulomb_godbyneeds(freq, coulomb)
    CASE (pade_approx)
      CALL plot_coulomb_pade(freq, coulomb)
    END SELECT

  END SUBROUTINE plot_coulomb

  !> Plot the screened Coulomb interaction along the real frequency axis for a
  !! given frequency mesh using the Pade approximation.
  SUBROUTINE plot_coulomb_pade(freq, coulomb)

    USE constants,          ONLY: RYTOEV
    USE freqbins_module,    ONLY: freqbins_type, freqbins_symm
    USE io_global,          ONLY: stdout
    USE kinds,              ONLY: dp

    !> the frequency grid used for the calculation
    TYPE(freqbins_type), INTENT(IN) :: freq

    !> the screened coulomb interaction for this frequency
    COMPLEX(dp),         INTENT(IN) :: coulomb(:,:,:)

    !> number of G vectors in correlation grid
    INTEGER num_g_corr

    !> counter on the G vectors
    INTEGER ig, igp

    !> the frequency grid used
    COMPLEX(dp), ALLOCATABLE :: omega_in(:)

    !> the frequency at which the Pade approximation is evaluated
    COMPLEX(dp) :: omega_out

    !> counter on the frequencies
    INTEGER ifreq

    !> Pade coefficients
    COMPLEX(dp), ALLOCATABLE :: coeff(:)

    !> symmetrized version of the coulomb array
    COMPLEX(dp), ALLOCATABLE :: coulomb_sym(:,:,:)

    !> Pade approximant to coulomb at a specific frequency
    COMPLEX(dp) pade

    num_g_corr = SIZE(coulomb, 1)
    IF (SIZE(coulomb, 2) /= num_g_corr) &
      CALL errore(__FILE__, "coulomb matrix should be square for each frequency", 1)
    IF (SIZE(coulomb, 3) > freq%num_freq()) &
      CALL errore(__FILE__, "coulomb matrix has more frequencies than frequency type", 1)

    ! copy coulomb to array for symmetrization
    ALLOCATE(coulomb_sym(num_g_corr, num_g_corr, freq%num_freq()))
    coulomb_sym(:,:,:SIZE(coulomb, 3)) = coulomb

    ! use symmetry to extend the frequency mesh
    CALL freqbins_symm(freq, omega_in, coulomb_sym)

    ! allocate arrays to contain Pade coefficients
    ALLOCATE(coeff(freq%num_freq()))

    WRITE(stdout, '(5x,a)') 'Plotting Pade approximation'

    ! evalute Pade approximation for all G and G'
    DO igp = 1, num_g_corr
      DO ig = 1, num_g_corr

        ! evaluate the coefficients
        CALL pade_coeff(freq%num_freq(), omega_in, coulomb_sym(ig, igp, :), coeff)
 
        ! plot the frequency dependent Coulomb interaction
        WRITE(stdout, '(17x,a,18x,a,i6,a,i6,a)') 'omega', 'W(', ig, ', ', igp, ', omega)'

        DO ifreq = 1, freq%num_coul()

          ! evaluate Pade approximation
          omega_out = freq%coul(ifreq)
          CALL pade_eval(freq%num_freq(), omega_in, coeff, omega_out, pade)

          WRITE(stdout, '(5x,2f15.8,2x,2f15.8)') omega_out * RYTOEV, pade

        END DO ! ifreq

        WRITE(stdout, '(a)')
        
      END DO ! ig
    END DO ! igp

    WRITE(stdout, '(5x,a)') 'End of Pade approximation'
    WRITE(stdout, '(a)')

  END SUBROUTINE plot_coulomb_pade

  !> Plot the screened Coulomb interaction along the real frequency axis for a
  !! given frequency mesh using the Godby-Needs plasmon-pole model.
  SUBROUTINE plot_coulomb_godbyneeds(freq, coulomb)

    USE constants,          ONLY: RYTOEV
    USE freqbins_module,    ONLY: freqbins_type, freqbins_symm
    USE godby_needs_module, ONLY: godby_needs_coeffs, godby_needs_model
    USE io_global,          ONLY: stdout
    USE kinds,              ONLY: dp

    !> the frequency grid used for the calculation
    TYPE(freqbins_type), INTENT(IN) :: freq

    !> the screened coulomb interaction for this frequency
    COMPLEX(dp),         INTENT(IN) :: coulomb(:,:,:)

    !> number of G vectors in correlation grid
    INTEGER num_g_corr

    !> counter on the G vectors
    INTEGER ig, igp

    !> the frequency at which the Pade approximation is evaluated
    COMPLEX(dp) :: omega_out

    !> counter on the frequencies
    INTEGER ifreq

    !> copy of the coulomb array
    COMPLEX(dp), ALLOCATABLE :: coulomb_(:,:,:)

    num_g_corr = SIZE(coulomb, 1)
    IF (SIZE(coulomb, 2) /= num_g_corr) &
      CALL errore(__FILE__, "coulomb matrix should be square for each frequency", 1)
    IF (SIZE(coulomb, 3) > 2) &
      CALL errore(__FILE__, "coulomb matrix should have two frequencies for PP model", 1)

    ! copy coulomb to array
    ALLOCATE(coulomb_(num_g_corr, num_g_corr, 2))
    coulomb_ = coulomb

    ! calculate Godby-Needs coefficients (in-place replacement)
    CALL godby_needs_coeffs(AIMAG(freq%solver(2)), coulomb_) 

    WRITE(stdout, '(5x,a)') 'Plotting Godby-Needs approximation'

    ! evalute Pade approximation for all G and G'
    DO igp = 1, num_g_corr
      DO ig = 1, num_g_corr

        ! plot the frequency dependent Coulomb interaction
        WRITE(stdout, '(17x,a,18x,a,i6,a,i6,a)') 'omega', 'W(', ig, ', ', igp, ', omega)'

        DO ifreq = 1, freq%num_coul()

          ! evaluate Godby-Needs approximation
          omega_out = freq%coul(ifreq)

          WRITE(stdout, '(5x,2f15.8,2x,2f15.8)') omega_out * RYTOEV, &
            godby_needs_model(omega_out, coulomb_(ig, igp, :))

        END DO ! ifreq

        WRITE(stdout, '(a)')
        
      END DO ! ig
    END DO ! igp

    WRITE(stdout, '(5x,a)') 'End of Godby-Needs approximaton'
    WRITE(stdout, '(a)')

  END SUBROUTINE plot_coulomb_godbyneeds

END MODULE plot_coulomb_module

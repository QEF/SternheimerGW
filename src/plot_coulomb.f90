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
!> This module can be used to plot the screened Coulomb interaction along the
!! real frequency axis.
MODULE plot_coulomb_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC plot_coulomb

CONTAINS

  !> Plot the screened Coulomb interaction for a given frequency grid.
  !!
  !! Note that this routines uses the same routines to perform the analytic
  !! continuation than the code will use later on. In this way the plotted
  !! function provide an actual test of the values that you will obtain later.
  SUBROUTINE plot_coulomb(model_coul, thres, grid, freq, coulomb)

    USE analytic_module,    ONLY: analytic_coeff, analytic_eval
    USE constants,          ONLY: RYTOEV
    USE freqbins_module,    ONLY: freqbins_type
    USE io_global,          ONLY: stdout
    USE kinds,              ONLY: dp
    USE sigma_grid_module,  ONLY: sigma_grid_type

    !> the screening model used for the analytic continuation
    INTEGER,               INTENT(IN) :: model_coul

    !> the threshold used for some of the screening models
    REAL(dp),              INTENT(IN) :: thres

    !> the FFT grids on which the screened Coulomb interaction is evaluated
    TYPE(sigma_grid_type), INTENT(IN) :: grid

    !> the frequency grid used for the calculation
    TYPE(freqbins_type),   INTENT(IN) :: freq

    !> the screened coulomb interaction for these frequencies
    COMPLEX(dp),           INTENT(IN) :: coulomb(:,:,:)

    !> copy of the screened coulomb interaction
    COMPLEX(dp), ALLOCATABLE :: scrcoul_in(:,:,:)

    !> result of the analytic continuation
    COMPLEX(dp), ALLOCATABLE :: scrcoul_out(:,:,:)

    !> work array
    COMPLEX(dp), ALLOCATABLE :: work(:,:)

    !> symmetry map between G order (we use the identity)
    INTEGER, ALLOCATABLE :: gmapsym(:)

    !> number of G vectors
    INTEGER num_g

    !> loop variables for G vectors
    INTEGER ig, igp

    !> loop variable for frequencies
    INTEGER ifreq

    !
    ! sanity test of the input
    !
    num_g = SIZE(coulomb, 1)
    IF (SIZE(coulomb, 2) /= num_g) &
      CALL errore(__FILE__, "coulomb must be square matrix", 1)
    IF (grid%corr_fft%ngm /= num_g) &
      CALL errore(__FILE__, "FFT grid and coulomb matrix inconsistent", 1)
    IF (grid%corr_fft%ngm /= grid%corr_par_fft%ngm) &
      CALL errore(__FILE__, "image parallelism not implemented for plotting", 1)

    !
    ! evaluate coefficients for analytic continuation
    !
    ! create a copy of coulomb
    ALLOCATE(scrcoul_in(num_g, num_g, freq%num_freq()))
    scrcoul_in(:, :, 1:SIZE(coulomb, 3)) = coulomb

    ! evaluate coefficients for screening model
    CALL analytic_coeff(model_coul, thres, freq, scrcoul_in)

    !
    ! perform analytic continuation to the desired output mesh
    !
    ! create array for the results
    ALLOCATE(scrcoul_out(num_g, num_g, freq%num_coul()))

    ! create symmetry map
    ALLOCATE(gmapsym(num_g))
    gmapsym = [(ig, ig = 1, num_g)]

    ! perform analytic continuation
    DO ifreq = 1, freq%num_coul()
      CALL analytic_eval(gmapsym, grid, freq, scrcoul_in, freq%coul(ifreq), work)
      scrcoul_out(:,:,ifreq) = work
    END DO ! ifreq

    !
    ! print the results
    !
    WRITE(stdout, '(5x,a)') 'Plotting analytic continuation'

    ! loop over all elements of matrix
    DO igp = 1, num_g
      DO ig = 1, num_g

        ! plot the frequency dependent Coulomb interaction
        WRITE(stdout, '(17x,a,18x,a,i6,a,i6,a)') 'omega', 'W(', ig, ', ', igp, ', omega)'
        DO ifreq = 1, freq%num_coul()
          WRITE(stdout, '(5x,2f15.8,2x,2f15.8)') freq%coul(ifreq) * RYTOEV, scrcoul_out(ig, igp, ifreq)
        END DO ! ifreq
        WRITE(stdout, '(a)')

      END DO ! ig
    END DO ! igp

    WRITE(stdout, '(5x,a)') 'End of analytic continuation'
    WRITE(stdout, '(a)')
    
  END SUBROUTINE plot_coulomb

END MODULE plot_coulomb_module

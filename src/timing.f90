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
!> This module specifies how the timing of different parts of the SGW code
!! is measured.
!!
!! The timing follows the general layout of the code
!! # initializing all required quantities
!! # evaluate the screened Coulomb potential W
!! # evaluate the electronic Green's function G
!! # convolute G and W to obtain correlation part of Sigma
!! # evaluate the exchange part of Sigma
!! # evaluate matrix element of Sigma with the wave functions
!!
MODULE timing_module

  IMPLICIT NONE

  !> label for the clock measuring the initialization time
  CHARACTER(*), PARAMETER :: time_setup = 'setup'

  !> label for the clock measuring the calculation of W
  CHARACTER(*), PARAMETER :: time_coulomb = 'coulomb'

  !> label for the clock measuring the calculation of G
  CHARACTER(*), PARAMETER :: time_green = 'green'

  !> label for the clock measuring the convolution of G and W
  CHARACTER(*), PARAMETER :: time_sigma_c = 'sigma_c'

  !> label for the clock measuring the exchange part of Sigma
  CHARACTER(*), PARAMETER :: time_sigma_x = 'sigma_x'

  !> label for the clock measuring the evaluation of the matrix elements
  CHARACTER(*), PARAMETER :: time_matel = 'mat_el'

CONTAINS

  !> Print the measure timing of the SGW run in a nice format.
  !!
  !!
  SUBROUTINE timing_print_clock()

    USE io_global, ONLY: stdout

    ! empty line to separate it from the rest of the output
    WRITE(stdout, *)

    ! info line
    WRITE(stdout, *) 'Timing of the code:'

    !
    ! Overview over the different parts
    !

    ! print the time needed for the setup
    CALL print_clock(time_setup)

    ! print the time needed to calculate W
    CALL print_clock(time_coulomb)

  END SUBROUTINE timing_print_clock

END MODULE timing_module

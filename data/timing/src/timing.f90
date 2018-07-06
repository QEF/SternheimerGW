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
!> This module specifies how the timing of different parts of the SternheimerGW code
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
  CHARACTER(*), PARAMETER :: time_green = 'Green'

  !> label for the clock measuring the convolution of G and W
  CHARACTER(*), PARAMETER :: time_sigma_c = 'Sigma_c'

  !> label for the clock measuring the exchange part of Sigma
  CHARACTER(*), PARAMETER :: time_sigma_x = 'Sigma_x'

  !> label for the clock measuring the evaluation of the matrix elements
  CHARACTER(*), PARAMETER :: time_matel = 'Sigma matel'

  !
  ! split screened Coulomb into parts
  !

  !> label for the clock measuring the nscf calculation (W part)
  CHARACTER(*), PARAMETER :: time_coul_nscf = 'coul nscf'

  !> label for the clock measuring the linear solver (W part)
  CHARACTER(*), PARAMETER :: time_coul_solver = 'coul solver'

  !> label for the clock measuring the matrix inversion (W part)
  CHARACTER(*), PARAMETER :: time_coul_invert = 'coul invert'

  !> label for the clock measuring the communication (W part)
  CHARACTER(*), PARAMETER :: time_coul_comm = 'coul communication'

  !> label for the clock measuring the file I/O (W part)
  CHARACTER(*), PARAMETER :: time_coul_io = 'coul IO'

  !> label for the clock measuring setting up the symmetry
  CHARACTER(*), PARAMETER :: time_coul_symm = 'coul symm'

  !> label for the clock measuring unfolding the matrix
  CHARACTER(*), PARAMETER :: time_coul_unfold = 'coul unfold'

  !
  ! split Sigma_c into parts
  !

  !> label for the clock measuring the time to setup quantities for Sigma
  CHARACTER(*), PARAMETER :: time_sigma_setup = 'setup sigma'

  !> label for the clock measuring the time to construct W in real frequency
  CHARACTER(*), PARAMETER :: time_construct_w = 'construct W'

  !> label for the clock measuring the time to multiply G and W
  CHARACTER(*), PARAMETER :: time_GW_product = 'Sigma = G*W'

  !> label for the clock measuring the communication time
  CHARACTER(*), PARAMETER :: time_sigma_comm = 'Sigma communication'

  !> label for the clock measuring the time to store Sigma to disk
  CHARACTER(*), PARAMETER :: time_sigma_io = 'Sigma IO'

  !
  ! general routines
  !
  !> label for the total time spent to evaluate the linear operator
  CHARACTER(*), PARAMETER :: time_linear_op = 'linear operator'

  !> label for the total time spent calculating forward 6d FFTs
  CHARACTER(*), PARAMETER :: time_fwfft6 = '6d FFT r->G'

  !> label for the total time spent calculating inverse 6d FFTs
  CHARACTER(*), PARAMETER :: time_invfft6 = '6d FFT G->r'

CONTAINS

  !> Print the measured timing of the SternheimerGW run in a nice format.
  !!
  !! First print a summary of the main parts
  !! -setup
  !! -calculation of W
  !! -calculation of G
  !! -convolution of G and W
  !! -exchange
  !! -matrix elements of Sigma
  !! then resolve these parts into the important contributions.
  !!
  SUBROUTINE timing_print_clock()

    USE io_global, ONLY: stdout

    ! empty line to separate it from the rest of the output
    WRITE(stdout, *)

    ! info line
    WRITE(stdout,'(a)') 'Timing of the code:'

    !
    ! Overview over the different parts
    !

    ! print the time needed for the setup
    CALL print_clock(time_setup)

    ! print the time needed to calculate W
    CALL print_clock(time_coulomb)

    ! print the time needed to calculate Sigma_c
    CALL print_clock(time_sigma_c)

    ! print the time needed to calculate Sigma_x
    CALL print_clock(time_sigma_x)

    ! print the time needed to evaluate the matrix elements
    CALL print_clock(time_matel)

    !
    ! Detailed part of screened Coulomb interaction
    !
    ! empty line to separate it from the rest of the output
    WRITE(stdout, *)

    ! info line
    WRITE(stdout,'(a)') 'Needed for screened Coulomb interaction'

    ! print the time needed for the nscf calculations
    CALL print_clock(time_coul_nscf)

    ! print the time needed for the linear solver
    CALL print_clock(time_coul_solver)

    ! print the time needed for the matrix inversion
    CALL print_clock(time_coul_invert)

    ! print the time needed for the communication
    CALL print_clock(time_coul_comm)

    ! print the time needed for the file I/O
    CALL print_clock(time_coul_io)

    ! print the time needed for setting up the symmetry
    CALL print_clock(time_coul_symm)

    ! print the time needed for unfolding the matrix
    CALL print_clock(time_coul_unfold)

    !
    ! Detailed part of correlation part of Sigma
    !
    ! empty line to separate it from the rest of the output
    WRITE(stdout, *)

    ! info line
    WRITE(stdout,'(a)') "Needed for correlation part of Sigma"

    ! print the time needed to setup things for Sigma
    CALL print_clock(time_sigma_setup)

    ! print the time needed to construct W on real frequency mesh
    CALL print_clock(time_construct_w)

    ! print the time needed to calculate G
    CALL print_clock(time_green)

    ! print the time necessary to convolute G and W
    CALL print_clock(time_GW_product)

    ! print the time needed to communicate Sigma across the processes
    CALL print_clock(time_sigma_comm)

    ! print the time needed to write Sigma to disk
    CALL print_clock(time_sigma_io)

    !
    ! Some general routines
    !
    ! empty line to separate it from the rest of the output
    WRITE(stdout, *)

    ! info line
    WRITE(stdout,'(a)') "Timing for some important routines"

    ! print the time spent to evaluate the linear operator
    CALL print_clock(time_linear_op)

    ! print the time spent in the FFT routines
    CALL print_clock(time_fwfft6)
    CALL print_clock(time_invfft6)

  END SUBROUTINE timing_print_clock

END MODULE timing_module

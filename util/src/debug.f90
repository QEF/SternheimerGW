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
!> This module provides various utility routines to facilitate the debugging
!! of the code.
MODULE debug_module

  IMPLICIT NONE

#if defined (__DEBUG)
  LOGICAL, PARAMETER :: debug_set = .TRUE.
#else
  LOGICAL, PARAMETER :: debug_set = .FALSE.
#endif

  !> This type defines the options of the code which can be debugged.
  !!
  !! @note when introducing additional flags the default should be false
  TYPE debug_type

    !> debug linear solver for the Green's function
    LOGICAL :: solver_green = .FALSE.

    !> debug correlation self energy
    LOGICAL :: sigma_corr = .FALSE.

    !> the unit into which the debug log is written
    INTEGER :: note = 6

  END TYPE debug_type

  !> this is a wrapper for the ieee_arithmetic and/or gfortran's isnan
  !! it also extends the functionality to complex quantities
  INTERFACE test_nan
    MODULE PROCEDURE test_nan_real, test_nan_complex
  END INTERFACE test_nan

CONTAINS

  !> test if the input value is NaN
  ELEMENTAL FUNCTION test_nan_real(x) RESULT(is_nan)

    USE kinds, ONLY: dp

! gfortran only supports ieee_arithmetic since v5.0
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: isnan => ieee_is_nan
#endif

    !> the tested variable
    REAL(dp), INTENT(IN) :: x

    !> is the tested variable NaN
    LOGICAL is_nan

    is_nan = isnan(x)

  END FUNCTION test_nan_real

  !> test if the input value is NaN
  ELEMENTAL FUNCTION test_nan_complex(x) RESULT(is_nan)

    USE kinds, ONLY: dp

    !> the tested variable
    COMPLEX(dp), INTENT(IN) :: x

    !> is the tested variable NaN
    LOGICAL is_nan

    is_nan = test_nan_real(REAL(x)) .OR. test_nan_real(AIMAG(x))

  END FUNCTION test_nan_complex

  !> broadcast the debug type
  SUBROUTINE mp_bcast_debug(debug, source, comm)

    USE mp, ONLY: mp_bcast

    !> debug type to be communicated
    TYPE(debug_type), INTENT(INOUT) :: debug

    !> index of process containing the filled type
    INTEGER, INTENT(IN) :: source

    !> communicator across which the type is distributed
    INTEGER, INTENT(IN) :: comm

    CALL mp_bcast(debug%solver_green, source, comm)
    CALL mp_bcast(debug%sigma_corr,   source, comm)
    CALL mp_bcast(debug%note,         source, comm)

  END SUBROUTINE mp_bcast_debug

END MODULE debug_module

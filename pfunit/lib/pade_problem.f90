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
!> This module provides routines to test the Pade approximation
MODULE pade_problem_module

  IMPLICIT NONE

CONTAINS

  !> evaluate a function in the complex plane on a circle around the origin
  SUBROUTINE pade_problem_evaluate(radius, num_point, func, val)

    USE constants, ONLY: tpi
    USE kinds,     ONLY: dp

    !> radius of the circle in the complex plane
    REAL(dp), INTENT(IN) :: radius

    !> number of points constructed on the circle
    INTEGER,  INTENT(IN) :: num_point

    !> function pointer defining the function in the complex plane
    INTERFACE
      FUNCTION func(zz)
        USE kinds, ONLY: dp
        !> the position in the complex plane
        COMPLEX(dp), INTENT(IN) :: zz
        !> the function value at this point
        COMPLEX(dp) func
      END FUNCTION func
    END INTERFACE

    !> the values of the function on the circle
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: val(:)

    !> counter on the points in reciprocal space
    INTEGER ipoint

    !> phase factor between points
    COMPLEX(dp) phase

    !> position at which the function is evaluated
    COMPLEX(dp) zz

    !> complex constant of 2 * pi * i
    COMPLEX(dp), PARAMETER :: c2PiI = CMPLX(0.0_dp, tpi, KIND=dp)

    ! determine phase factor
    phase = EXP(c2PiI / num_point)

    ! initialize position
    zz = radius

    ! create output array
    ALLOCATE(val(num_point))
    !
    DO ipoint = 1, num_point
      !
      ! evaluate the function on the circle
      val(ipoint) = func(zz)
      !
      ! update position on the circle
      zz = zz * phase
      !
    END DO ! ipoint

  END SUBROUTINE pade_problem_evaluate

  !> wrapper of complex exponential function
  FUNCTION pade_exp(zz)

    USE kinds, ONLY: dp

    !> argument of the exponential function
    COMPLEX(dp), INTENT(IN) :: zz

    !> exponential of argument
    COMPLEX(dp) pade_exp

    pade_exp = EXP(zz)

  END FUNCTION pade_exp

  !> wrapper of complex cosine function
  FUNCTION pade_cos(zz)

    USE kinds, ONLY: dp

    !> argument of the cosine function
    COMPLEX(dp), INTENT(IN) :: zz

    !> cosine of the argument
    COMPLEX(dp) pade_cos

    pade_cos = COS(zz)

  END FUNCTION pade_cos

END MODULE pade_problem_module

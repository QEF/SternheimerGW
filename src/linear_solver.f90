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
!> This module provides a specialized linear solver for the problems that occur
!! in a Sternheimer \f$GW\f$ calculation.
!!
!! We need a solver capable of solving the linear equation \f$(A + \sigma I) x = b\f$.
!! Here, A is a linear operator, \f$\sigma\f$ is a shift, x is the solutions of the
!! problem, and b is the right-hand side of the problem. One can expand the operator
!! to matrices and solve for multiple right-hand sides at once.
MODULE linear_solver_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC linear_solver, linear_solver_config

  !> stores the configuration of the linear solver
  TYPE linear_solver_config

    !> the threshold when the calculation is considered converged
    REAL(dp) threshold

  END TYPE linear_solver_config

CONTAINS

  !! We use the following algorithm continuously updating the solution and the
  !! residual error of the solution.
  SUBROUTINE linear_solver(config, AA, sigma, BB, XX)

    !> configuration of the linear solver
    TYPE(linear_solver_config), INTENT(IN) :: config

    !> Function pointer that applies the linear operator to a vector.
    INTERFACE
      SUBROUTINE AA(sigma, xx, Ax)
        USE kinds, ONLY: dp
        !> The shift of this system.
        COMPLEX(dp), INTENT(IN)  :: sigma
        !> The input vector.
        COMPLEX(dp), INTENT(IN)  :: xx(:)
        !> The operator applied to the vector.
        COMPLEX(dp), INTENT(OUT) :: Ax(:)
      END SUBROUTINE AA
    END INTERFACE

    !> Right hand side of the linear equation.
    COMPLEX(dp), INTENT(IN)  :: BB(:,:)

    !> Shift \f$\sigma\f$ in the linear operator.
    COMPLEX(dp), INTENT(IN)  :: sigma(:)

    !> On output: the solution of the linear system
    COMPLEX(dp), INTENT(OUT) :: XX(:,:,:)

  END SUBROUTINE linear_solver

END MODULE linear_solver_module

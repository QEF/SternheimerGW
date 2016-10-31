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
!> Choose the linear solver suited for the linear problem.
!!
!! We solve multiple shifted linear systems
!! \f{equation}{
!!   (A + \sigma I) x = b
!! \f}
!! by reusing part of the Krylov subspace. If one linear solver struggles solving
!! the linear problem, we switch to another one, unless the user requests not to
!! do so.
MODULE select_solver_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  !> use the BiCGstab solver with multishift
  INTEGER, PARAMETER :: bicgstab_multi = 1

  !> use the BiCGstab solver without multishift
  INTEGER, PARAMETER :: bicgstab_no_multi = 2

  !> use the solver specialized for SGW
  INTEGER, PARAMETER :: sgw_linear_solver = 3

  !> configuration of the linear solver
  TYPE select_solver_type

    !> priority in which the linear solvers are chosen
    INTEGER, ALLOCATABLE :: priority(:)

  END TYPE select_solver_type

CONTAINS

  !> select the linear solver based on user preference and convergence
  SUBROUTINE select_solver(config, AA, bb, sigma, xx, ierr)

    !> configuration of the linear solver
    TYPE(select_solver_type), INTENT(IN) :: config

    !> Function pointer that applies the linear operator to a vector.
    INTERFACE
      SUBROUTINE AA(sigma, xx, Ax)
        USE kinds, ONLY: dp
        !> The shift of this system.
        COMPLEX(dp), INTENT(IN)  :: sigma
        !> the vector to which the linear operator is applied
        COMPLEX(dp), INTENT(IN)  :: xx(:)
        !! the vector after applying the linear operator
        COMPLEX(dp), INTENT(OUT) :: Ax(:)
      END SUBROUTINE AA
    END INTERFACE

    !> Right hand side of the linear equation.
    COMPLEX(dp), INTENT(IN)  :: bb(:)

    !> Shift \f$\sigma\f$ in the linear operator.
    COMPLEX(dp), INTENT(IN)  :: sigma(:)

    !> On output: the solution of the linear system
    COMPLEX(dp), INTENT(OUT) :: xx(:,:)

    !> the error code - should be zero for successful convergence
    INTEGER,     INTENT(OUT) :: ierr

    !> counter on the available linear solvers
    INTEGER isolver

    ! set error code and try solvers until it is cleared
    ierr = 1
    DO isolver = 1, SIZE(config%priority)

      SELECT CASE (config%priority(isolver))

      CASE (bicgstab_multi)

      CASE (bicgstab_no_multi)

      CASE (sgw_linear_solver)

      END SELECT ! solver

      ! if successfully converged, we are done here
      IF (ierr == 0) EXIT

    END DO ! isolver

  END SUBROUTINE select_solver

END MODULE select_solver_module

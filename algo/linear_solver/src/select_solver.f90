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
  INTEGER,  PARAMETER :: bicgstab_multi = 1

  !> use the BiCGstab solver without multishift
  INTEGER,  PARAMETER :: bicgstab_no_multi = 2

  !> use the solver specialized for SternheimerGW
  INTEGER,  PARAMETER :: sgw_linear_solver = 3

  !> configuration of the linear solver
  TYPE select_solver_type

    !> priority in which the linear solvers are chosen
    INTEGER, ALLOCATABLE :: priority(:)

    !> maximum number of iterations allowed
    INTEGER  :: max_iter = 10000

    !> the convergence threshold
    REAL(dp) :: threshold = 1e-4

    !> the number of MR steps for BiCGstab
    INTEGER  :: bicg_lmax = 4

  END TYPE select_solver_type

CONTAINS

  !> select the linear solver based on user preference and convergence
  SUBROUTINE select_solver(config, AA, bb, sigma, xx, ierr)

    USE bicgstab_module,      ONLY: bicgstab, bicgstab_type
    USE io_global,            ONLY: stdout
    USE linear_solver_module, ONLY: linear_solver, linear_solver_config

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

    !> counter on the shifts
    INTEGER ishift

    !> configuration of the BiCGstab solver
    TYPE(bicgstab_type) bicgstab_config

    !> configuration for the SternheimerGW linear solver
    TYPE(linear_solver_config) sgw_solver_config

    !
    ! sanity test of the linear solver
    !
    IF (.NOT.ALLOCATED(config%priority)) THEN
      CALL errore(__FILE__, "priority of the solvers not specified", 1)
    END IF

    ! set error code and try solvers until it is cleared
    ierr = 1
    DO isolver = 1, SIZE(config%priority)

      ! notify user that a different solver will be tried
      IF (isolver > 1) THEN
        WRITE(stdout, '(a)') 'First choice of solver did not converge, try a different one'
      END IF

      !
      ! select the next solver in the priority list and solve the linear problem
      !
      SELECT CASE (config%priority(isolver))

      CASE (bicgstab_multi)
        !
        ! converge using the BiCGstab solver with multishift
        bicgstab_config = bicg_config(config)
        CALL bicgstab(bicgstab_config, AA, bb, sigma, xx, ierr)

      CASE (bicgstab_no_multi)
        !
        ! converge using the BiCGstab solver without multishift
        bicgstab_config = bicg_config(config)
        DO ishift = 1, SIZE(sigma)
          CALL bicgstab(bicgstab_config, AA, bb, sigma, xx, ierr)
        END DO ! ishift

      CASE (sgw_linear_solver)
        !
        ! converge using the SternheimerGW linear solver
        sgw_solver_config = solver_config(config)
        CALL linear_solver(sgw_solver_config, AA, bb, sigma, xx, ierr)

      END SELECT ! solver

      ! if successfully converged, we are done here
      IF (ierr == 0) EXIT

    END DO ! isolver

  END SUBROUTINE select_solver

  !> generate the configuration for the BiCGstab solver
  FUNCTION bicg_config(config)

    USE bicgstab_module, ONLY: bicgstab_type

    !> the general configuration of all linear solvers
    TYPE(select_solver_type), INTENT(IN) :: config

    !> the particular configuration of the BiCGstab solver
    TYPE(bicgstab_type) bicg_config

    !
    ! overwrite the variable in the BiCGstab solver with the ones provided
    ! in the configuration
    !
    bicg_config%max_iter  = config%max_iter
    bicg_config%threshold = config%threshold
    bicg_config%lmax      = config%bicg_lmax

  END FUNCTION bicg_config

  !> generate the configuration for the linear solver
  FUNCTION solver_config(config)

    USE linear_solver_module, ONLY: linear_solver_config

    !> the general configuration of all linear solvers
    TYPE(select_solver_type), INTENT(IN) :: config

    !> the particular configuration of the SternheimerGW linear solver
    TYPE(linear_solver_config) solver_config

    !
    ! overwrite the variable in the linear solver configuration with the
    ! ones provided in the configuration
    !
    solver_config%max_iter  = config%max_iter
    solver_config%threshold = config%threshold

  END FUNCTION solver_config

END MODULE select_solver_module

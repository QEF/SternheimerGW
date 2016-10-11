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
    REAL(dp) :: threshold = 1e-4_dp

    !> maximum number of iterations before aborting
    INTEGER  :: max_iter = 10000

  END TYPE linear_solver_config

  !> store the Krylov subspace that is already generated
  TYPE linear_solver_subspace
  END TYPE linear_solver_subspace

CONTAINS

  !! We use the following algorithm continuously updating the solution and the
  !! residual error of the solution.
  SUBROUTINE linear_solver(config, AA, sigma, BB, XX)

    !> configuration of the linear solver
    TYPE(linear_solver_config), INTENT(IN) :: config

    !> Function pointer that applies the linear operator to a vector.
    INTERFACE
      SUBROUTINE AA(sigma, work)
        USE kinds, ONLY: dp
        !> The shift of this system.
        COMPLEX(dp), INTENT(IN)   :: sigma
        !> *on input* - the vector to which the linear operator is applied <br>
        !! *on output* - the vector after applying the linear operator
        COMPLEX(dp), INTENT(INOUT) :: work(:,:)
      END SUBROUTINE AA
    END INTERFACE

    !> Right hand side of the linear equation.
    COMPLEX(dp), INTENT(IN) :: BB(:,:)

    !> Shift \f$\sigma\f$ in the linear operator.
    COMPLEX(dp), INTENT(IN) :: sigma(:)

    !> On output: the solution of the linear system
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: XX(:,:,:)

    !> the size of the linear problem
    INTEGER vec_size

    !> number of linear problems solved at the same time
    INTEGER num_problem

    !> number of shifts for which we solve
    INTEGER num_shift

    !> counter on the shifts
    INTEGER ishift

    !> counter on the iterations
    INTEGER iter

    !> flag indicating the convergence of the linear solver
    LOGICAL conv

    !> store the Krylov subspace that is already generated
    TYPE(linear_solver_subspace) subspace

    !> work array used for multiple purposes throughout the calculation
    COMPLEX(dp), ALLOCATABLE :: work(:,:)

    ! set helper variables
    vec_size = SIZE(BB, 1)
    num_problem = SIZE(BB, 2)
    num_shift = SIZE(sigma)

    ! create array for solution
    ALLOCATE(XX(vec_size, num_problem, num_shift))

    !! 1. recover previous Krylov subspace information
    CALL linear_solver_recover_subspace(config, vec_size, subspace)

    ! loop over all shifts
    DO ishift = 1, num_shift

      conv = .FALSE.

      ! loop until converged
      DO iter = 1, config%max_iter

        !! 2. determine residual for given shift and right-hand sides 
        ! work contains the residual
        CALL linear_solver_residual(config, sigma(ishift), BB, subspace, work, conv)

        ! exit the loop once convergence is achieved
        IF (conv) EXIT

        !! 3. apply linear operator to residual vector
        ! work contains A applied to the residual
        CALL AA(sigma(ishift), work)

        !! 4. expand the Krylov subspace
        CALL linear_solver_expand_subspace(config, work, subspace)

      END DO ! iter

      ! check for convergence
      IF (.NOT.conv) THEN
        CALL errore(__FILE__, "linear solver did not converge", 1)
      END IF
      
      !! 5. we repeat step 2 - 4 until convergence is achieved, then we determine
      !!    the solution to the linear problem from the subspace information
      CALL linear_solver_obtain_result(config, sigma(ishift), BB, subspace, XX(:, :, ishift))

    END DO ! ishift

  END SUBROUTINE linear_solver

  !> Recover Krylov subspace from memory, file, or generate empty one.
  !!
  !! Depending on the configuration the Krylov subspace is either read from file,
  !! from memory, or a new type is generated.
  SUBROUTINE linear_solver_recover_subspace(config, vec_size, subspace)

    !> configuration of the linear solver
    TYPE(linear_solver_config),   INTENT(IN)  :: config

    !> the dimensionality of the problem
    INTEGER,                      INTENT(IN)  :: vec_size

    !> read the Krylov subspace from previous iteration or generate new one
    TYPE(linear_solver_subspace), INTENT(OUT) :: subspace

  END SUBROUTINE linear_solver_recover_subspace

  !> Determine the residual of the linear equation.
  SUBROUTINE linear_solver_residual(config, sigma, BB, subspace, residual, conv)

    !> configuration of the linear solver
    TYPE(linear_solver_config),   INTENT(IN)  :: config

    !> the shift of the linear problem
    COMPLEX(dp),                  INTENT(IN)  :: sigma

    !> the right hand side of the problem
    COMPLEX(dp),                  INTENT(IN)  :: BB(:,:)

    !> the Krylov subspace generated so far
    TYPE(linear_solver_subspace), INTENT(IN)  :: subspace

    !> the residual error of the solution
    COMPLEX(dp), ALLOCATABLE,     INTENT(OUT) :: residual(:,:)

    !> did the residual error converge
    LOGICAL,                      INTENT(OUT) :: conv

  END SUBROUTINE linear_solver_residual

  !> Expand the subspace to increase the precision of the linear solver
  SUBROUTINE linear_solver_expand_subspace(config, new_vector, subspace)

    !> configuration of the linear solver
    TYPE(linear_solver_config),   INTENT(IN)    :: config

    !> a set of new vectors used to expand the subspace
    COMPLEX(dp),                  INTENT(IN)    :: new_vector(:,:)

    !> the Krylov subspace to be expanded
    TYPE(linear_solver_subspace), INTENT(INOUT) :: subspace

  END SUBROUTINE linear_solver_expand_subspace

  !> Determine the result in
  SUBROUTINE linear_solver_obtain_result(config, sigma, BB, subspace, XX)

    !> configuration of the linear solver
    TYPE(linear_solver_config),   INTENT(IN)  :: config

    !> the shift of the linear problem
    COMPLEX(dp),                  INTENT(IN)  :: sigma

    !> the right hand side of the problem
    COMPLEX(dp),                  INTENT(IN)  :: BB(:,:)

    !> the Krylov subspace generated so far
    TYPE(linear_solver_subspace), INTENT(IN)  :: subspace

    !> the solution of the linear problem
    COMPLEX(dp),                  INTENT(OUT) :: XX(:,:)

  END SUBROUTINE linear_solver_obtain_result

END MODULE linear_solver_module

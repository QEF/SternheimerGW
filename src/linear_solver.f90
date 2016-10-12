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

    !> the relative threshold when the calculation is considered converged
    REAL(dp) :: threshold = 1e-4_dp

    !> maximum number of iterations before aborting
    INTEGER  :: max_iter = 10000

    !> the absolute threshold used to check when the calculation is converged
    REAL(dp), ALLOCATABLE :: abs_threshold(:)

  END TYPE linear_solver_config

  !> Store the Krylov subspace that is already generated.
  !!
  !! This type contains the trial vectors \f$v_i\f$ and the resulting vectors
  !! \f$w_i =  A v_i\f$ after application of the linear operator. We can then
  !! construct a basis from the \f$w_i\f$ and express the right-hand side in
  !! it. Note that the basis has to be orthogonalized again when a shift is
  !! employed.
  TYPE linear_solver_subspace

    !> the last shift that was used
    COMPLEX(dp) sigma

    !> the trial vectors \f$v_i\f$
    COMPLEX(dp), ALLOCATABLE :: vv(:,:)

    !> the basis vectors \f$w_i\f$
    COMPLEX(dp), ALLOCATABLE :: ww(:,:)

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

      !! 2. the subspace is orthogonalized for the current shift
      CALL linear_solver_orthogonal_subspace(config, sigma(ishift), subspace)

      ! loop until converged
      DO iter = 1, config%max_iter

        !! 3. determine residual for given shift and right-hand sides 
        ! work contains the residual
        CALL linear_solver_residual(config, BB, subspace, work, conv)

        ! exit the loop once convergence is achieved
        IF (conv) EXIT

        !! 4. apply linear operator to residual vector
        ! work contains A applied to the residual
        CALL AA(sigma(ishift), work)

        !! 5. expand the Krylov subspace
        CALL linear_solver_expand_subspace(config, work, subspace)

      END DO ! iter

      ! check for convergence
      IF (.NOT.conv) THEN
        CALL errore(__FILE__, "linear solver did not converge", 1)
      END IF
      
      !! 6. we repeat step 3 - 5 until convergence is achieved, then we determine
      !!    the solution to the linear problem from the subspace information
      CALL linear_solver_obtain_result(config, sigma(ishift), BB, subspace, XX(:, :, ishift))

    END DO ! ishift

  END SUBROUTINE linear_solver

  !> Determine the absolute thresholds used from the relative one.
  !!
  !! The absolute threshold is the relative threshold multiplied with the norm
  !! of the right-hand side of the equation.
  SUBROUTINE linear_solver_threshold(BB, rel_threshold, abs_threshold)

    !> the right hand side of the equation
    COMPLEX(dp), INTENT(IN) :: BB(:,:)

    !> the relative threshold to check for convergence
    REAL(dp),    INTENT(IN) :: rel_threshold

    !> the absolute threshold to check for convergence
    REAL(dp),    INTENT(OUT), ALLOCATABLE :: abs_threshold(:)

    !> the dimensionality of the problem
    INTEGER vec_size

    !> the number of right hand side equations
    INTEGER num_problem

    !> counter on the right hand side of the linear problem
    INTEGER iproblem

    !> LAPACK function to evaluate the 2-norm
    REAL(dp), EXTERNAL :: DNRM2

    ! initialize helper variables
    ! note: factor 2 because DNRM2 takes real instead of complex
    vec_size    = 2 * SIZE(BB, 1)
    num_problem = SIZE(BB, 2)

    ALLOCATE(abs_threshold(num_problem))

    DO iproblem = 1, num_problem
      !
      abs_threshold(iproblem) = DNRM2(vec_size, BB(:, iproblem), 1)
      !
    END DO ! iproblem

  END SUBROUTINE linear_solver_threshold

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

  !> orthogonalize the basis of the subspace for the current shift
  !!
  !! The original vectors are constructed with a \f$\sigma_{\text{old}}\f$
  !! \f{equation}{
  !!   (A + \sigma_{\text{old}} I) \vert v_i \rangle = \vert w_i^{\text{old}} \rangle
  !! \f}
  !! and we want to construct a basis for the new shift \f$\sigma_{\text{new}}\f$
  !! \f{equation}{
  !!   (A + \sigma_{\text{new}} I) \vert v_i \rangle = \vert w_i^{\text{new}} \rangle~.
  !! \f}
  !! From the difference of these two equations, we obtain
  !! \f{equation}{
  !!   \vert w_i^{\text{new}} \rangle = \vert w_i^{\text{old}} \rangle +
  !!     (\sigma_{\text{new}} - \sigma_{\text{old}}) \vert v_i \rangle~.
  !! \f}
  !! Then the \f$\vert w_i \rangle\f$ that we want to use as basis functions are
  !! not orthonormal anymore, so that we use a Gram-Schmidt orthonormalization.
  !! In this step, we keep track of the transformation of the \f$\vert v_i\rangle\f$
  !! so that these and the basis functions are still related by the linear operator.
  SUBROUTINE linear_solver_orthogonal_subspace(config, sigma, subspace)

    USE gram_schmidt_module, ONLY: gram_schmidt
    USE kinds,               ONLY: dp

    !> configuration of the linear solver
    TYPE(linear_solver_config),   INTENT(IN)    :: config

    !> the shift of the linear problem
    COMPLEX(dp),                  INTENT(IN)    :: sigma

    !> the Krylov subspace to be orthogonalized
    TYPE(linear_solver_subspace), INTENT(INOUT) :: subspace

    !> the difference of old and new shift
    COMPLEX(dp) diff_sigma

    ! determine the difference of current and last sigma
    diff_sigma = sigma - subspace%sigma

    ! update the w_i = w_i + diff_sigma * v_i
    CALL ZAXPY(SIZE(subspace%ww), diff_sigma, subspace%vv, 1, subspace%ww, 1)

    ! orthogonalize the new w_i keeping the linear relation intact
    CALL gram_schmidt(1, subspace%ww, subspace%vv)

    ! store the current sigma in the subspace
    subspace%sigma = sigma

  END SUBROUTINE linear_solver_orthogonal_subspace

  !> Determine the residual of the linear equation.
  !!
  !! The right hand \f$\vert b \rangle\f$ side of the equation can be expanded
  !! in the basis functions \f$\vert w_i \rangle\f$ of the subspace. The
  !! residual error is given as
  !! \f{equation}{
  !!   \vert r \rangle = \vert b \rangle - \sum_i \vert w_i \rangle \langle w_i \vert b \rangle~.
  !! \f}
  !! If the norm of \f$\vert r\rangle\f$ decreases below a given threshold the
  !! calculation is converged.
  SUBROUTINE linear_solver_residual(config, BB, subspace, residual, conv)

    !> configuration of the linear solver
    TYPE(linear_solver_config),   INTENT(IN)  :: config

    !> the right hand side of the problem
    COMPLEX(dp),                  INTENT(IN)  :: BB(:,:)

    !> the Krylov subspace generated so far
    TYPE(linear_solver_subspace), INTENT(IN)  :: subspace

    !> the residual error of the solution
    COMPLEX(dp), ALLOCATABLE,     INTENT(OUT) :: residual(:,:)

    !> did the residual error converge
    LOGICAL,                      INTENT(OUT) :: conv

    !> the dimensionality of the problem
    INTEGER vec_size

    !> the number of right hand side terms
    INTEGER num_problem

    !> counter on the right hand size terms
    INTEGER iproblem

    !> the number of basis functions in the subspace
    INTEGER num_basis

    !> counter on the basis function
    INTEGER ibasis

    !> the norm of the residual
    REAL(dp) norm

    !> LAPACK function to evaluate the 2-norm
    REAL(dp), EXTERNAL :: DNRM2

    !> overlap of right-hand side and basis vector
    COMPLEX(dp) overlap

    !> LAPACK routine to evaluate dot product
    COMPLEX(dp), EXTERNAL :: ZDOTC

    ! set helper variables
    vec_size    = SIZE(BB, 1)
    num_problem = SIZE(BB, 2)
    num_basis   = SIZE(subspace%ww, 2)

    ! convergence flag is initially set, but cleared if any problem did not converge
    conv = .TRUE.

    ! copy B to residual
    ALLOCATE(residual(vec_size, num_problem))
    CALL ZCOPY(SIZE(BB), BB, 1, residual, 1)

    ! loop over right-hand side
    DO iproblem = 1, num_problem

      ! loop over basis functions
      DO ibasis = 1, num_basis

        ! r = b - (w_i, b) w_i
        overlap = ZDOTC(vec_size, BB(:,iproblem), 1, subspace%ww(:,ibasis), 1)
        residual(:,iproblem) = residual(:,iproblem) - overlap * subspace%ww(:,ibasis)

      END DO ! ibasis

      ! check convergence (factor 2 because of complex)
      IF (conv) THEN
        norm = DNRM2(2 * vec_size, residual(:,iproblem), 1)
        conv = conv .AND. (norm < config%abs_threshold(iproblem))
      END IF

    END DO ! iproblem

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

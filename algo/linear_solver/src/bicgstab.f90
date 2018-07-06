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
!> Shifted BiCGStab(l) solver for linear equation \f$(A + \sigma I) x = b\f$.
!!
!! Implements the shifted BiCGStab(l) algorithm according to the paper of
!! Frommer, Computing **70**, 87 (2003). The general idea is that the matrix
!! \f$A\f$ and \f$A + \sigma I\f$ span the same Krylov subspace. Hence, we can
!! solve the linear equation of all linear problems at the cost of a single one.
!!
!! Notice that the variables \f$\theta\f$, \f$\phi\f$, and \f$\xi\f$ are stored
!! as the inverse of the ones defined in the algorithm to improve the numeric
!! stability.
MODULE bicgstab_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC bicgstab, bicgstab_type

  !> Configuration of the BiCGstab solver
  TYPE bicgstab_type

    !> the number of MR steps before resetting
    INTEGER  :: lmax = 4

    !> the convergence threshold
    REAL(dp) :: threshold = 1e-4_dp

    !> the maximum number of iterations
    INTEGER  :: max_iter = 10000

  END TYPE bicgstab_type

  !> Variables and arrays used for the seed system.
  TYPE seed_system_type

    !> search direction
    COMPLEX(dp), ALLOCATABLE :: uu(:,:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

    !> residual vector
    COMPLEX(dp), ALLOCATABLE :: rr(:,:)

    !> \f$\tilde r_0\f$ of Frommer's algorithm.
    COMPLEX(dp), ALLOCATABLE :: tilde_r0(:)

    !> the shift of this system
    COMPLEX(dp) sigma

    !> rho contains dot product of residuals
    COMPLEX(dp) rho

    !> \f$\rho_{\text{old}}\f$ of Frommer's algorithm
    COMPLEX(dp) rho_old

    !> \f$\alpha\f$ of Frommer's algorithm
    COMPLEX(dp) alpha

    !> \f$\alpha_{\text{old}}\f$ of Frommer's algorithm
    COMPLEX(dp) alpha_old

    !> beta determines step size
    COMPLEX(dp) beta

    !> \f$\omega\f$ of Frommer's algorithm
    COMPLEX(dp) omega
  
    !> \f$\gamma\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: gamma(:)

    !> \f$\gamma'\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: gamma_p(:)

    !> \f$\gamma''\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: gamma_pp(:)

  END TYPE seed_system_type

  !> Variables and arrays used for the shifted system.
  TYPE shift_system_type

    !> search direction
    COMPLEX(dp), ALLOCATABLE :: uu(:,:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

    !> the shift of this system
    COMPLEX(dp) sigma

    !> The inverse of \f$\phi_{\text{old}}^\sigma\f$ of Frommer's algorithm
    COMPLEX(dp) inv_phi_old

    !> The inverse of \f$\phi^\sigma\f$ of Frommer's algorithm
    COMPLEX(dp) inv_phi

    !> The inverse of \f$\phi_{\text{new}}^\sigma\f$ of Frommer's algorithm
    COMPLEX(dp) inv_phi_new

    !> The inverse of \f$\vartheta^\sigma\f$ of Frommer's algorithm
    COMPLEX(dp) inv_theta

    !> \f$\alpha^\sigma\f$ of Frommer's algorithm
    COMPLEX(dp) alpha

    !> \f$\beta^\sigma\f$ of Frommer's algorithm
    COMPLEX(dp) beta

    !> \f$\mu_{ij}\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: mu(:)

    !> \f$\gamma\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: gamma(:)

    !> \f$\gamma'\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: gamma_p(:)

    !> \f$\gamma''\f$ of Frommer's algorithm
    COMPLEX(dp), ALLOCATABLE :: gamma_pp(:)

  END TYPE shift_system_type

CONTAINS

  !> Main driver routine of the algorithm.
  !!
  !! This subroutine implements the *Algorithm 2* of Frommer's paper.
  !!
  SUBROUTINE bicgstab(config, AA, bb, sigma, xx, ierr)

    USE debug_module, ONLY: test_nan
    USE io_global,    ONLY: stdout

    USE norm_module, ONLY: norm
    !> the configuration of the BiCGstab solver
    TYPE(bicgstab_type), INTENT(IN) :: config

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
    COMPLEX(dp), INTENT(IN)  :: bb(:)

    !> Shift \f$\sigma\f$ in the linear operator. The first element of the array
    !! will be used as seed system. The other ones as shifted systems.
    COMPLEX(dp), INTENT(IN)  :: sigma(:)

    !> On output: the solution of the linear system
    COMPLEX(dp), INTENT(OUT) :: xx(:,:)

    !> error code returned if the solver did not converge
    INTEGER,     INTENT(OUT) :: ierr

    !> Size of the right-hand vector.
    INTEGER vec_size

    !> loop over shifted systems
    INTEGER ishift

    !> Counter for number of iterations.
    INTEGER iter

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the seed system.
    TYPE(seed_system_type)  seed_system

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the shifted system.
    TYPE(shift_system_type) shift_system(SIZE(sigma) - 1)

    ! determine dimension of vector
    vec_size = SIZE(bb)

    !
    ! sanity check of the input
    !
    IF (SIZE(xx, 1) /= vec_size) &
      CALL errore(__FILE__, "right-hand side and solution must have same size", 1)
    IF (SIZE(xx, 2) /= SIZE(sigma)) &
      CALL errore(__FILE__, "we need one solution vector per shift for the result", 1)
    ierr = 0

    !
    ! initialization seed system
    !
    CALL init_seed(config%lmax, bb, sigma(1), seed_system)

    !
    ! initialization shifted systems
    !
    CALL init_shift(config%lmax, vec_size, sigma, shift_system)

    ! loop until solution is found
    DO iter = 1, config%max_iter

      !
      ! perform BiCG part (Algorithm 3)
      !
      CALL bicg_part(config%lmax, AA, seed_system, shift_system)

      ! stop loop if result is converged
      IF (converged(config%threshold, seed_system%rr(:,0))) EXIT

      !
      ! perform MR part (Algorithm 4)
      !
      CALL mr_part(config%lmax, seed_system, shift_system)

      ! stop loop if result is converged
      IF (converged(config%threshold, seed_system%rr(:,0))) EXIT

    END DO ! iter

    IF (iter > config%max_iter) THEN
      WRITE(stdout, '(a,1x,i0,1x,a)') "WARNING: BiCGstab algorithm did not converge in", &
                                      config%max_iter, "iterations."
      ierr = 1
    END IF

    !
    ! copy the result to the output array
    !
    CALL ZCOPY(vec_size, seed_system%xx, 1, xx(:,1), 1)
    DO ishift = 2, SIZE(sigma)
      CALL ZCOPY(vec_size, shift_system(ishift - 1)%xx, 1, xx(:,ishift), 1)
    END DO ! ishift

    ! check if output is sane
    IF (ANY(test_nan(xx))) THEN
      WRITE(stdout, *) 'WARNING: BiCGstab algorithm resulted in one or more NaN'
      ierr = 2
    END IF

    !
    ! destroy the allocated array in the types
    !
    CALL destroy_seed(seed_system)
    CALL destroy_shift(shift_system)

  END SUBROUTINE bicgstab

  !> Check if the residual is below the given threshold.
  FUNCTION converged(threshold, residual)

    USE norm_module, ONLY: norm

    !> If the norm of the residual drops below this threshold, the system is
    !! considered as converged.
    REAL(dp),    INTENT(IN) :: threshold

    !> The residual in the seed system of the multishift problem.
    COMPLEX(dp), INTENT(IN) :: residual(:)

    !> returns true, if the norm of the residual is smaller than the threshold
    LOGICAL converged

    !> norm of the residual
    REAL(dp) norm_residual

    ! check residual of seed system
    norm_residual = norm(residual)

    ! if the norm is smaller than the threshold the system is converged
    converged = (norm_residual < threshold)

  END FUNCTION converged

  !> Initialize the seed system.
  SUBROUTINE init_seed(lmax, bb, sigma, seed_system)

    !> Dimensionality of the GMRES algorithm.
    INTEGER,     INTENT(IN) :: lmax

    !> Initial residual (right hand side of equation).
    COMPLEX(dp), INTENT(IN) :: bb(:)

    !> The shift of the system that we will use as the seed system.
    COMPLEX(dp), INTENT(IN) :: sigma

    !> On output contains the initialized seed system.
    TYPE(seed_system_type), INTENT(OUT) :: seed_system

    !> Size of the vectors to initialize.
    INTEGER vec_size

    ! determine size of vector
    vec_size = SIZE(bb)

    ! allocate arrays of given size
    ALLOCATE(seed_system%uu(vec_size, 0:lmax))
    ALLOCATE(seed_system%xx(vec_size))
    ALLOCATE(seed_system%rr(vec_size, 0:lmax))
    ALLOCATE(seed_system%tilde_r0(vec_size))
    ALLOCATE(seed_system%gamma(lmax))
    ALLOCATE(seed_system%gamma_p(lmax))
    ALLOCATE(seed_system%gamma_pp(lmax))

    ! initialize arrays
    ! initialize only u0, because the other ones will be set in the algorithm 
    seed_system%uu(:,0) = 0
    seed_system%xx = 0

    ! initial residual is r0 = (b - A x) = b, because initial x = 0
    ! initialize only r0, because the other ones will be set in the algorithm 
    seed_system%rr(:,0) = bb
    seed_system%tilde_r0 = bb

    ! init the variables
    seed_system%sigma     = sigma
    seed_system%rho_old   = 1.0
    seed_system%alpha_old = 1.0
    seed_system%alpha     = 0.0
    seed_system%omega     = 1.0

  END SUBROUTINE init_seed

  !> Deallocate the arrays in the seed system type.
  SUBROUTINE destroy_seed(seed_system)

    !> Deallocate the arrays of this system.
    TYPE(seed_system_type), INTENT(INOUT) :: seed_system

    ! deallocate the arrays
    DEALLOCATE(seed_system%uu)
    DEALLOCATE(seed_system%xx)
    DEALLOCATE(seed_system%rr) 
    DEALLOCATE(seed_system%tilde_r0) 
    DEALLOCATE(seed_system%gamma)
    DEALLOCATE(seed_system%gamma_p)
    DEALLOCATE(seed_system%gamma_pp)

  END SUBROUTINE destroy_seed

  !> Initialize the shifted systems.
  SUBROUTINE init_shift(lmax, vec_size, sigma, shift_system)

    !> Dimensionality of the GMRES algorithm.
    INTEGER,     INTENT(IN) :: lmax

    !> Size of the vectors to initialize.
    INTEGER, INTENT(IN) :: vec_size

    !> Shifts \f$\sigma\f$ relative to the seed system.
    COMPLEX(dp), INTENT(IN) :: sigma(:)

    !> On output contains the initialized shifted systems.
    TYPE(shift_system_type), INTENT(OUT) :: shift_system(SIZE(sigma) - 1)

    !> counter over shifted systems
    INTEGER ishift

    !> dimension of the array mu
    INTEGER mu_size

    !> binomials \f$\begin{pmatrix}j \\ i\end{pmatrix}\f$
    REAL(dp),    ALLOCATABLE :: binomial(:)

    !> array to store \f$\sigma^n\f$ for \f$n = 0, \ldots, l_{\text{max}}\f$
    COMPLEX(dp), ALLOCATABLE :: sigma_pow(:)

    !> counter variables to evaluate mu
    INTEGER ii, jj, offset, ij

    ! determine size of array mu
    mu_size = lmax * (lmax + 1) / 2

    ! allocate helper arrays
    ALLOCATE(binomial(mu_size))
    ALLOCATE(sigma_pow(0:lmax - 1))

    !                     / j \        j!
    ! evaluate binomials (     ) = -----------
    !                     \ i /    i! (j - i)!
    DO jj = 0, lmax - 1

      ! the ij index starts at this value
      offset = jj * (jj + 1) / 2 + 1

      DO ii = 0, jj

        ! evaluate index in array
        ij = offset + ii

        ! trivial case i = 0: binomial = 1
        IF (ii == 0) THEN
          binomial(ij) = 1

        ! in general, the following is valid
        !  / j \     /   j   \  j - i + 1
        ! (     ) = (         ) ---------
        !  \ i /     \ i - 1 /      i
        ELSE
          binomial(ij) = (binomial(ij - 1) * (jj - ii + 1)) / ii

        END IF

      END DO ! i
    END DO ! j

    ! loop over all shifted systems
    DO ishift = 1, SIZE(sigma) - 1

      ! allocate arrays of given size
      ALLOCATE(shift_system(ishift)%uu(vec_size, 0:lmax))
      ALLOCATE(shift_system(ishift)%xx(vec_size))
      ALLOCATE(shift_system(ishift)%gamma(lmax))
      ALLOCATE(shift_system(ishift)%gamma_p(lmax))
      ALLOCATE(shift_system(ishift)%gamma_pp(lmax))

      ! allocate array for mu
      ALLOCATE(shift_system(ishift)%mu(mu_size))

      ! initialize arrays to 0
      ! initialize only u0, because the other ones will be set in the algorithm 
      shift_system(ishift)%uu(:,0) = 0
      shift_system(ishift)%xx = 0

      ! initialize the variables
      shift_system(ishift)%inv_phi_old = 1.0
      shift_system(ishift)%inv_phi     = 1.0
      shift_system(ishift)%inv_theta   = 1.0
      ! subtract the shift of the initial system
      shift_system(ishift)%sigma       = sigma(ishift + 1) - sigma(1)

      ! construct sigma_pow array
      sigma_pow(0) = 1.0
      DO ii = 1, lmax - 1
        sigma_pow(ii) = shift_system(ishift)%sigma * sigma_pow(ii - 1)
      END DO ! i

      ! construct mu_ij
      DO jj = 0, lmax - 1

        ! the ij index starts at this value
        offset = jj * (jj + 1) / 2 + 1

        DO ii = 0, jj

          ! evaluate index in array
          ij = offset + ii

          !          / j \       j - i
          ! mu_ij = (     ) sigma
          !          \ i /
          shift_system(ishift)%mu(ij) = binomial(ij) * sigma_pow(jj - ii)

        END DO ! i
      END DO ! j

    END DO ! ishift

    DEALLOCATE(binomial)
    DEALLOCATE(sigma_pow)

  END SUBROUTINE init_shift

  !> Deallocate the arrays in the shifted system types.
  SUBROUTINE destroy_shift(shift_system)

    !> Deallocate the arrays in all of this types.
    TYPE(shift_system_type), INTENT(INOUT) :: shift_system(:)

    !> counter over shifted systems
    INTEGER ishift

    ! loop over shifted systems
    DO ishift = 1, SIZE(shift_system)

      ! deallocate arrays
      DEALLOCATE(shift_system(ishift)%uu)
      DEALLOCATE(shift_system(ishift)%xx)
      DEALLOCATE(shift_system(ishift)%mu)
      DEALLOCATE(shift_system(ishift)%gamma)
      DEALLOCATE(shift_system(ishift)%gamma_p)
      DEALLOCATE(shift_system(ishift)%gamma_pp)

    END DO ! ishift

  END SUBROUTINE destroy_shift

  !> The BiCG part of the BiCGstab algorithm.
  !!
  !! This subroutine implements the *Algorithm 3* of Frommer's paper.
  !!
  SUBROUTINE bicg_part(lmax, AA, seed_system, shift_system)

    !> Dimensionality of the GMRES algorithm.
    INTEGER, INTENT(IN) :: lmax

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

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the seed system. Updated
    !! to the next best guess at the end of the routine.
    TYPE(seed_system_type),  INTENT(INOUT) :: seed_system

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the shifted system. Updated
    !! to the next best guess at the end of the routine.
    TYPE(shift_system_type), INTENT(INOUT), TARGET :: shift_system(:)

    !> pointer to the currently active shifted system
    TYPE(shift_system_type), POINTER :: active

    !> Size of the vectors to initialize.
    INTEGER vec_size

    !> Number of shifted systems
    INTEGER num_shift

    !> Counter for the shifted systems
    INTEGER ishift

    !> counter for previous iterations
    INTEGER ii

    !> counter for GMRES iteration
    INTEGER jj

    !> helper variable to store factor
    COMPLEX(dp) factor

    !> complex 1
    COMPLEX(dp), PARAMETER :: one = CMPLX(1, 0)

    ! BLAS dot product routine
    COMPLEX(dp), EXTERNAL :: ZDOTU

    ! determine vector size
    vec_size = SIZE(seed_system%tilde_r0)

    ! determine number of shifted systems
    num_shift = SIZE(shift_system)

    !
    ! Frommer's algorithm 3
    ! Note: the L?? refers to the line numbers given in Frommer's paper
    !       we abbreviate the superscript sigma with ^
    !
    ! L3: rho_old = - omega * rho_old
    seed_system%rho_old = - seed_system%omega * seed_system%rho_old
    !
    ! L4: loop over GMRES iterations
    DO jj = 0, lmax - 1
      !
      ! seed system
      !
      ! L6: rho = (r_j, ~r_0),
      seed_system%rho = ZDOTU(vec_size, seed_system%rr(:,jj), 1, seed_system%tilde_r0, 1)
      !     beta = alpha * rho / rho_old,
      seed_system%beta = seed_system%alpha * seed_system%rho / seed_system%rho_old
      !     rho_old = rho
      seed_system%rho_old = seed_system%rho
      !
      ! L7: loop over previous steps
      DO ii = 0, jj
        !
        ! L8: u_i = r_i - beta u_i
        CALL ZSCAL(vec_size, -seed_system%beta, seed_system%uu(:,ii), 1)
        CALL ZAXPY(vec_size, one, seed_system%rr(:,ii), 1, seed_system%uu(:,ii), 1)
        !
      END DO ! i
      !
      ! L10: u_j+1 = A u_j
      CALL AA(seed_system%sigma, seed_system%uu(:,jj), seed_system%uu(:, jj + 1))
      !
      ! L11: alpha = rho / (u_j+1, ~r_0)
      seed_system%alpha = seed_system%rho &
                        / ZDOTU(vec_size, seed_system%uu(:, jj + 1), 1, &
                                seed_system%tilde_r0, 1)
      !
      ! shifted system
      !
      ! loop over all shifted systems
      DO ishift = 1, num_shift
        !
        active => shift_system(ishift)
        !
        ! L13: phi_new^ = (1 + alpha sigma) phi^ + alpha beta / alpha_old (phi_old^ - phi^)
        ! <=>                                        (1/phi^)
        ! (1/phi_new^) = ------------------------------------------------------------------------
        !                (1 + alpha sigma) + alpha beta / alpha_old [(1/phi^) / (1/phi_old^) - 1]
        active%inv_phi_new = active%inv_phi / (one + seed_system%alpha * active%sigma &
                           + seed_system%alpha * seed_system%beta / seed_system%alpha_old &
                           * (active%inv_phi / active%inv_phi_old - one))
        !      beta^ = (phi_old^ / phi^)**2 beta
        active%beta = (active%inv_phi / active%inv_phi_old)**2 * seed_system%beta
        !      alpha^ = (phi^ / phi_new^) alpha
        active%alpha = (active%inv_phi_new / active%inv_phi) * seed_system%alpha
        !
        ! evaluate 1 / (theta^ phi^) used in L15 and L17
        factor = active%inv_theta * active%inv_phi
        !
        ! L14: loop over previous steps
        DO ii = 0, jj
          !
          ! L15: u_i^ = r_i / (theta^ phi^) - beta^ u_i^
          CALL ZSCAL(vec_size, -active%beta, active%uu(:,ii), 1)
          CALL ZAXPY(vec_size, factor, seed_system%rr(:,ii), 1, active%uu(:,ii), 1)
          !
        END DO ! i
        !
        ! L17: x^ = x^ + alpha^ u0^
        CALL ZAXPY(vec_size, active%alpha, active%uu(:,0), 1, active%xx, 1)
        !
        ! note - we shifted update of alpha after the loop to allow
        !        for systems with multiple shifts
        ! L18: phi_old^ = phi^
        active%inv_phi_old = active%inv_phi
        !      phi^ = phi_new^
        active%inv_phi = active%inv_phi_new
        !
        ! note: we can use factor here, because phi_old = phi in L18
        ! L19: u_j+1^ = r_j / (theta^ phi_old^)
        CALL ZCOPY(vec_size, seed_system%rr(:,jj), 1, active%uu(:, jj + 1), 1)
        CALL ZSCAL(vec_size, factor, active%uu(:, jj + 1), 1)
        !
      END DO ! ishift
      !
      ! this update must be done here instead of with the rest of
      ! line 18, because it acts on the seed_system and would affect
      ! the other shifted systems otherwise
      ! L18: alpha_old = alpha
      seed_system%alpha_old = seed_system%alpha
      !
      ! L21: loop over previous iterations
      DO ii = 0, jj
        !
        ! L22: r_i = r_i - alpha u_i+1
        CALL ZAXPY(vec_size, -seed_system%alpha, seed_system%uu(:, ii + 1), 1, &
                   seed_system%rr(:,ii), 1)
        !
      END DO ! i
      !
      ! L24: r_j+1 = A r_j
      CALL AA(seed_system%sigma, seed_system%rr(:,jj), seed_system%rr(:, jj + 1))
      !
      ! L25: x = x + alpha u0
      CALL ZAXPY(vec_size, seed_system%alpha, seed_system%uu(:,0), 1, seed_system%xx, 1)
      !
      ! now update shifted systems again
      DO ishift = 1, num_shift
        !
        active => shift_system(ishift)
        !
        ! L26: u_j+1^ = (u_j+1^ - r_j / (theta^ phi^)) / alpha^ - sigma u_j^
        factor = -active%inv_theta * active%inv_phi
        CALL ZAXPY(vec_size, factor, seed_system%rr(:,jj), 1, active%uu(:, jj + 1), 1)
        CALL ZSCAL(vec_size, 1.0 / active%alpha, active%uu(:, jj + 1), 1)
        CALL ZAXPY(vec_size, -active%sigma, active%uu(:,jj), 1, active%uu(:, jj + 1), 1)
        !
      END DO ! ishift
      !
    END DO ! j

  END SUBROUTINE bicg_part

  !> The MR part of the BiCGstab algorithm.
  !!
  !! This subroutine implements the *Algorithm 4* of Frommer's paper.
  !!
  SUBROUTINE mr_part(lmax, seed_system, shift_system)

    !> Dimensionality of the GMRES algorithm.
    INTEGER, INTENT(IN) :: lmax

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the seed system. Updated
    !! to the next best guess at the end of the routine.
    TYPE(seed_system_type),  INTENT(INOUT) :: seed_system

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the shifted system. Updated
    !! to the next best guess at the end of the routine.
    TYPE(shift_system_type), INTENT(INOUT), TARGET :: shift_system(:)

    !> pointer to the currently active shifted system
    TYPE(shift_system_type), POINTER :: active

    !> Size of the vectors to initialize.
    INTEGER vec_size

    !> Number of shifted systems
    INTEGER num_shift

    !> Counter for the shifted systems
    INTEGER ishift

    !> counter for previous iterations
    INTEGER ii

    !> counter for GMRES iteration
    INTEGER jj

    !> helper to access the right element
    INTEGER offset, ij

    !> Temporary storage for prefactors
    COMPLEX(dp) factor

    !> The inverse of the variable \f$\xi\f$ of Frommer's algorithm.
    COMPLEX(dp) inv_xi

    !> The variable \f$\psi\f$ of Frommer's algorithm.
    COMPLEX(dp) psi

    !> The vector \f$\nu\f$ of Frommer's algorithm.
    COMPLEX(dp), ALLOCATABLE :: nu(:)

    !> The array \f$\tau\f$ of Frommer's algorithm.
    COMPLEX(dp), ALLOCATABLE :: tau(:)

    ! BLAS dot product routine
    COMPLEX(dp), EXTERNAL :: ZDOTU

    ! allocate vectors of appropriate size
    ALLOCATE(nu(lmax))
    ALLOCATE(tau((lmax + 1) * (lmax + 2) / 2))

    ! determine vector size
    vec_size = SIZE(seed_system%tilde_r0)

    ! determine number of shifted systems
    num_shift = SIZE(shift_system)

    !
    ! Frommer's algorithm 4
    ! Note: the L?? refers to the line numbers given in Frommer's paper
    !       we abbreviate the superscript sigma with ^
    !
    ! modified Gram-Schmidt orthogonalization
    !
    ! L1: loop over GMRES iterations
    DO jj = 1, lmax
      !
      offset = jj * (jj + 1) / 2 + 1
      !
      ! L2: loop over previous iterations
      DO ii = 1, jj - 1
        !
        ij = offset + ii
        !
        ! L3: tau_ij = (r_j, r_i) / nu_i
        tau(ij) = ZDOTU(vec_size, seed_system%rr(:,jj), 1, seed_system%rr(:,ii), 1) / nu(ii)
        !     r_j = r_j - tau_ij r_i
        CALL ZAXPY(vec_size, -tau(ij), seed_system%rr(:,ii), 1, seed_system%rr(:,jj), 1)
        !
      END DO ! i
      !
      ! L5: nu_j = (r_j, r_j)
      nu(jj) = ZDOTU(vec_size, seed_system%rr(:,jj), 1, seed_system%rr(:,jj), 1)
      !     gamma_j' = (r_0, r_j) / nu(j)
      seed_system%gamma_p(jj) = ZDOTU(vec_size, seed_system%rr(:,0), 1, seed_system%rr(:,jj), 1) / nu(jj)
      !
    END DO ! j
    !
    ! L7: gamma_l = gamma_l'
    seed_system%gamma(lmax) = seed_system%gamma_p(lmax)
    !     omega = gamma_l
    seed_system%omega = seed_system%gamma(lmax)
    !
    ! L8: sum over tau gamma
    DO jj = lmax - 1, 1, -1
      !
      ! L9: gamma_j = gamma_j' - sum_{i = j + 1}^l tau_ji gamma_i
      seed_system%gamma(jj) = seed_system%gamma_p(jj)
      DO ii = jj + 1, lmax
        ij = ii * (ii + 1) / 2 + jj + 1
        seed_system%gamma(jj) = seed_system%gamma(jj) - tau(ij) * seed_system%gamma(ii)
      END DO ! i
      !
    END DO ! j
    !
    ! L11: sum over tau gamma
    DO jj = 1, lmax - 1
      !
      ! L12: gamma_j" = gamma_j+1 + sum_{i = j + 1}^{l-1} tau_ji gamma_i+1
      seed_system%gamma_pp(jj) = seed_system%gamma(jj + 1)
      DO ii = jj + 1, lmax - 1
        ij = ii * (ii + 1) / 2 + jj + 1
        seed_system%gamma_pp(jj) = seed_system%gamma_pp(jj) + tau(ij) * seed_system%gamma(ii + 1)
      END DO ! i
    END DO ! j
    !
    ! update step
    ! L14: x = x + gamma_1 r_0
    CALL ZAXPY(vec_size, seed_system%gamma(1), seed_system%rr(:,0), 1, seed_system%xx, 1)
    !      u_0 = u_0 - gamma_l u_l
    CALL ZAXPY(vec_size, -seed_system%gamma(lmax), seed_system%uu(:,lmax), 1, seed_system%uu(:,0), 1)
    !
    ! L15: loop over GMRES iterations
    DO jj = 1, lmax - 1
      !
      ! L16: u_0 = u_0 - gamma_j u_j
      CALL ZAXPY(vec_size, -seed_system%gamma(jj), seed_system%uu(:,jj), 1, seed_system%uu(:,0), 1)
      !
      ! L17: x = x + gamma_j" r_j
      CALL ZAXPY(vec_size, seed_system%gamma_pp(jj), seed_system%rr(:,jj), 1, seed_system%xx, 1)
      !
    END DO ! j
    !
    ! shifted systems
    !
    DO ishift = 1, num_shift
      !
      active => shift_system(ishift)
      !
      ! L20: call Horner's scheme to obtain psi^ and gamma_j^
      CALL horner_scheme(seed_system%gamma, active%sigma, active%gamma, psi)
      !
      ! L21: xi = theta^ phi^
      inv_xi = active%inv_theta * active%inv_phi
      !      theta^ = theta^ psi
      active%inv_theta = active%inv_theta / psi 
      !
      ! L22: loop over GMRES iterations
      DO jj = 1, lmax
        !
        ! L23: gamma_j'^ = sum_{i = j}^l mu_{j - 1, i - 1} gamma_i^
        active%gamma_p(jj) = 0
        DO ii = jj, lmax
          ij = (ii - 1) * ii / 2 + jj
          active%gamma_p(jj) = active%gamma_p(jj) + active%mu(ij) * active%gamma(ii)
        END DO ! i
        !
      END DO ! j
      !
      ! L25: loop over GMRES iterations
      DO jj = 1, lmax - 1
        !
        ! L26: gamma_j"^ = gamma_j+1'^ + sum_{i = j + 1}^{l - 1} tau_ji gamma_i+1'^
        active%gamma_pp(jj) = active%gamma_p(jj + 1)
        DO ii = jj + 1, lmax - 1
          ij = ii * (ii + 1) / 2 + jj + 1
          active%gamma_pp(jj) = active%gamma_pp(jj) + tau(ij) * active%gamma_p(ii + 1)
        END DO ! i
        !
      END DO ! j
      !
      ! note: the residuum update is done later, outside of the loop
      ! L28: x^ = x^ + gamma_1'^ / xi r_0
      factor = active%gamma_p(1) * inv_xi
      CALL ZAXPY(vec_size, factor, seed_system%rr(:,0), 1, active%xx, 1)
      !      u_0^ = u_0^ - gamma_l u_l^
      ! note: the latter equation is incorrect in Frommer's paper
      !       where u_l instead of u_l^ is used
      CALL ZAXPY(vec_size, -seed_system%gamma(lmax), active%uu(:,lmax), 1, &
                 active%uu(:,0), 1)
      !
      ! L29: loop over GMRES iterations
      DO jj = 1, lmax - 1
        !
        ! L30: u_0^ = u_0^ - gamma_j u_j^
        CALL ZAXPY(vec_size, -seed_system%gamma(jj), active%uu(:,jj), 1, &
                   active%uu(:,0), 1)
        !
        ! L31: x^ = x^ + (gamma_j"^ / xi) r_j
        factor = active%gamma_pp(jj) * inv_xi
        CALL ZAXPY(vec_size, factor, seed_system%rr(:,jj), 1, active%xx, 1)
        !
      END DO ! j
      !
      ! L34: u_0^ = u_0^ / psi
      factor = 1 / psi
      CALL ZSCAL(vec_size, factor, active%uu(:,0), 1)
      !
    END DO ! ishift
    !
    ! now do the delayed residuum update
    ! note that the effect of L28 is to shift the upper boundary of the loop
    !
    ! L29: loop over GMRES iterations
    DO jj = 1, lmax
      !
      ! L32: r_0 = r_0 - gamma_j' r_j
      CALL ZAXPY(vec_size, -seed_system%gamma_p(jj), seed_system%rr(:,jj), 1, &
                 seed_system%rr(:,0), 1)
      !
    END DO ! j

    ! free memory
    DEALLOCATE(tau)
    DEALLOCATE(nu)

  END SUBROUTINE mr_part

  !> The Horner's scheme part of the BiCGstab algorithm.
  !!
  !! This subroutine implements the *Algorithm 5* of Frommer's paper.
  !!
  !! Computes \f$\gamma^\sigma_j$ and $\psi^\sigma$ for the shifted
  !! system from the values in the seed system and the shift 
  !! \f$\sigma\f$.
  !!
  SUBROUTINE horner_scheme(seed_gamma, sigma, shift_gamma, psi)

    !> The \f$gamma\f$ values in the seed system.
    COMPLEX(dp), INTENT(IN)  :: seed_gamma(:)

    !> The shift relative to the seed system.
    COMPLEX(dp), INTENT(IN)  :: sigma

    !> The \f$\gamma^\sigma\f$ values in the shifted system.
    COMPLEX(dp), INTENT(OUT) :: shift_gamma(:)

    !> The \f$\psi^\sigma\f$ value for the shifted system.
    COMPLEX(dp), INTENT(OUT) :: psi

    !> dimension of the gamma arrays
    INTEGER lmax

    !> loop variables
    INTEGER ii, jj

    !> determine size of the array
    lmax = SIZE(seed_gamma)

    !
    ! Frommer's algorithm 5
    ! Note: we abbreviate the superscript sigma with ^
    !
    ! gamma_l^ = -gamma_l
    shift_gamma(lmax) = -seed_gamma(lmax)
    !
    ! loop over GMRES iterations
    DO jj = lmax - 1, 1, -1
      !
      ! gamma_j^ = -sigma gamma_j+1^ - gamma_j
      shift_gamma(jj) = -sigma * shift_gamma(jj + 1) - seed_gamma(jj)
      !
    END DO ! j
    !
    ! psi^ = -sigma gamma_1^ + 1
    psi = -sigma * shift_gamma(1) + 1
    !
    DO ii = 1, lmax - 1
      DO jj = lmax - 1, ii, -1
        !
        ! gamma_j^ = -sigma gamma_j+1^ + gamma_j^
        shift_gamma(jj) = -sigma * shift_gamma(jj + 1) + shift_gamma(jj)
        !
      END DO ! j
    END DO ! i
    !
    ! loop over GMRES iterations
    DO jj = 1, lmax
      !
      ! gamma_j^ = -gamma_j^ / psi^
      shift_gamma(jj) = -shift_gamma(jj) / psi
      !
    END DO ! j

  END SUBROUTINE horner_scheme

END MODULE bicgstab_module

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
!> Shifted BiCGStab(l) solver for linear equation \f$(A + \sigma I) x = b\f$.
!!
!! Implements the shifted BiCGStab(l) algorithm according to the paper of
!! Fromme, Computing **70**, 87 (2003). The general idea is that the matrix
!! \f$A\f$ and \f$A + \sigma I\f$ span the same Krylov subspace. Hence, we can
!! solve the linear equation of all linear problems at the cost of a single one.
!!
MODULE bicgstab_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC bicgstab

  !> Variables and arrays used for the seed system.
  TYPE seed_system_type

    !> search direction
    COMPLEX(dp), ALLOCATABLE :: u0(:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

    !> residual vector
    COMPLEX(dp), ALLOCATABLE :: r0(:)

    !> \f$\tilde r_0\f$ of Fromme's algorithm.
    COMPLEX(dp), ALLOCATABLE :: tilde_r0(:)

    !> rho contains dot product of residuals
    COMPLEX(dp) rho

    !> \f$\rho_{\text{old}}\f$ of Fromme's algorithm
    COMPLEX(dp) rho_old

    !> \f$\alpha\f$ of Fromme's algorithm
    COMPLEX(dp) alpha

    !> \f$\alpha_{\text{old}}\f$ of Fromme's algorithm
    COMPLEX(dp) alpha_old

    !> beta determines step size
    COMPLEX(dp) beta

    !> \f$\omega\f$ of Fromme's algorithm
    COMPLEX(dp) omega
  
  END TYPE seed_system_type

  !> Variables and arrays used for the shifted system.
  TYPE shift_system_type

    !> search direction
    COMPLEX(dp), ALLOCATABLE :: u0(:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

    !> the shift of this system
    COMPLEX(dp) sigma

    !> \f$\phi_{\text{old}}^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) phi_old

    !> \f$\phi^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) phi

    !> \f$\phi_{\text{new}}^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) phi_new

    !> \f$\vartheta^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) theta

    !> \f$\alpha^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) alpha

    !> \f$\beta^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) beta

    !> \f$\mu_{ij}\f$ of Fromme's algorithm
    COMPLEX(dp), ALLOCATABLE :: mu(:)

  END TYPE shift_system_type

CONTAINS

  !> Main driver routine of the algorithm.
  !!
  !! This subroutine implements the *Algorithm 2* of Fromme's paper.
  !!
  SUBROUTINE bicgstab(lmax, AA, bb, sigma)

    !> Dimensionality of the GMRES algorithm.
    INTEGER,     INTENT(IN) :: lmax

    !> Function pointer that applies the linear operator to a vector.
    INTERFACE
      SUBROUTINE AA(xx, Ax)
        USE kinds, ONLY: dp
        !> The input vector.
        COMPLEX(dp), INTENT(IN)  :: xx(:)
        !> The operator applied to the vector.
        COMPLEX(dp), INTENT(OUT) :: Ax(:)
      END SUBROUTINE AA
    END INTERFACE

    !> Right hand side of the linear equation.
    COMPLEX(dp), INTENT(IN) :: bb(:)

    !> Shifts \f$\sigma\f$ relative to the seed system.
    COMPLEX(dp), INTENT(IN) :: sigma(:)

    !> Size of the right-hand vector.
    INTEGER vec_size

    !> Counter for number of iterations.
    INTEGER iter

    !> Maximum number of iterations.
    INTEGER, PARAMETER :: max_iter = 10000

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the seed system.
    TYPE(seed_system_type)  seed_system

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the shifted system.
    TYPE(shift_system_type) shift_system(SIZE(sigma))

    ! determine dimension of vector
    vec_size = SIZE(bb)

    !
    ! initialization seed system
    !
    CALL init_seed(bb, seed_system)

    !
    ! initialization shifted systems
    !
    CALL init_shift(lmax, vec_size, sigma, shift_system)

    ! loop until solution is found
    DO iter = 1, max_iter

      !
      ! perform BiCG part (Algorithm 3)
      !
      CALL bicg_part(lmax, AA, seed_system, shift_system)

      !
      ! perform MR part (Algorithm 4)
      !

    END DO ! iter

    IF (iter > max_iter) THEN
      CALL errore(__FILE__, "BiCGstab algorithm did not converge in given&
                           & number of iterations", max_iter)
    END IF

    !
    ! destroy the allocated array in the types
    !
    CALL destroy_seed(seed_system)
    CALL destroy_shift(shift_system)

  END SUBROUTINE bicgstab

  !> Initialize the seed system.
  SUBROUTINE init_seed(bb, seed_system)

    !> Initial residual (right hand side of equation).
    COMPLEX(dp), INTENT(IN) :: bb(:)

    !> On output contains the initialized seed system.
    TYPE(seed_system_type), INTENT(OUT) :: seed_system

    !> Size of the vectors to initialize.
    INTEGER vec_size

    ! determine size of vector
    vec_size = SIZE(bb)

    ! allocate arrays of given size
    ALLOCATE(seed_system%u0(vec_size))
    ALLOCATE(seed_system%xx(vec_size))
    ALLOCATE(seed_system%r0(vec_size))
    ALLOCATE(seed_system%tilde_r0(vec_size))

    ! initialize arrays
    seed_system%u0 = 0
    seed_system%xx = 0

    ! initial residual is r0 = (b - A x) = b, because initial x = 0
    seed_system%r0 = bb
    seed_system%tilde_r0 = bb

    ! init the variables
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
    DEALLOCATE(seed_system%u0)
    DEALLOCATE(seed_system%xx)
    DEALLOCATE(seed_system%r0) 
    DEALLOCATE(seed_system%tilde_r0) 

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
    TYPE(shift_system_type), INTENT(OUT) :: shift_system(SIZE(sigma))

    !> counter over shifted systems
    INTEGER ishift

    !> dimension of the array mu
    INTEGER mu_size

    !> binomials \f$\begin{pmatrix}j - 1\\ i - 1\end{pmatrix}\f$
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

    !                     / j - 1 \        (j - 1)!
    ! evaluate binomials (         ) = -----------------
    !                     \ i - 1 /    (i - 1)! (j - i)!
    ! note: shift by one compared to Fromme, because Fortran starts indices at 1
    DO jj = 1, lmax

      ! the ij index starts at this value
      offset = jj * (jj - 1) / 2

      DO ii = 1, jj

        ! evaluate index in array
        ij = offset + ii

        ! trivial case i = 1: binomial = 1
        IF (ii == 1) THEN
          binomial(ij) = 1

        ! in general, the following is valid
        !  / j - 1 \     / j - 1 \  j - i + 1
        ! (         ) = (         ) ---------
        !  \ i - 1 /     \ i - 2 /    i - 1
        ELSE
          binomial(ij) = binomial(ij - 1) * (jj - ii + 1) / (ii - 1)

        END IF

      END DO ! i
    END DO ! j

    ! loop over all shifted systems
    DO ishift = 1, SIZE(sigma)

      ! allocate arrays of given size
      ALLOCATE(shift_system(ishift)%u0(vec_size))
      ALLOCATE(shift_system(ishift)%xx(vec_size))

      ! allocate array for mu
      ALLOCATE(shift_system(ishift)%mu(mu_size))

      ! initialize arrays to 0
      shift_system(ishift)%u0 = 0
      shift_system(ishift)%xx = 0

      ! initialize the variables
      shift_system(ishift)%sigma   = sigma(ishift)
      shift_system(ishift)%phi_old = 1.0
      shift_system(ishift)%phi     = 1.0
      shift_system(ishift)%theta   = 1.0

      ! construct sigma_pow array
      sigma_pow(0) = 1.0
      DO ii = 1, lmax - 1
        sigma_pow(ii) = sigma(ishift) * sigma_pow(ii - 1)
      END DO ! i

      ! construct mu_ij
      DO jj = 1, lmax

        ! the ij index starts at this value
        offset = jj * (jj - 1) / 2

        DO ii = 1, jj

          ! evaluate index in array
          ij = offset + ii

          !          /j - 1\       j - i
          ! mu_ij = (       ) sigma
          !          \i - 1/
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
      DEALLOCATE(shift_system(ishift)%u0)
      DEALLOCATE(shift_system(ishift)%xx)
      DEALLOCATE(shift_system(ishift)%mu)

    END DO ! ishift

  END SUBROUTINE destroy_shift

  !> The BiCG part of the BiCGstab algorithm.
  !!
  !! This subroutine implements the *Algorithm 3* of Fromme's paper.
  !!
  SUBROUTINE bicg_part(lmax, AA, seed_system, shift_system)

    !> Dimensionality of the GMRES algorithm.
    INTEGER, INTENT(IN) :: lmax

    !> Function pointer that applies the linear operator to a vector.
    INTERFACE
      SUBROUTINE AA(xx, Ax)
        USE kinds, ONLY: dp
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

    !> residuals of the seed system
    COMPLEX(dp), ALLOCATABLE :: rr(:,:)

    !> search direction of all systems
    !! index 0: seed; index 1-n: shifted
    COMPLEX(dp), ALLOCATABLE, TARGET :: uu_all(:,:,:)

    !> search direction of current system
    COMPLEX(dp), POINTER :: uu(:,:)

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

    ! BLAS routines
    COMPLEX(dp), EXTERNAL :: ZDOTC

    ! determine vector size
    vec_size = SIZE(seed_system%r0)

    ! determine number of shifted systems
    num_shift = SIZE(shift_system)

    ! allocate necessary arrays
    ALLOCATE(uu_all(vec_size, 0:lmax, 0:num_shift))
    ALLOCATE(rr(vec_size, 0:lmax))

    ! copy input to array
    uu_all(:, 0, 0) = seed_system%u0
    rr(:, 0) = seed_system%r0
    DO ishift = 1, ishift
      active => shift_system(ishift)
      uu_all(:, 1, ishift) = active%u0
    END DO ! ishift

    !
    ! Fromme's algorithm 3
    ! Note: the L?? refers to the line numbers given in Fromme's paper
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
      uu => uu_all(:,:,0)
      !
      ! L6: rho = (r_j, ~r_0),
      seed_system%rho = ZDOTC(vec_size, rr(:,jj), 1, seed_system%tilde_r0, 1)
      !     beta = alpha * rho / rho_old,
      seed_system%beta = seed_system%alpha * seed_system%rho / seed_system%rho_old
      !     rho_old = rho
      seed_system%rho_old = seed_system%rho
      !
      ! L7: loop over previous steps
      DO ii = 0, jj
        !
        ! L8: u_i = r_i - beta u_i
        CALL ZSCAL(vec_size, -seed_system%beta, uu(:,ii), 1)
        CALL ZAXPY(vec_size, one, rr(:,ii), 1, uu(:,ii), 1)
        !
      END DO ! i
      !
      ! L10: u_j+1 = A u_j
      CALL AA(uu(:,jj), uu(:, jj + 1))
      !
      ! L11: alpha = rho / (u_j+1, ~r_0)
      seed_system%alpha = seed_system%rho &
                        / ZDOTC(vec_size, uu(:, jj + 1), 1, seed_system%tilde_r0, 1)
      !
      ! shifted system
      !
      ! loop over all shifted systems
      DO ishift = 1, num_shift
        !
        active => shift_system(ishift)
        uu => uu_all(:,:,ishift)
        !
        ! L13: phi_new^ = (1 + alpha sigma) phi^ + alpha beta / alpha_old (phi_old^ - phi^)
        active%phi_new = (one + seed_system%alpha * active%sigma) * active%phi &
                       + seed_system%alpha * seed_system%beta / seed_system%alpha_old &
                       * (active%phi_old - active%phi)
        !      beta^ = (phi_old^ / phi^)**2 beta
        active%beta = (active%phi_old / active%phi)**2 * seed_system%beta
        !      alpha^ = (phi_old^ / phi^) alpha
        active%alpha = (active%phi_old / active%phi) * seed_system%alpha
        !
        ! evaluate 1 / (theta^ phi^) used in L15 and L17
        factor = 1.0 / (active%theta * active%phi)
        !
        ! L14: loop over previous steps
        DO ii = 0, jj
          !
          ! L15: u_i^ = r_i / (theta^ phi^) - beta^ u_i^
          CALL ZSCAL(vec_size, -active%beta, uu(:,ii), 1)
          CALL ZAXPY(vec_size, factor, rr(:,ii), 1, uu(:,ii), 1)
          !
        END DO ! i
        !
        ! L17: x^ = x^ + alpha^ u0^
        CALL ZAXPY(vec_size, active%alpha, uu(:,0), 1, active%xx, 1)
        !
        ! note - we shifted update of alpha after the loop to allow
        !        for systems with multiple shifts
        ! L18: phi_old^ = phi^
        active%phi_old = active%phi
        !      phi^ = phi_new^
        active%phi = active%phi_new
        !
        ! note: we can use factor here, because phi_old = phi in L18
        ! L19: u_j+1^ = r_j / (theta^ phi_old^)
        CALL ZCOPY(vec_size, rr(:,jj), 1, uu(:, jj + 1), 1)
        CALL ZSCAL(vec_size, factor, uu(:, jj + 1), 1)
        !
      END DO ! ishift
      !
      ! back to seed system
      uu => uu_all(:,:,0)
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
        CALL ZAXPY(vec_size, -seed_system%alpha, uu(:, ii + 1), 1, rr(:,ii), 1)
        !
      END DO ! i
      !
      ! L24: r_j+1 = A r_j
      CALL AA(rr(:,jj), rr(:, jj + 1))
      !
      ! L25: x = x + alpha u0
      CALL ZAXPY(vec_size, seed_system%alpha, uu(:,0), 1, seed_system%xx, 1)
      !
      ! now update shifted systems again
      DO ishift = 1, num_shift
        !
        active => shift_system(ishift)
        uu => uu_all(:,:,ishift)
        !
        ! L26: u_j+1^ = (u_j+1^ - r_j / (theta^ phi^)) / alpha^ - sigma u_j^
        factor = -1.0 / (active%theta * active%phi)
        CALL ZAXPY(vec_size, factor, rr(:,jj), 1, uu(:, jj + 1), 1)
        CALL ZSCAL(vec_size, 1.0 / active%alpha, uu(:, jj + 1), 1)
        CALL ZAXPY(vec_size, -active%sigma, uu(:,jj), 1, uu(:, jj + 1), 1)
        !
      END DO ! ishift
      !
    END DO ! j

    ! free memory
    DEALLOCATE(uu)
    DEALLOCATE(rr)

  END SUBROUTINE bicg_part

END MODULE bicgstab_module

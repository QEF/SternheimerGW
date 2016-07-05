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

    !> \f$\rho_{\text{old}}\f$ of Fromme's algorithm
    COMPLEX(dp) rho_old

    !> \f$\alpha_{\text{old}}\f$ of Fromme's algorithm
    COMPLEX(dp) alpha_old

    !> \f$\alpha\f$ of Fromme's algorithm
    COMPLEX(dp) alpha

    !> \f$\omega\f$ of Fromme's algorithm
    COMPLEX(dp) omega
  
  END TYPE seed_system_type

  !> Variables and arrays used for the shifted system.
  TYPE shift_system_type

    !> search direction
    COMPLEX(dp), ALLOCATABLE :: u0(:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

    !> \f$\phi_{\text{old}}^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) phi_old

    !> \f$\phi^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) phi

    !> \f$\vartheta^\sigma\f$ of Fromme's algorithm
    COMPLEX(dp) theta

    !> \f$\mu_{ij}\f$ of Fromme's algorithm
    COMPLEX(dp), ALLOCATABLE :: mu(:)

  END TYPE shift_system_type

CONTAINS

  !> Main driver routine of the algorithm.
  !!
  !! This subroutine implements the *Algorithm 2* of Fromme's paper.
  !!
  SUBROUTINE bicgstab(lmax, bb, sigma)

    !> Dimensionality of the GMRES algorithm.
    INTEGER,     INTENT(IN) :: lmax

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
      CALL bicg_part(lmax, seed_system, shift_system)

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
  SUBROUTINE bicg_part(lmax, seed_system, shift_system)

    !> Dimensionality of the GMRES algorithm.
    INTEGER, INTENT(IN) :: lmax

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the seed system. Updated
    !! to the next best guess at the end of the routine.
    TYPE(seed_system_type),  INTENT(INOUT) :: seed_system

    !> Contains the current best guess for the solution of the linear
    !! system and the corresponding residual in the shifted system. Updated
    !! to the next best guess at the end of the routine.
    TYPE(shift_system_type), INTENT(INOUT) :: shift_system(:)

  END SUBROUTINE bicg_part

END MODULE bicgstab_module

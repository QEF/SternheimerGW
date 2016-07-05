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
!> Shifted BiCGStab(l) solver for linear equation \f$(A + \lambda I) x = b\f$.
!!
!! Implements the shifted BiCGStab(l) algorithm according to the paper of
!! Fromme, Computing **70**, 87 (2003). The general idea is that the matrix
!! \f$A\f$ and \f$A + \lambda I\f$ span the same Krylov subspace. Hence, we can
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
    COMPLEX(dp), ALLOCATABLE :: uu(:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

    !> residual vector
    COMPLEX(dp), ALLOCATABLE :: rr(:)

  END TYPE seed_system_type

  !> Variables and arrays used for the shifted system.
  TYPE shift_system_type

    !> search direction
    COMPLEX(dp), ALLOCATABLE :: uu(:)

    !> approximate solution to the linear system
    COMPLEX(dp), ALLOCATABLE :: xx(:)

  END TYPE shift_system_type

CONTAINS

  !> Main driver routine of the algorithm.
  !!
  !! This subroutine implements the *Algorithm 2* of Fromme's paper.
  !!
  SUBROUTINE bicgstab(lmax, bb)

    !> Dimensionality of the GMRES algorithm.
    INTEGER,     INTENT(IN) :: lmax

    !> Right hand side of the linear equation.
    COMPLEX(dp), INTENT(IN) :: bb(:)

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
    TYPE(shift_system_type) shift_system(1)

    ! determine dimension of vector
    vec_size = SIZE(bb)

    !
    ! initialization seed system
    !
    CALL init_seed(bb, seed_system)

    !
    ! initialization shifted systems
    !
    CALL init_shift(vec_size, shift_system)

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
    ALLOCATE(seed_system%uu(vec_size))
    ALLOCATE(seed_system%xx(vec_size))
    ALLOCATE(seed_system%rr(vec_size))

    ! initialize arrays
    seed_system%uu = 0
    seed_system%xx = 0

    ! initial residual is r = (b - A x) = b, because initial x = 0
    seed_system%rr = bb

  END SUBROUTINE init_seed

  !> Deallocate the arrays in the seed system type.
  SUBROUTINE destroy_seed(seed_system)

    !> Deallocate the arrays of this system.
    TYPE(seed_system_type), INTENT(INOUT) :: seed_system

    ! deallocate the arrays
    DEALLOCATE(seed_system%uu)
    DEALLOCATE(seed_system%xx)
    DEALLOCATE(seed_system%rr) 

  END SUBROUTINE destroy_seed

  !> Initialize the shifted systems.
  SUBROUTINE init_shift(vec_size, shift_system)

    !> Size of the vectors to initialize.
    INTEGER, INTENT(IN) :: vec_size

    !> On output contains the initialized shifted systems.
    TYPE(shift_system_type), INTENT(OUT) :: shift_system(:)

    !> counter over shifted systems
    INTEGER ishift

    ! loop over all shifted systems
    DO ishift = 1, SIZE(shift_system)

      ! allocate arrays of given size
      ALLOCATE(shift_system(ishift)%uu(vec_size))
      ALLOCATE(shift_system(ishift)%xx(vec_size))

      ! initialize arrays to 0
      shift_system(ishift)%uu = 0
      shift_system(ishift)%xx = 0

    END DO ! ishift

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

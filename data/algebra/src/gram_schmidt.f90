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
!> Implements the modified Gram-Schmidt orthonormalization
MODULE gram_schmidt_module

  IMPLICIT NONE

CONTAINS

  !> Use the modified Gram-Schmidt algorithm to generate a set of orthonormal vectors
  !!
  !! This subroutine is written to allow for extension of a previously generated set
  !! of vectors by additional ones added to the back. It also allows to transform a
  !! second set of vectors with the same transformation at the same time.
  SUBROUTINE gram_schmidt(first, basis, vector)

    USE kinds,       ONLY: dp
    USE norm_module, ONLY: znorm => norm

    !> the first vector that is not orthonormal to the rest
    !! (set to 1 to orthonormalize all)
    INTEGER,     INTENT(IN)    :: first

    !> the vectors intended to form an orthonormal basis
    COMPLEX(dp), INTENT(INOUT) :: basis(:,:)

    !> the vectors which should be transformed according to the same transformation
    COMPLEX(dp), INTENT(INOUT), OPTIONAL :: vector(:,:)

    !> the dimensionality of the problem
    INTEGER vec_size

    !> the number of basis functions
    INTEGER num_basis

    !> counter on the basis functions
    INTEGER ibasis, jbasis

    !> flag set when an additional vector set is provided
    LOGICAL transform

    !> the result of the dot product
    COMPLEX(dp) norm

    !> complex value of 1
    COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

    !> LAPACK function to evaluate the 2-norm
    REAL(dp),    EXTERNAL :: DNRM2

    !> LAPACK function to evaluate the dot product
    COMPLEX(dp), EXTERNAL :: ZDOTC

    ! initialize helper variables
    vec_size  = SIZE(basis, 1)
    num_basis = SIZE(basis, 2)
    transform = PRESENT(vector)

    ! sanity check
    IF (vec_size < num_basis) &
      CALL errore(__FILE__, "number of vectors larger than dimension", 1)
    IF (transform) THEN
      IF (ANY(SHAPE(basis) /= SHAPE(vector))) &
        CALL errore(__FILE__, "basis and vector array inconsistent", 1)
    END IF

    !!
    !! 1. Orthogonalize to the existing basis \f$\vert w_i\rangle \f$
    !!    \f{equation}{
    !!      \vert v_i \rangle = \vert v_i \rangle - \sum_j \vert w_j \rangle \langle w_j \vert v_i \rangle
    !!    \f}
    !!
    DO ibasis = first, num_basis
      !
      DO jbasis = 1, first - 1
        !
        ! v_i = v_i - (w_j, v_i) w_j
        norm = ZDOTC(vec_size, basis(:, jbasis), 1, basis(:, ibasis), 1)
        CALL ZAXPY(vec_size, -norm, basis(:, jbasis), 1, basis(:, ibasis), 1)
        !
        ! t_i = t_i - (w_j, v_i) t_j
        IF (transform) CALL ZAXPY(vec_size, -norm, vector(:, jbasis), 1, vector(:, ibasis), 1)
        !
      END DO
      !
    END DO ! ibasis

    !!
    !! 2. loop over all new vectors \f$\vert v_i\rangle\f$
    !!
    DO ibasis = first, num_basis
      !!
      !! 3. normalize the current vector
      norm = one / znorm(basis(:, ibasis))
      CALL ZSCAL(vec_size, norm, basis(:, ibasis), 1)
      !
      IF (transform) CALL ZSCAL(vec_size, norm, vector(:, ibasis), 1)
      !!
      !! 4. orthogonalize all other vector with respect to current one
      !!
      DO jbasis = ibasis + 1, num_basis
        !
        ! v_j = v_j - (v_i, v_j) v_i
        norm = ZDOTC(vec_size, basis(:, ibasis), 1, basis(:, jbasis), 1)
        CALL ZAXPY(vec_size, -norm, basis(:, ibasis), 1, basis(:, jbasis), 1)
        !
        ! t_j = t_j - (v_i, v_j) v_i
        IF (transform) CALL ZAXPY(vec_size, -norm, vector(:, ibasis), 1, vector(:, jbasis), 1)
        !
      END DO ! jbasis
      !
    END DO ! ibasis

  END SUBROUTINE gram_schmidt

END MODULE gram_schmidt_module

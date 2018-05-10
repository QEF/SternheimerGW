! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

MODULE lapack_module

  USE constant_module, ONLY: dp

  IMPLICIT NONE

  INTEGER, PARAMETER :: determine = -1
  INTEGER, PARAMETER :: no_error = 0

  PRIVATE
  PUBLIC svd, eigenvalue, rational_vector, no_error

  INTERFACE svd
    MODULE PROCEDURE svd_value, svd_full
  END INTERFACE svd

  TYPE svd_type

    COMPLEX(dp), ALLOCATABLE :: copy(:,:)
    INTEGER, ALLOCATABLE :: int_(:)
    REAL(dp), ALLOCATABLE :: real_(:)
    COMPLEX(dp), ALLOCATABLE :: complex_(:)

    REAL(dp), ALLOCATABLE :: sigma(:)
    COMPLEX(dp), ALLOCATABLE :: U_matrix(:,:)
    COMPLEX(dp), ALLOCATABLE :: V_dagger_matrix(:,:)

    INTEGER info

  END TYPE svd_type

  TYPE rational_vector
    COMPLEX(dp), ALLOCATABLE :: numerator(:)
    COMPLEX(dp), ALLOCATABLE :: denominator(:)
  END TYPE rational_vector

  TYPE eigenvalue_type
    INTEGER :: size
    COMPLEX(dp), ALLOCATABLE :: A_matrix(:,:)
    COMPLEX(dp), ALLOCATABLE :: B_matrix(:,:)
    REAL(dp), ALLOCATABLE :: real_(:)
    COMPLEX(dp), ALLOCATABLE :: complex_(:)
    COMPLEX(dp), ALLOCATABLE :: numerator(:)
    COMPLEX(dp), ALLOCATABLE :: denominator(:)
    INTEGER info
  END TYPE eigenvalue_type

CONTAINS

  SUBROUTINE svd_value(matrix, sigma, info)

    COMPLEX(dp), INTENT(IN) :: matrix(:,:)
    REAL(dp), INTENT(OUT), ALLOCATABLE :: sigma(:)
    INTEGER, INTENT(OUT) :: info

    COMPLEX(dp), ALLOCATABLE :: U_dummy(:,:), V_dummy(:,:) 

    CHARACTER, PARAMETER :: jobz = 'N'

    CALL svd_general(jobz, matrix, sigma, U_dummy, V_dummy, info)

  END SUBROUTINE svd_value

  SUBROUTINE svd_full(matrix, sigma, U_matrix, V_dagger_matrix, info)

    COMPLEX(dp), INTENT(IN) :: matrix(:,:)
    REAL(dp), INTENT(OUT), ALLOCATABLE :: sigma(:)
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: U_matrix(:,:)
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: V_dagger_matrix(:,:)
    INTEGER, INTENT(OUT) :: info

    CHARACTER, PARAMETER :: jobz = 'A'

    CALL svd_general(jobz, matrix, sigma, U_matrix, V_dagger_matrix, info)

  END SUBROUTINE svd_full

  SUBROUTINE svd_general(jobz, matrix, sigma, U_matrix, V_dagger_matrix, info)

    CHARACTER, INTENT(IN) :: jobz
    COMPLEX(dp), INTENT(IN) :: matrix(:,:)
    REAL(dp), INTENT(OUT), ALLOCATABLE :: sigma(:)
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: U_matrix(:,:)
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: V_dagger_matrix(:,:)
    INTEGER, INTENT(OUT) :: info

    TYPE(svd_type) :: work

    IF (trivial_case(matrix)) THEN

      CALL set_default(matrix, work)

    ELSE

      CALL allocate_for_svd(jobz, matrix, work)

      CALL perform_svd(jobz, work)

    END IF

    CALL copy_svd_to_output(work, sigma, U_matrix, V_dagger_matrix, info)

  END SUBROUTINE svd_general

  LOGICAL FUNCTION trivial_case(matrix)

    COMPLEX(dp), INTENT(IN) :: matrix(:,:)

    trivial_case = SIZE(matrix) == 0

  END FUNCTION trivial_case

  SUBROUTINE set_default(matrix, work)

    COMPLEX(dp), INTENT(IN) :: matrix(:,:)
    TYPE(svd_type), INTENT(OUT) :: work

    ALLOCATE(work%sigma(0))
    CALL identity_matrix(SIZE(matrix, 1), work%U_matrix)
    CALL identity_matrix(SIZE(matrix, 2), work%V_dagger_matrix)
    work%info = no_error

  END SUBROUTINE set_default

  SUBROUTINE identity_matrix(array_size, matrix)

    USE constant_module, ONLY: c_zero, c_one

    INTEGER, INTENT(IN) :: array_size
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: matrix(:,:)

    INTEGER ii

    ALLOCATE(matrix(array_size, array_size))
    matrix = c_zero

    DO ii = 1, array_size
      matrix(ii, ii) = c_one
    END DO

  END SUBROUTINE identity_matrix

  SUBROUTINE allocate_for_svd(jobz, matrix, work)

    CHARACTER, INTENT(IN) :: jobz
    COMPLEX(dp), INTENT(IN) :: matrix(:,:)
    TYPE(svd_type), INTENT(OUT) :: work

    INTEGER :: min_row_col, max_row_col, opt_size

    min_row_col = MINVAL(SHAPE(matrix))
    max_row_col = MAXVAL(SHAPE(matrix))

    ALLOCATE(work%sigma(min_row_col))
    ALLOCATE(work%copy(SIZE(matrix, 1), SIZE(matrix, 2)))
    CALL ZCOPY(SIZE(matrix), matrix, 1, work%copy, 1)

    ALLOCATE(work%int_(8 * min_row_col))

    SELECT CASE (jobz)
    CASE ('N')
      ALLOCATE(work%real_(7 * min_row_col))
      ALLOCATE(work%U_matrix(1, 1))
      ALLOCATE(work%V_dagger_matrix(1, 1))

    CASE ('A')
      ALLOCATE(work%real_(min_row_col * (2 * max_row_col + 3 * min_row_col + 7)))
      ALLOCATE(work%U_matrix(SIZE(matrix, 1), SIZE(matrix, 1)))
      ALLOCATE(work%V_dagger_matrix(SIZE(matrix, 2), SIZE(matrix, 2)))

    END SELECT

    opt_size = svd_work_size(jobz, work)
    ALLOCATE(work%complex_(opt_size))

  END SUBROUTINE allocate_for_svd

  INTEGER FUNCTION svd_work_size(jobz, work) RESULT(res)

    CHARACTER, INTENT(IN) :: jobz
    TYPE(svd_type), INTENT(IN) :: work

    COMPLEX(dp) opt_size
 
    CALL ZGESDD(jobz, SIZE(work%copy, 1), SIZE(work%copy, 2), work%copy, SIZE(work%copy, 1), work%sigma, &
                work%U_matrix, SIZE(work%U_matrix, 1), work%V_dagger_matrix, SIZE(work%V_dagger_matrix, 1), &
                opt_size, determine, work%real_, work%int_, work%info)

    res = NINT(ABS(opt_size))

  END FUNCTION svd_work_size

  SUBROUTINE perform_svd(jobz, work)

    CHARACTER, INTENT(IN) :: jobz
    TYPE(svd_type), INTENT(INOUT) :: work

    IF (work%info /= no_error) RETURN

    CALL ZGESDD(jobz, SIZE(work%copy, 1), SIZE(work%copy, 2), work%copy, SIZE(work%copy, 1), work%sigma, &
                work%U_matrix, SIZE(work%U_matrix, 1), work%V_dagger_matrix, SIZE(work%V_dagger_matrix, 1), &
                work%complex_, SIZE(work%complex_), work%real_, work%int_, work%info)

  END SUBROUTINE perform_svd

  SUBROUTINE copy_svd_to_output(work, sigma, U_matrix, V_dagger_matrix, info)

    TYPE(svd_type), INTENT(INOUT) :: work
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: sigma(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: U_matrix(:,:), V_dagger_matrix(:,:)
    INTEGER, INTENT(OUT) :: info

    info = work%info
    CALL MOVE_ALLOC(work%sigma, sigma)
    CALL MOVE_ALLOC(work%U_matrix, U_matrix)
    CALL MOVE_ALLOC(work%V_dagger_matrix, V_dagger_matrix)

  END SUBROUTINE copy_svd_to_output


  ! A x = lambda B x
  SUBROUTINE eigenvalue(A_matrix, B_matrix, lambda, info)
    !
    COMPLEX(dp), INTENT(IN) :: A_matrix(:,:), B_matrix(:,:)
    TYPE(rational_vector), INTENT(OUT) :: lambda
    INTEGER, INTENT(OUT) :: info
    !
    CHARACTER, PARAMETER :: jobv = 'N'
    TYPE(eigenvalue_type) work
    !
    CALL allocate_for_eigenvalue(jobv, A_matrix, B_matrix, work, info)
    IF (info /= 0) RETURN
    CALL perform_eigenvalue(jobv, work, info)
    IF (info /= 0) RETURN
    CALL copy_eigenvalue_to_output(work, lambda)
    !
  END SUBROUTINE eigenvalue


  SUBROUTINE allocate_for_eigenvalue(jobv, A_matrix, B_matrix, work, info)
    !
    USE array_module, ONLY: allocate_copy_from_to
    !
    CHARACTER, INTENT(IN) :: jobv
    COMPLEX(dp), INTENT(IN) :: A_matrix(:,:), B_matrix(:,:)
    TYPE(eigenvalue_type), INTENT(OUT) :: work
    INTEGER, INTENT(OUT) :: info
    !
    INTEGER opt_size
    !
    work%size = SIZE(A_matrix, 1)
    CALL allocate_copy_from_to(A_matrix, work%A_matrix)
    CALL allocate_copy_from_to(B_matrix, work%B_matrix)
    ALLOCATE(work%numerator(work%size), work%denominator(work%size))
    opt_size = eigenvalue_work_size(jobv, work, info)
    IF (info /= no_error) RETURN
    ALLOCATE(work%complex_(opt_size), work%real_(8 * work%size))
    ! 
  END SUBROUTINE allocate_for_eigenvalue


  INTEGER FUNCTION eigenvalue_work_size(jobv, work, info) RESULT(res)
    !
    CHARACTER, INTENT(IN) :: jobv
    TYPE(eigenvalue_type), INTENT(IN) :: work
    INTEGER, INTENT(OUT) :: info
    !
    COMPLEX(dp) opt_size, dummy
    !
    CALL ZGGEV(jobv, jobv, work%size, work%A_matrix, work%size, work%B_matrix, work%size, &
               work%numerator, work%denominator, dummy, 1, dummy, 1, &
               opt_size, determine, work%real_, info)
    res = NINT(ABS(opt_size))
    !
  END FUNCTION eigenvalue_work_size


  SUBROUTINE perform_eigenvalue(jobv, work, info)
    !
    CHARACTER, INTENT(IN) :: jobv
    TYPE(eigenvalue_type), INTENT(IN) :: work
    INTEGER, INTENT(OUT) :: info
    !
    COMPLEX(dp) dummy
    !
    CALL ZGGEV(jobv, jobv, work%size, work%A_matrix, work%size, work%B_matrix, work%size, &
               work%numerator, work%denominator, dummy, 1, dummy, 1, &
               work%complex_, SIZE(work%complex_), work%real_, info)
    !
  END SUBROUTINE perform_eigenvalue


  SUBROUTINE copy_eigenvalue_to_output(work, lambda)
    !
    TYPE(eigenvalue_type), INTENT(INOUT) :: work
    TYPE(rational_vector), INTENT(OUT) :: lambda
    !
    CALL MOVE_ALLOC(work%numerator, lambda%numerator)
    CALL MOVE_ALLOC(work%denominator, lambda%denominator)
    !
  END SUBROUTINE copy_eigenvalue_to_output

END MODULE lapack_module

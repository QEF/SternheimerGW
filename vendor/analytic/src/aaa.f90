! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

MODULE aaa_module

  USE constant_module, ONLY: dp

  IMPLICIT NONE

  INTEGER, PARAMETER :: no_error = 0
  INTEGER, PARAMETER :: input_error = 1

  INTEGER, PARAMETER :: no_restriction = -1

  TYPE aaa_approx
    COMPLEX(dp), ALLOCATABLE :: position(:)
    COMPLEX(dp), ALLOCATABLE :: value(:)
    COMPLEX(dp), ALLOCATABLE :: weight(:)
  END TYPE aaa_approx

  TYPE pole_residual
    COMPLEX(dp), ALLOCATABLE :: pole(:)
    COMPLEX(dp), ALLOCATABLE :: residual(:)
  END TYPE pole_residual

  TYPE aaa_work_type
    REAL(dp) thres
    INTEGER max_point
    COMPLEX(dp), ALLOCATABLE :: zz(:), ff(:), fit(:)
    COMPLEX(dp), ALLOCATABLE :: loewner_matrix(:,:)
    LOGICAL, ALLOCATABLE :: is_support_point(:)
  END TYPE aaa_work_type

  TYPE mask_type
    INTEGER num_row, num_col
    LOGICAL, ALLOCATABLE :: matrix(:,:)
  END TYPE mask_type

  TYPE distance_type
    INTEGER index
    REAL(dp) value
  END TYPE distance_type

  PRIVATE
  PUBLIC aaa_approx, aaa_generate, aaa_evaluate, aaa_pole_residual, &
    pole_residual, mask_type, no_error, input_error, no_restriction, &
    select_point_with_maximum_error, construct_cauchy_matrix, &
    construct_loewner_matrix, extract_submatrix, and_mask_matrix

CONTAINS

  SUBROUTINE aaa_generate(thres, max_point, zz, ff, analytic_cont, info)
    !
    REAL(dp), INTENT(IN) :: thres
    INTEGER, INTENT(IN) :: max_point
    COMPLEX(dp), INTENT(IN) :: zz(:), ff(:)
    TYPE(aaa_approx), INTENT(OUT) :: analytic_cont
    INTEGER, INTENT(OUT) :: info
    !
    TYPE(aaa_work_type) :: work
    !
    CALL check_input_for_error(zz, ff, info)
    IF (info /= no_error) RETURN
    CALL setup_work_type(thres, max_point, zz, ff, work)
    CALL determine_analytic_cont(work, analytic_cont, info)
    !
  END SUBROUTINE aaa_generate


  SUBROUTINE aaa_evaluate(analytic_cont, zz, ff)
    !
    USE constant_module, ONLY: eps14
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    COMPLEX(dp), INTENT(IN) :: zz(:)
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: ff(:)
    !
    INTEGER ii
    TYPE(distance_type) dist
    !
    CALL evaluate_analytic_cont(analytic_cont, zz, ff)
    DO ii = 1, SIZE(zz)
      dist = distance_of_closest_support_point(zz(ii), analytic_cont%position)
      IF (dist%value < eps14) THEN
        ff(ii) = analytic_cont%value(dist%index)
      END IF
    END DO
    !
  END SUBROUTINE aaa_evaluate


  SUBROUTINE aaa_pole_residual(analytic_cont, pole_res, info)
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    TYPE(pole_residual), INTENT(OUT) :: pole_res
    INTEGER, INTENT(OUT) :: info
    !
    CALL find_pole(analytic_cont, pole_res%pole, info)
    IF (info /= no_error) RETURN
    CALL calculate_residual(analytic_cont, pole_res%pole, pole_res%residual)
    !
  END SUBROUTINE aaa_pole_residual


  SUBROUTINE check_input_for_error(zz, ff, info)
    !
    COMPLEX(dp), INTENT(IN) :: zz(:), ff(:)
    INTEGER, INTENT(OUT) :: info
    !
    info = no_error
    IF (SIZE(zz) /= SIZE(ff)) THEN
      info = input_error
    END IF
    !
  END SUBROUTINE check_input_for_error


  SUBROUTINE setup_work_type(relative_threshold, max_point, zz, ff, work)
    !
    USE array_module,  ONLY: allocate_init_to, allocate_copy_from_to
    USE assert_module, ONLY: assert
    !
    REAL(dp), INTENT(IN) :: relative_threshold
    INTEGER, INTENT(IN) :: max_point
    COMPLEX(dp), INTENT(IN) :: zz(:), ff(:)
    TYPE(aaa_work_type), INTENT(OUT) :: work
    !
    INTEGER array_size
    COMPLEX(dp) avg_ff
    !
    array_size = SIZE(ff)
    CALL assert(array_size > 0, "cannot calculate average for empty array")
    avg_ff = SUM(ff) / array_size
    !
    work%thres = absolute_threshold(relative_threshold, ff)
    IF (max_point == no_restriction) THEN
      work%max_point = array_size
    ELSE
      work%max_point = max_point
    END IF
    CALL allocate_copy_from_to(zz, work%zz)
    CALL allocate_copy_from_to(ff, work%ff)
    CALL allocate_init_to(array_size, avg_ff, work%fit)
    CALL allocate_init_to(array_size, .FALSE., work%is_support_point) 
    CALL construct_loewner_matrix(zz, ff, work%loewner_matrix)
    !
  END SUBROUTINE setup_work_type


  REAL(dp) FUNCTION absolute_threshold(relative_threshold, ff)
    !
    REAL(dp), INTENT(IN) :: relative_threshold
    COMPLEX(dp), INTENT(IN) :: ff(:)
    !
    absolute_threshold = relative_threshold * MAXVAL(ABS(ff))
    !
  END FUNCTION absolute_threshold


  SUBROUTINE determine_analytic_cont(work, analytic_cont, info)
    !
    TYPE(aaa_work_type), INTENT(INOUT) :: work
    TYPE(aaa_approx), INTENT(OUT) :: analytic_cont
    INTEGER, INTENT(OUT) :: info
    !
    DO
      CALL update_support_point(work)
      CALL create_current_best_analytic_cont(work, analytic_cont, info)
      IF (info /= no_error) RETURN
      CALL update_fit(analytic_cont, work)
      IF (converged(work).OR.max_step_reached(analytic_cont, work)) EXIT
    END DO
    !
  END SUBROUTINE determine_analytic_cont


  SUBROUTINE update_support_point(work)
    !
    TYPE(aaa_work_type), INTENT(INOUT) :: work
    INTEGER new_point
    !
    new_point = select_point_with_maximum_error(work%ff, work%fit)
    work%fit(new_point) = work%ff(new_point)
    work%is_support_point(new_point) = .TRUE.
    !
  END SUBROUTINE update_support_point


  SUBROUTINE create_current_best_analytic_cont(work, analytic_cont, info)
    !
    USE array_module, ONLY: allocate_copy_from_to
    !
    TYPE(aaa_work_type), INTENT(IN) :: work
    TYPE(aaa_approx), INTENT(OUT) :: analytic_cont
    INTEGER, INTENT(OUT) :: info
    !
    CALL allocate_copy_from_to(PACK(work%zz, work%is_support_point), analytic_cont%position)
    CALL allocate_copy_from_to(PACK(work%ff, work%is_support_point), analytic_cont%value)
    CALL determine_weight(work, analytic_cont%weight, info)
    !
  END SUBROUTINE create_current_best_analytic_cont


  SUBROUTINE determine_weight(work, weight, info)
    !
    USE array_module, ONLY: allocate_copy_from_to
    USE lapack_module, ONLY: svd
    !
    TYPE(aaa_work_type), INTENT(IN) :: work
    COMPLEX(dp), ALLOCATABLE :: weight(:)
    INTEGER, INTENT(OUT) :: info
    !
    TYPE(mask_type) mask_t
    COMPLEX(dp), ALLOCATABLE :: loewner_submatrix(:,:)
    REAL(dp), ALLOCATABLE :: sigma(:)
    COMPLEX(dp), ALLOCATABLE :: U_matrix(:,:), V_dagger_matrix(:,:)
    !
    CALL AND_mask_matrix(.NOT.work%is_support_point, work%is_support_point, mask_t)
    CALL extract_submatrix(work%loewner_matrix, mask_t, loewner_submatrix)
    CALL svd(loewner_submatrix, sigma, U_matrix, V_dagger_matrix, info)
    IF (info /= no_error) RETURN
    CALL allocate_copy_from_to(CONJG(V_dagger_matrix(SIZE(V_dagger_matrix, 1),:)), weight)
    !
  END SUBROUTINE determine_weight


  SUBROUTINE update_fit(analytic_cont, work)
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    TYPE(aaa_work_type), INTENT(INOUT) :: work
    !
    COMPLEX(dp), ALLOCATABLE :: ff(:)
    !
    CALL evaluate_analytic_cont(analytic_cont, work%zz, ff)
    WHERE (.NOT.work%is_support_point)
      work%fit = ff
    END WHERE
    !
  END SUBROUTINE update_fit


  SUBROUTINE evaluate_analytic_cont(analytic_cont, zz, ff)
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    COMPLEX(dp), INTENT(IN) :: zz(:)
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: ff(:)
    !
    COMPLEX(dp), ALLOCATABLE :: cauchy_matrix(:,:)
    COMPLEX(dp), ALLOCATABLE :: numerator(:), denominator(:)
    !
    ALLOCATE(ff(SIZE(zz)), numerator(SIZE(zz)), denominator(SIZE(zz)))
    CALL construct_cauchy_matrix(zz, analytic_cont%position, cauchy_matrix)
    numerator = MATMUL(cauchy_matrix, analytic_cont%weight * analytic_cont%value)
    denominator = MATMUL(cauchy_matrix, analytic_cont%weight)
    ff = numerator / denominator
    !
  END SUBROUTINE evaluate_analytic_cont


  TYPE(distance_type) FUNCTION distance_of_closest_support_point(zz, support_point) RESULT(dist)
    !
    COMPLEX(dp), INTENT(IN) :: zz
    COMPLEX(dp), INTENT(IN) :: support_point(:)
    !
    REAL(dp), ALLOCATABLE :: distance(:)
    !
    ALLOCATE(distance(SIZE(support_point)))
    distance = ABS(zz - support_point)
    dist%index = MINLOC(distance, DIM = 1)
    dist%value = distance(dist%index)
    !
  END FUNCTION distance_of_closest_support_point


  LOGICAL FUNCTION converged(work)
    !
    TYPE(aaa_work_type), INTENT(IN) :: work
    !
    converged = ALL(ABS(work%fit - work%ff) <= work%thres)
    !
  END FUNCTION converged


  LOGICAL FUNCTION max_step_reached(analytic_cont, work)
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    TYPE(aaa_work_type), INTENT(IN) :: work
    !
    max_step_reached = SIZE(analytic_cont%position) >= work%max_point
    !
  END FUNCTION max_step_reached


  INTEGER FUNCTION select_point_with_maximum_error(actual, fit) RESULT (point)
    !
    USE assert_module,   ONLY: assert
    !
    COMPLEX(dp), INTENT(IN) :: actual(:), fit(:)
    !
    CALL assert(SIZE(actual) == SIZE(fit), "compare arrays or same dimensionality.")
    point = MAXLOC(ABS(actual - fit), 1)
    !
  END FUNCTION select_point_with_maximum_error


  ! Cauchy matrix: a_ij = 1 / (x_i - y_j)
  PURE SUBROUTINE construct_Cauchy_matrix(xx, yy, cauchy_matrix)
    !
    USE constant_module, ONLY: one, eps14
    !
    COMPLEX(dp), INTENT(IN) :: xx(:), yy(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: cauchy_matrix(:,:)
    !
    INTEGER ix, iy
    COMPLEX(dp) denominator
    !
    ALLOCATE(cauchy_matrix(SIZE(xx), SIZE(yy)))
    DO iy = 1, SIZE(yy)
      DO ix = 1, SIZE(xx)
        denominator = xx(ix) - yy(iy)
        IF (ABS(denominator) <= eps14) denominator = eps14
        cauchy_matrix(ix, iy) = one / denominator
      END DO
    END DO
    !
  END SUBROUTINE construct_Cauchy_matrix


  ! Loewner matrix: a_ij = (f_i - f_j) / (z_i _ z_j)
  SUBROUTINE construct_Loewner_matrix(zz, ff, loewner_matrix)
    !
    USE assert_module,   ONLY: assert
    USE constant_module, ONLY: c_zero
    !
    COMPLEX(dp), INTENT(IN) :: zz(:), ff(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: loewner_matrix(:,:)
    !
    INTEGER ii, jj
    !
    CALL assert(SIZE(zz) == SIZE(ff), "function points and values must have same dimension.")
    ALLOCATE(loewner_matrix(SIZE(zz), SIZE(zz)))
    DO jj = 1, SIZE(zz)
      loewner_matrix(jj, jj) = c_zero
      DO ii = jj + 1, SIZE(zz) 
        loewner_matrix(ii, jj) = (ff(ii) - ff(jj)) / (zz(ii) - zz(jj))
        loewner_matrix(jj, ii) = loewner_matrix(ii, jj) 
      END DO
    END DO 
    !
  END SUBROUTINE construct_Loewner_matrix


  SUBROUTINE extract_submatrix(matrix, mask, submatrix)
    !
    USE assert_module,   ONLY: assert
    !
    COMPLEX(dp), INTENT(IN) :: matrix(:,:)
    TYPE(mask_type), INTENT(IN) :: mask
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: submatrix(:,:)
    !
    CALL assert(ALL(SHAPE(mask%matrix) == SHAPE(matrix)), "matrix and mask must have consistent shape")
    !
    ALLOCATE(submatrix(mask%num_row, mask%num_col))
    submatrix = RESHAPE(PACK(matrix, mask%matrix), SHAPE(submatrix)) 
    !
  END SUBROUTINE extract_submatrix


  SUBROUTINE AND_mask_matrix(mask_row, mask_col, mask)
    !
    LOGICAL, INTENT(IN) :: mask_row(:), mask_col(:)
    TYPE(mask_type), INTENT(OUT) :: mask
    !
    INTEGER ii, jj
    !
    mask%num_row = COUNT(mask_row)
    mask%num_col = COUNT(mask_col)
    !
    ALLOCATE(mask%matrix(SIZE(mask_row), SIZE(mask_col)))
    FORALL(ii = 1:SIZE(mask_row), jj = 1:SIZE(mask_col))
      mask%matrix(ii, jj) = mask_row(ii).AND.mask_col(jj)
    END FORALL
    !
  END SUBROUTINE AND_mask_matrix


  SUBROUTINE find_pole(analytic_cont, pole, info)
    !
    USE lapack_module, ONLY: rational_vector, eigenvalue
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: pole(:)
    INTEGER, INTENT(OUT) :: info
    !
    COMPLEX(dp), ALLOCATABLE :: A_matrix(:,:), B_matrix(:,:)
    TYPE(rational_vector) lambda
    !
    !
    CALL setup_eigenvalue_problem(analytic_cont, A_matrix, B_matrix)
    CALL eigenvalue(A_matrix, B_matrix, lambda, info)
    IF (info /= no_error) RETURN
    CALL pole_from_noninfinite_eigenvalue(lambda, pole)
    !
  END SUBROUTINE find_pole


  SUBROUTINE setup_eigenvalue_problem(analytic_cont, A_matrix, B_matrix)
    !
    USE array_module,    ONLY: allocate_init_to
    USE constant_module, ONLY: c_zero, c_one
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: A_matrix(:,:), B_matrix(:,:)
    !
    INTEGER size_, ii
    !
    size_ = SIZE(analytic_cont%position) + 1
    CALL allocate_init_to([size_, size_], c_zero, A_matrix)
    CALL allocate_init_to([size_, size_], c_zero, B_matrix)
    DO ii = 2, size_
      A_matrix(1, ii) = c_one
      A_matrix(ii, ii) = analytic_cont%position(ii - 1)
      A_matrix(ii, 1) = analytic_cont%weight(ii - 1)
      B_matrix(ii, ii) = c_one
    END DO
    !
  END SUBROUTINE setup_eigenvalue_problem

  SUBROUTINE pole_from_noninfinite_eigenvalue(lambda, pole)
    !
    USE lapack_module, ONLY: rational_vector
    !
    TYPE(rational_vector), INTENT(IN) :: lambda
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: pole(:)
    !
    INTEGER num_pole, ii, ipole
    !
    num_pole = COUNT(eigenvalue_noninfinite(lambda%denominator))
    ALLOCATE(pole(num_pole))
    ipole = 1
    DO ii = 1, SIZE(lambda%denominator)
      IF (eigenvalue_noninfinite(lambda%denominator(ii))) THEN
        pole(ipole) = lambda%numerator(ii) / lambda%denominator(ii)
        ipole = ipole + 1
      END IF
    END DO
    !
  CONTAINS
    !
    ELEMENTAL LOGICAL FUNCTION eigenvalue_noninfinite(denominator)
      !
      USE constant_module, ONLY: eps14
      !
      COMPLEX(dp), INTENT(IN) :: denominator
      !
      eigenvalue_noninfinite = ABS(denominator) > eps14
      !
    END FUNCTION eigenvalue_noninfinite
    !
  END SUBROUTINE pole_from_noninfinite_eigenvalue


  SUBROUTINE calculate_residual(analytic_cont, pole, residual)
    !
    USE constant_module, ONLY: eps6, imag, c_zero
    !
    TYPE(aaa_approx), INTENT(IN) :: analytic_cont
    COMPLEX(dp), INTENT(IN) :: pole(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: residual(:)
    !
    COMPLEX(dp), PARAMETER :: shift(4) = [eps6 + c_zero, eps6 * imag, -eps6 + c_zero, -eps6 * imag]
    COMPLEX(dp), ALLOCATABLE :: zz_near_pole(:), ff_near_pole(:)
    !
    CALL generate_point_near_pole(pole, shift, zz_near_pole)
    CALL aaa_evaluate(analytic_cont, zz_near_pole, ff_near_pole)
    CALL average_residual(shift, ff_near_pole, residual)
    !
  END SUBROUTINE calculate_residual


  SUBROUTINE generate_point_near_pole(pole, shift, zz_near_pole)
    !
    COMPLEX(dp), INTENT(IN) :: pole(:), shift(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: zz_near_pole(:)
    INTEGER ii, first, last
    !
    last = 0
    ALLOCATE(zz_near_pole(SIZE(shift) * SIZE(pole)))
    DO ii = 1, SIZE(pole)
      first = last + 1
      last = ii * SIZE(shift)
      zz_near_pole(first:last) = pole(ii) + shift
    END DO
    !
  END SUBROUTINE generate_point_near_pole


  SUBROUTINE average_residual(shift, ff_near_pole, residual)
    !
    COMPLEX(dp), INTENT(IN) :: shift(:), ff_near_pole(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: residual(:)
    INTEGER ii, first, last
    !
    last = 0
    ALLOCATE(residual(SIZE(ff_near_pole) / SIZE(shift)))
    DO ii = 1, SIZE(residual)
      first = last + 1
      last = ii * SIZE(shift)
      residual(ii) = SUM(ff_near_pole(first:last) * shift) / SIZE(shift)
    END DO
    !
  END SUBROUTINE average_residual

END MODULE aaa_module

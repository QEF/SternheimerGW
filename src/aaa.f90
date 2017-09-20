!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2017
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
!> Provide the routines necessary to perform a analytic continuation using
!! the AAA method.
!!
!! Reference: Y. Nakatsukasa, O. Sete, L.N. Trefethen, https://arxiv.org/abs/1612.00337
MODULE aaa_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC aaa_coeff, aaa_eval

CONTAINS

!> The *core* AAA algorithm presented in Figure 4.1.
!!
!! Reimplement the AAA algorithm in Fortran. Because the original code uses
!! uppercase and lowercase letters for some variables (not differentiated in
!! Fortran), we use a *u* and *l* suffix to distinguish between them.
FUNCTION aaa_coeff(zu, fu, tol, mmax) RESULT (coeff)

  USE kinds,         ONLY: dp
  USE lapack_module, ONLY: svd
  USE norm_module,   ONLY: norm

  !> A vector of sample points.
  COMPLEX(dp), INTENT(IN) :: zu(:)

  !> A vector of data values at these points.
  COMPLEX(dp), INTENT(IN) :: fu(:)

  !> Relative tolerance, set to 1e-13 if omitted.
  REAL(dp), INTENT(IN), OPTIONAL :: tol

  !> Maximum degree of the polynomials, set to 100 if omitted.
  INTEGER,  INTENT(IN), OPTIONAL :: mmax

  !> The coefficients of the barycentric representation. The resulting array
  !! will have three components (in this order): The position \f$z_i\f$, the
  !! value at the position \f$f_i\f$, and the weight of the position \f$w_i\f$.
  !! The AAA approximation is then given by
  !! \f{equation}{
  !!   r(z) = \sum_i \frac{w_i f_i}{z - z_i} \middle/ \sum_i \frac{w_i}{z - z_i}
  !! \f}
  COMPLEX(dp), ALLOCATABLE :: coeff(:,:)

  !> Absolute tolerance used (input tolerance or default times function vector).
  REAL(dp) thres

  !> Maximum degree of polynomials.
  INTEGER mmax_

  !> next position used to improve the rational approximation
  INTEGER jj

  !> loop variable for polynomial degree
  INTEGER mm

  !> number of points in array
  INTEGER num_point

  !> index vector
  INTEGER,     ALLOCATABLE :: index_vector(:)

  !> vector of support points
  COMPLEX(dp), ALLOCATABLE :: zl(:)

  !> vector of function values
  COMPLEX(dp), ALLOCATABLE :: fl(:)

  !> vector of weights
  REAL(dp),    ALLOCATABLE :: wl(:)

  !> current estimate for rational approximation at support points
  COMPLEX(dp), ALLOCATABLE :: rational(:)

  !> Cauchy matrix
  COMPLEX(dp), ALLOCATABLE :: cauchy(:,:)

  !> Loewner matrix
  COMPLEX(dp), ALLOCATABLE :: loewner(:,:)

  !> singular values of decomposition
  REAL(dp),    ALLOCATABLE :: sigma(:)

  !> left matrix of SVD
  COMPLEX(dp), ALLOCATABLE :: umat(:,:)

  !> right matrix of SVD
  COMPLEX(dp), ALLOCATABLE :: vmat(:,:)

  !> weight vector
  COMPLEX(dp), ALLOCATABLE :: weight(:)

  !> coefficients of numerator
  COMPLEX(dp), ALLOCATABLE :: numerator(:)

  !> coefficients of denominator
  COMPLEX(dp), ALLOCATABLE :: denominator(:)

  !> error of approximation
  REAL(dp) err

  !
  ! use input values or set to default if input is not present
  !
  thres = norm(fu, infinity=.true.)
  IF (PRESENT(tol)) THEN
    thres = thres * tol
  ELSE
    thres = thres * 1e-13
  END IF
  !
  IF (PRESENT(mmax)) THEN
    mmax_ = mmax
  ELSE
    mmax_ = MIN(100, SIZE(zu))
  END IF

  !
  ! sanity check of input
  !
  num_point = SIZE(zu)
  IF (num_point == 0) THEN
    CALL errore(__FILE__, "cannot fit function to empty array", 1)
  END IF
  IF (SIZE(fu) /= num_point) THEN
    CALL errore(__FILE__, "size of sample points and function values inconsistent", 1)
  END IF
  IF (thres < 0.0_dp) THEN
    CALL errore(__FILE__, "negative tolerance not allowed", 1)
  END IF
  IF (mmax_ <= 0) THEN
    CALL errore(__FILE__, "degree of polynomial must be positive", 1)
  END IF

  !
  ! initialization
  !
  ALLOCATE(zl(num_point))
  ALLOCATE(fl(num_point))
  ALLOCATE(wl(num_point))
  ALLOCATE(index_vector(num_point))
  ALLOCATE(rational(num_point))
  ALLOCATE(cauchy(num_point, num_point))
  ALLOCATE(loewner(num_point, num_point))
  ALLOCATE(weight(num_point))
  ALLOCATE(numerator(num_point))
  ALLOCATE(denominator(num_point))
  !
  index_vector = [(mm, mm = 1, num_point)]
  !
  rational = SUM(fu) / num_point
  !
  cauchy = 0.0_dp

  !
  ! main loop
  !
  DO mm = 1, mmax_
    !
    ! select next support point (largest error to approximation)
    jj = MAXLOC(ABS(fu - rational), 1)
    !
    ! update support points, data values
    zl(mm) = zu(jj)
    fl(mm) = fu(jj)
    !
    ! update index vector
    CALL remove_element(jj, index_vector)
    !
    ! next column of Cauchy matrix
    where (zu /= zu(jj))
      cauchy(:, mm) = 1.0 / (zu - zu(jj))
    end where
    !
    ! Loewner matrix
    loewner(:, mm) = (fu - fl(mm)) * cauchy(:, mm)
    !
    ! SVD
    CALL svd(loewner(index_vector, :mm), sigma, umat, vmat)
    !
    ! weight vector = min sing vector
    ! note: SVD returns V^H
    weight(:mm) = CONJG(vmat(mm,:))
    !
    ! numerator and denominator
    numerator = MATMUL(cauchy(:, :mm), weight(:mm) * fl(:mm))
    denominator = MATMUL(cauchy(:, :mm), weight(:mm))
    !
    ! rational approximation
    rational = fu
    rational(index_vector) = numerator(index_vector) / denominator(index_vector)
    !
    ! stop if converged
    err = norm(fu - rational, infinity=.true.)
    IF (err <= thres) EXIT
    !
  END DO ! mm

  ! copy result in coeff array
  num_point = MIN(mm, mmax_)
  ALLOCATE(coeff(num_point, 3))
  !
  coeff(:, 1) = zl(:num_point)
  coeff(:, 2) = fl(:num_point)
  coeff(:, 3) = weight(:num_point)

END FUNCTION aaa_coeff

!> Evaluate the AAA approximation at a certain frequency.
!!
!! \f{equation}{
!!   r(z) = \sum_i \frac{w_i f_i}{z - z_i} \middle/ \frac{w_i}{z - z_i}
!! \f}
FUNCTION aaa_eval(coeff, freq) RESULT (res)

  USE constants, ONLY: eps12
  USE kinds,     ONLY: dp

  !> coefficients array - first dimension \f$z_i\f$, second dimension \f$f_i\f$,
  !!                      and third dimension \f$w_i\f$
  COMPLEX(dp), INTENT(IN) :: coeff(:, :)

  !> frequency at which the AAA approximation is evaluated
  COMPLEX(dp), INTENT(IN) :: freq

  !> value of the AAA approximation at the given frequency
  COMPLEX(dp) res

  !> loop variable
  INTEGER ii

  !> accumulate numerator
  COMPLEX(dp) num

  !> accumalate denominator
  COMPLEX(dp) den

  !> common term w_i / (z - z_i) in numerator and denominator
  COMPLEX(dp) frac

  ! initialize numerator and denominator
  num = 0.0_dp
  den = 0.0_dp

  ! sanity check - at least one frequency and three coefficients per frequency
  IF (SIZE(coeff, 1) == 0) &
    CALL errore(__FILE__, "at least one frequency required", 1)
  IF (SIZE(coeff, 2) /= 3) &
    CALL errore(__FILE__, "three coefficients (z_i, f_i, w_i) expected", 1)

  ! loop over coefficients
  DO ii = 1, SIZE(coeff, 1)

    ! skip points with 0 weight
    IF (ABS(coeff(ii, 3)) <= eps12) CONTINUE

    ! check for trivial case
    IF (ABS(coeff(ii, 1) - freq) <= eps12) THEN
      res = coeff(ii, 2)
      RETURN
    END IF

    ! evaluate common fraction w_i / (z - z_i)
    frac = coeff(ii, 3) / (freq - coeff(ii, 1))

    ! accumulate sums
    num = num + coeff(ii, 2) * frac
    den = den + frac

  END DO ! ii

  res = num / den

END FUNCTION aaa_eval

!> helper function to remove a certain element from an array
SUBROUTINE remove_element(element, array)

  !> value removed from the array
  INTEGER, INTENT(IN) :: element

  !> *on input* array with all elements <br>
  !! *on output* array with selected elements removed
  INTEGER, INTENT(INOUT), ALLOCATABLE :: array(:)

  !> mask indicating which elements to remove
  LOGICAL, ALLOCATABLE :: mask(:)

  !> number of elements after removal
  INTEGER num_elem

  !> temporary array to copy the data
  INTEGER, ALLOCATABLE :: copy_array(:)

  ! trivial case - array empty or not allocated
  IF (.NOT.ALLOCATED(array)) THEN
    ALLOCATE(array(0))
    RETURN
  END IF
  IF (SIZE(array) == 0) THEN
    RETURN
  END IF

  ! find matching elements
  ALLOCATE(mask(SIZE(array)))
  mask = (array /= element)
  num_elem = COUNT(mask)

  ! trivial case - no matching elements
  IF (num_elem == SIZE(array)) THEN
    RETURN
  END IF

  ! create smaller array
  ALLOCATE(copy_array(num_elem))

  ! remove elements from array
  copy_array = PACK(array, mask)

  ! move allocation to output array
  CALL MOVE_ALLOC(copy_array, array)

END SUBROUTINE remove_element

END MODULE aaa_module

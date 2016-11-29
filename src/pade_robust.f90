!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file have been taken from the chebfun code.
! See http://www.chebfun.org/ for Chebfun information.
! 
! Copyright (C) 2010 - 2016 
! The University of Oxford and The Chebfun Developers,
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
!> Provides the routines to evaluate the Pade approximation to a function.
MODULE pade_module

  IMPLICIT NONE

  PRIVATE

  PUBLIC pade_robust

  !> LAPACK flag to determine work size
  INTEGER, PARAMETER :: determine = -1

CONTAINS

  !> Pade approximation to a function.
  !!
  !! Constructs a Pade approximant to a function using the robust algorithm from
  !! [1] based on the SVD.
  !!
  !! This code is included in the Chebfun distribution for the convenience of
  !! readers of _Approximation Theory and Approximation Practice_, but it is not
  !! actually a Chebfun code. A Chebfun analogue is CHEBPADE.
  !!
  !! <h4>References:</h4>
  !! [1] P. Gonnet, S. Guettel, and L. N. Trefethen, "ROBUST PADE APPROXIMATION 
  !!     VIA SVD", SIAM Rev., 55:101-117, 2013.
  !!
  SUBROUTINE pade_robust(radius, func, deg_num, deg_den, coeff_num, coeff_den, tol_coeff, tol_fft)

    USE constants,   ONLY: eps14
    USE kinds,       ONLY: dp
    USE norm_module, ONLY: norm

    !> The radius of the circle in the complex plane.
    REAL(dp),    INTENT(IN) :: radius

    !> The values of the function evaluated on a circle in the complex plane.
    !! The derivatives of the functions are computed via FFT.
    COMPLEX(dp), INTENT(IN) :: func(:)

    !> The degree of the numerator (must be positive), potentially modified if
    !! Pade approximation has vanishing terms
    INTEGER,     INTENT(INOUT) :: deg_num

    !> The degree of the denominator (must be positive), potentially modified if
    !! Pade approximation has vanishing terms
    INTEGER,     INTENT(INOUT) :: deg_den

    !> The Pade coefficient vector of the numerator
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: coeff_num(:)

    !> The Pade coefficient vector of the denominator
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: coeff_den(:)

    !> The optional **tol_coeff** argument specifies the relative tolerance for
    !! the coefficients; if omitted, it defaults to 1e-14. Set to 0 to turn off
    !! robustness.
    REAL(dp),    OPTIONAL,    INTENT(IN)  :: tol_coeff

    !> The optional **tol_fft** argument specifies the relative tolerance for
    !! the FFT. Coefficients smaller than it are set to 0; if omitted it defaults
    !! to the value of **tol_coeff**.
    REAL(dp),    OPTIONAL,    INTENT(IN)  :: tol_fft

    !> number of converged elements
    INTEGER rho

    !> number of leading or trailing zeros
    INTEGER lam

    !> local variable for the relative tolerance
    REAL(dp) rel_tol

    !> absolute value for the tolerance
    REAL(dp) abs_tol

    !> local variable for the relative tolerance of the fft
    REAL(dp) rel_tol_fft

    !> singular values of matrix
    REAL(dp),    ALLOCATABLE :: sigma(:)

    !> the coefficients of the Taylor series
    COMPLEX(dp), ALLOCATABLE :: coeff(:)

    !> the first row of the Toeplitz matrix
    COMPLEX(dp), ALLOCATABLE :: row(:)

    !> the first column of the Toeplitz matrix
    COMPLEX(dp), ALLOCATABLE :: col(:)

    !> the constructed Toeplitz matrix
    COMPLEX(dp), ALLOCATABLE :: zmat(:,:)

    !> matrix for calculating the numerical rank
    COMPLEX(dp), ALLOCATABLE :: cmat(:,:)

    !> U matrix for singular value decomposition
    COMPLEX(dp), ALLOCATABLE :: umat(:,:)

    !> V matrix for singular value decomposition
    COMPLEX(dp), ALLOCATABLE :: vmat(:,:)

    !> D matrix is diagonal, so we can store it in a vector
    COMPLEX(dp), ALLOCATABLE :: dmat(:)

    !> Q matrix of the QR factorization
    COMPLEX(dp), ALLOCATABLE :: qmat(:,:)

    !> numeric precision (1 + eps) > 1
    REAL(dp),    PARAMETER :: eps = EPSILON(1.0_dp)

    !> complex constant of 1
    COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

    !> complex constant of 0
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    !> use infinity norm
    LOGICAL, PARAMETER :: inf = .TRUE.

    ! sanity check of the input
    IF (radius <= 0) &
      CALL errore(__FILE__, "radius in the complex plane must be > 0", 1)
    IF (deg_num < 0) &
      CALL errore(__FILE__, "degree of numerator must be positive", deg_num)
    IF (deg_den < 0) &
      CALL errore(__FILE__, "degree of denominator must be positive", deg_den)

    ! default value of tolerance 1e-14
    IF (PRESENT(tol_coeff)) THEN
      rel_tol = tol_coeff
    ELSE
      rel_tol = eps14
    END IF

    ! use same tolererance for FFT and coeff unless specified otherwise
    IF (PRESENT(tol_fft)) THEN
      rel_tol_fft = tol_fft
    ELSE
      rel_tol_fft = rel_tol
    END IF

    ! Compute coefficients
    CALL pade_derivative(radius, func, rel_tol_fft, deg_num + deg_den + 1, coeff)

    ! determine the absolute value of the tolerance
    abs_tol = rel_tol * norm(coeff)

    !
    ! Compute the Pade approximation.
    !
    IF (norm(coeff(1:deg_num + 1), inf) <= rel_tol * norm(coeff, inf)) THEN
      !
      ! special case - the function is 0
      !
      ! set numerator to 0
      deg_num = 0
      ALLOCATE(coeff_num(deg_num + 1))
      coeff_num = zero
      !
      ! set denominator to 1
      deg_den = 0
      ALLOCATE(coeff_den(deg_den + 1))
      coeff_den = one
      !
      RETURN
      !
    END IF
    !
    ! general case
    !
    ! First row/column of Toeplitz matrix.
    ALLOCATE(row(deg_den + 1))
    row = zero
    row(1) = coeff(1)
    !
    ALLOCATE(col(SIZE(coeff)))
    col = coeff;
    !
    ! Do diagonal hopping across block. 
    DO WHILE (.TRUE.)
      !
      ! Special case n == 0.
      IF (deg_den == 0) THEN
        !
        ALLOCATE(coeff_num(deg_num + 1))
        coeff_num = coeff(:deg_num + 1) 
        !
        ALLOCATE(coeff_den(deg_den + 1))
        coeff_den = one
        !
        EXIT
        !
      END IF
      !
      ! Form Toeplitz matrix.
      CALL toeplitz_nonsym(col(1:deg_num + deg_den + 1), row(1:deg_den + 1), zmat)
      !
      ! Compute numerical rank.
      ALLOCATE(cmat(deg_den, SIZE(zmat, 2))) 
      cmat = zmat(deg_num + 2:deg_num + deg_den + 1, :)
      CALL svd(cmat, sigma)
      rho = COUNT(sigma > abs_tol)
      !
      IF (rho == deg_den) EXIT
      !
      ! Decrease mn, n if rank-deficient.
      deg_num = deg_num - (deg_den - rho)
      deg_den = rho
      !
      DEALLOCATE(zmat, cmat)
      !
    END DO ! while
    !
    DEALLOCATE(coeff)
    !
    ! Hopping finished. Now compute b and a.
    IF (deg_num > 1) THEN
      !
      CALL svd(cmat, sigma, umat, vmat)
      !
      ! Null vector gives b.
      ALLOCATE(coeff_den(deg_den + 1))
      coeff_den = vmat(deg_den + 1, :)
      !
      ! Do final computation via reweighted QR for better zero preservation.
      ALLOCATE(dmat(SIZE(coeff_den)))
      dmat = ABS(coeff_den) + SQRT(eps)
      ! replace C -> (C D)^T
      CALL matmul_transpose(cmat, dmat)
      CALL qr(cmat, qmat)
      !
      ! Compensate for reweighting.
      coeff_den = dmat * qmat(:, deg_den + 1)
      coeff_den = coeff_den / norm(coeff_den)
      !
      ! Multiplying gives a.
      ALLOCATE(coeff_num(deg_num + 1))
      coeff_num = MATMUL(zmat(1:deg_num + 1, 1:deg_den + 1), coeff_den)
      !
      ! Count leading zeros of b.
      lam = first_nonzero(coeff_den, rel_tol) - 1
      !
      IF (lam > 0) THEN
        !
        ! Discard leading zeros of b and a.
        deg_num = deg_num - lam
        deg_den = deg_den - lam
        !
        ALLOCATE(coeff(deg_num + 1))
        coeff = coeff_num(lam + 1:)
        CALL MOVE_ALLOC(coeff, coeff_num)
        !
        ALLOCATE(coeff(deg_den + 1))
        coeff = coeff_den(lam + 1:)
        CALL MOVE_ALLOC(coeff, coeff_den)
        !
      END IF ! lam > 0
      !
      ! Discard trailing zeros of b.
      lam = last_nonzero(coeff_den, rel_tol) - 1
      !
      IF (lam /= deg_den) THEN
        !
        deg_den = lam
        !
        ALLOCATE(coeff(deg_den + 1))
        coeff = coeff_den(:deg_den + 1)
        CALL MOVE_ALLOC(coeff, coeff_den)
        !
      END IF ! lam /= deg_den
      !
    END IF ! deg_num > 1
    !
    ! Discard trailing zero coefficients in a.
    lam = last_nonzero(coeff_num, abs_tol) - 1
    !
    IF (lam /= deg_num) THEN
      !
      deg_num = lam
      !
      ALLOCATE(coeff(deg_num + 1))
      coeff = coeff_num(:deg_num + 1)
      CALL MOVE_ALLOC(coeff, coeff_num)
      !
    END IF ! lam /= deg_num
    !
    ! Normalize.
    coeff_num = coeff_num / coeff_den(1)
    coeff_den = coeff_den / coeff_den(1)

  END SUBROUTINE pade_robust

  !> determine the derivatives of the function
  !!
  !! we use a FFT to evaluate the derivative of the function
  SUBROUTINE pade_derivative(radius, func, tol, num_deriv, deriv)

    USE constants,   ONLY: eps14
    USE fft_scalar,  ONLY: cft_1z
    USE kinds,       ONLY: dp
    USE norm_module, ONLY: norm

    !> The radius of the circle in the complex plane.
    REAL(dp),    INTENT(IN) :: radius

    !> The values of the function evaluated on a circle in the complex plane.
    COMPLEX(dp), INTENT(IN) :: func(:)

    !> tolerance to identify small coeffficients
    REAL(dp),    INTENT(IN) :: tol

    !> The number of derivatives that should be generated
    INTEGER,     INTENT(IN) :: num_deriv

    !> The derivatives of the functions are computed via FFT.
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: deriv(:)

    !> use a backward FFT
    INTEGER,     PARAMETER :: backward = -1

    !> the number of FFTs done per call
    INTEGER,     PARAMETER :: num_fft = 1

    !> real constant of 1
    REAL(dp),    PARAMETER :: one = 1.0_dp

    !> complex constant of 0
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    !> the number of points in a FFT
    INTEGER num_point

    !> counter on the derivatives
    INTEGER ipoint

    !> rescale the derivatives if the radius is not 1.0
    REAL(dp) rescale

    !> absolute value of the tolerance
    REAL(dp) abs_tol

    !> work array for FFT
    COMPLEX(dp), ALLOCATABLE :: work(:)

    ! create array for FFT
    num_point = SIZE(func)
    ALLOCATE(work(num_point))

    ! evalute FFT of function
    ! work contains now the derivatives up to a factor
    CALL cft_1z(func, num_fft, num_point, num_point, backward, work)

    ! create array for the derivatives
    ALLOCATE(deriv(num_deriv))
    deriv = zero

    ! evaluate the derivatives (truncating or filling with zeros as needed)
    num_point = MIN(num_deriv, num_point)
    deriv(:num_point) = work(:num_point)

    ! rescale the derivatives by radius^(-order of derivative)
    IF (ABS(radius - one) > eps14) THEN
      !
      rescale = one
      DO ipoint = 2, num_point
        !
        rescale = rescale / radius
        deriv(ipoint) = deriv(ipoint) * rescale
        !
      END DO ! ipoint
    END IF ! radius /= 1

    ! Discard near-zero coefficients.
    abs_tol = tol * norm(deriv)
    WHERE(ABS(deriv) < abs_tol)
      deriv = zero
    END WHERE

    ! Remove imaginary rounding errors. (Make real functions real.)
    IF (norm(AIMAG(deriv), .TRUE.) < abs_tol) THEN
      deriv = REAL(deriv)
    END IF

  END SUBROUTINE pade_derivative

  !> create a nonsymetric Toeplitz matrix (mathlab-like behavior)
  SUBROUTINE toeplitz_nonsym(col, row, matrix)

    USE constants, ONLY: eps14
    USE kinds,     ONLY: dp

    !> the first column of the matrix
    COMPLEX(dp), INTENT(IN) :: col(:)

    !> the first row of the matrix
    COMPLEX(dp), INTENT(IN) :: row(:)

    !> the resulting Toeplitz matrix
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: matrix(:,:)

    !> the number of row and colums
    INTEGER num_row, num_col

    !> counter on row and colums
    INTEGER irow, icol

    num_row = SIZE(col)
    num_col = SIZE(row)
    ALLOCATE(matrix(num_row, num_col))

    ! trivial case - zero length array
    IF (num_row == 0 .OR. num_col == 0) RETURN
    
    ! sanity check of the input
    IF (ABS(col(1) - row(1)) > eps14) THEN
      WRITE(0,*) 'Warning: First element of input column does not match first &
                 &element of input row. Column wins diagonal conflict.'
    END IF

    ! create the Toeplitz matrix
    DO icol = 1, num_col
      DO irow = 1, num_row
        !
        IF (irow < icol) THEN
          ! use row for upper triangle
          matrix(irow, icol) = row(icol - irow + 1)
        ELSE
          ! use col for lower triangle and diagonal
          matrix(irow, icol) = col(irow - icol + 1)
        END IF
        !
      END DO ! irow
    END DO ! icol

  END SUBROUTINE toeplitz_nonsym

  !> wrapper for LAPACK singular value decomposition routine
  SUBROUTINE svd(matrix, sigma, umat, vmat)

    USE kinds, ONLY: dp

    !> The matrix for which the SVD \f$A = U \Sigma V^{\text{H}}\f$ is evaluated.
    COMPLEX(dp), INTENT(IN) :: matrix(:,:)

    !> singular values of the matrix \f$\Sigma\f$ (ascending)
    REAL(dp),    ALLOCATABLE, INTENT(OUT) :: sigma(:)

    !> Unitary matrix U (left)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL :: umat(:,:)

    !> Unitary matrix V (right), note returns \f$V^{\text{H}}\f$.
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL :: vmat(:,:)

    !> jobz parameter for LAPACK SVD
    CHARACTER(1) jobz

    !> number of rows M of the input matrix
    INTEGER num_row

    !> number of columns N of the input matrix
    INTEGER num_col

    !> miniumum of number of rows and number of columns
    INTEGER num_min

    !> number of elements in work array
    INTEGER num_work

    !> error flag returned by LAPACK
    INTEGER ierr

    !> optimal size of work
    COMPLEX(dp) opt_size

    !> integer work array for SVD
    INTEGER,     ALLOCATABLE :: iwork(:)

    !> real work array for SVD
    REAL(dp),    ALLOCATABLE :: rwork(:)

    !> copy of the input matrix, will be destroyed or overwritten by LAPACK call
    COMPLEX(dp), ALLOCATABLE :: amat(:,:)

    !> work array for SVD
    COMPLEX(dp), ALLOCATABLE :: work(:)

    !> complex constant of 0
    COMPLEX(dp), PARAMETER   :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    ! set helper variables
    num_row = SIZE(matrix, 1)
    num_col = SIZE(matrix, 2)
    num_min = MIN(num_row, num_col)

    ! sanity check - either both or none of U and V present
    IF ((PRESENT(umat).AND..NOT.PRESENT(vmat)).OR. &
        (PRESENT(vmat).AND..NOT.PRESENT(umat))) THEN
      CALL errore(__FILE__, "either both or none of the optional arguments must be present", 1)
    END IF

    ! allocate arrays for output
    ! U is M x M matrix
    IF (PRESENT(umat)) ALLOCATE(umat(num_row, num_row))
    ! V is N x N matrix
    IF (PRESENT(vmat)) ALLOCATE(vmat(num_col, num_col))
    ! Sigma has MIN(N, M) diagonal entries
    ALLOCATE(sigma(num_min))
    ! integer work array
    ALLOCATE(iwork(8 * num_min))
    ! real work array
    IF (PRESENT(umat)) THEN
      ALLOCATE(rwork(5 * num_min**2 + 7 * num_min))
    ELSE
      ! note: newer version of LAPACK require only 5 * num_min
      ALLOCATE(rwork(7 * num_min))
    END IF

    ! create copy of input matrix, because LAPACK destroys input
    ALLOCATE(amat(num_row, num_col))
    CALL ZCOPY(SIZE(amat), matrix, 1, amat, 1)

    IF (.NOT.PRESENT(umat)) THEN
      ! we evaluate only the singular values
      jobz = 'N'
    ELSE IF (num_row > num_col) THEN
      ! if the number of rows is larger than the number of columns, 
      ! we only evaluate the first N columns of U
      jobz = 'O'
    ELSE
      ! we evaluate all elements.
      jobz = 'A'
    END IF

    ! determine optimum work array size
    CALL ZGESDD(jobz, num_row, num_col, amat, num_row, sigma, umat, num_row, &
                vmat, num_col, opt_size, determine, rwork, iwork, ierr)
    CALL errore(__FILE__, "error calculating work size for SVD", ierr)

    ! create work array
    num_work = NINT(ABS(opt_size))
    ALLOCATE(work(num_work))

    ! perform SVD
    CALL ZGESDD(jobz, num_row, num_col, amat, num_row, sigma, umat, num_row, &
                vmat, num_col, work, num_work, rwork, iwork, ierr)
    CALL errore(__FILE__, "error calculating SVD", ierr)

    ! if jobz is 'O' output was written to A instead of U
    IF (jobz == 'O') THEN
      umat(:, num_col + 1:) = zero
      CALL ZCOPY(SIZE(amat), amat, 1, umat, 1)
    END IF

  END SUBROUTINE svd
 
  !> wrapper around the LAPACK QR factorization routines
  SUBROUTINE qr(amat, qmat)

    USE kinds, ONLY: dp

    !> the matrix for which a QR factorization is constructed
    COMPLEX(dp), INTENT(IN) :: amat(:,:)

    !> the matrix Q of the factorization
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: qmat(:,:)

    !> the number of rows in the matrix A
    INTEGER num_row

    !> the number of columns in the matrix A
    INTEGER num_col

    !> dimensionalty of the work array
    INTEGER num_work

    !> error flag raised by LAPACK routines
    INTEGER ierr

    !> optimal size of the work array
    COMPLEX(dp) opt_size

    !> array for the scalar reflectors
    COMPLEX(dp), ALLOCATABLE :: tau(:)

    !> work array for LAPACK
    COMPLEX(dp), ALLOCATABLE :: work(:)

    ! set helper variable
    num_row = SIZE(amat, 1)
    num_col = SIZE(amat, 2)

    ! for this implementation we assume that num_row >= num_col
    IF (num_row < num_col) &
      CALL errore(__FILE__, "the case num_row < num_col is not implemented", 1)
    
    ALLOCATE(tau(MIN(num_row, num_col)))

    ! copy A to Q because LAPACK will overwrite the input
    ALLOCATE(qmat(num_row, num_row))
    CALL ZCOPY(SIZE(amat), amat, 1, qmat, 1)

    ! determine size for work array
    CALL ZGEQRF(num_row, num_col, qmat, num_row, tau, opt_size, determine, ierr)
    CALL errore(__FILE__, "failed to determine size for work in QR factorization", ierr)

    ! create work array
    num_work = NINT(ABS(opt_size))
    ALLOCATE(work(num_work))

    ! perform the QR factorization
    CALL ZGEQRF(num_row, num_col, qmat, num_row, tau, work, num_work, ierr)
    CALL errore(__FILE__, "error in QR factorization", ierr)

    ! check if we need a larger work array
    CALL ZUNGQR(num_row, num_row, SIZE(tau), qmat, num_row, tau, opt_size, determine, ierr)
    CALL errore(__FILE__, "failed to determine size of work array", ierr)

    ! create larger work array if necessary
    IF (NINT(ABS(opt_size)) > num_work) THEN
      DEALLOCATE(work)
      num_work = NINT(ABS(opt_size))
      ALLOCATE(work(num_work))
    END IF

    ! now generate the Q matrix
    CALL ZUNGQR(num_row, num_row, SIZE(tau), qmat, num_row, tau, work, num_work, ierr)
    CALL errore(__FILE__, "error constructing Q matrix", ierr)

  END SUBROUTINE qr

  !> Evaluate \f$(C D)^{\text{T}}\f$, where D is diagonal.
  SUBROUTINE matmul_transpose(cmat, dmat)

    USE kinds, ONLY: dp

    !> *on input* the matrix C <br>
    !! *on output* the result \f$(C D)^{\text{T}}\f$
    COMPLEX(dp), INTENT(INOUT), ALLOCATABLE :: cmat(:,:)

    !> the diagonal elements of the matrix D
    COMPLEX(dp), INTENT(IN)  :: dmat(:)

    !> work array that will store the output
    COMPLEX(dp), ALLOCATABLE :: work(:,:)

    ! counter on the columns of C
    INTEGER icol

    ! create work array - transpose of C
    ALLOCATE(work(SIZE(cmat, 2), SIZE(cmat, 1)))

    DO icol = 1, SIZE(cmat, 2)
      !
      ! C_ij D_jj
      work(icol, :) = cmat(:, icol) * dmat(icol)
      !
    END DO ! icol

    ! copy work to output
    CALL MOVE_ALLOC(work, cmat)

  END SUBROUTINE matmul_transpose

  !> determine first nonzero element of complex array (up to a tolerance)
  FUNCTION first_nonzero(array, tol)

    USE kinds, ONLY: dp

    !> the array of which the nonzero elements are determined
    COMPLEX(dp), INTENT(IN) :: array(:)

    !> the tolerance up to which elements are considered as 0
    REAL(dp),    INTENT(IN) :: tol

    !> the position of the first nonzero element
    INTEGER first_nonzero

    !> a ordered array of integers
    INTEGER, ALLOCATABLE :: order(:)

    !> counter to initialize the order array
    INTEGER ii

    ! create a ordered integer array
    ALLOCATE(order(SIZE(array)))
    order = [ (ii, ii = 1, SIZE(array)) ]

    ! determine the first element for which the mask is true
    first_nonzero = MINLOC(order, 1, ABS(array) > tol)

  END FUNCTION first_nonzero

  !> determine last nonzero element of complex array (up to a tolerance)
  FUNCTION last_nonzero(array, tol)

    USE kinds, ONLY: dp

    !> the array of which the nonzero elements are determined
    COMPLEX(dp), INTENT(IN) :: array(:)

    !> the tolerance up to which elements are considered as 0
    REAL(dp),    INTENT(IN) :: tol

    !> the position of the last nonzero element
    INTEGER last_nonzero

    !> a ordered array of integers
    INTEGER, ALLOCATABLE :: order(:)

    !> counter to initialize the order array
    INTEGER ii

    ! create a ordered integer array
    ALLOCATE(order(SIZE(array)))
    order = [ (ii, ii = 1, SIZE(array)) ]

    ! determine the last element for which the mask is true
    last_nonzero = MAXLOC(order, 1, ABS(array) > tol)

  END FUNCTION last_nonzero

END MODULE pade_module

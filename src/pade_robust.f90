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
  SUBROUTINE pade_robust(radius, func, deg_num, deg_den, coeff_num, coeff_den, tol_in)

    USE constants, ONLY: eps14
    USE kinds,     ONLY: dp

    !> The radius of the circle in the complex plane.
    REAL(dp),    INTENT(IN) :: radius

    !> The values of the function evaluated on a circle in the complex plane.
    !! The derivatives of the functions are computed via FFT.
    COMPLEX(dp), INTENT(IN) :: func(:)

    !> The degree of the numerator (must be positive)
    INTEGER,     INTENT(IN) :: deg_num

    !> The degree of the denominator (must be positive)
    INTEGER,     INTENT(IN) :: deg_den

    !> The Pade coefficient vector of the numerator
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: coeff_num(:)

    !> The Pade coefficient vector of the denominator
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: coeff_den(:)

    !> The optional **tol** argument specifies the relative tolerance; if
    !! omitted, it defaults to 1e-14. Set to 0 to turn off robustness.
    REAL(dp),    OPTIONAL,    INTENT(IN)  :: tol_in

    !> local variable for the tolerance
    REAL(dp) tol

    ! default value of tolerance 1e-14
    IF (PRESENT(tol_in)) THEN
      tol = tol_in
    ELSE
      tol = eps14
    END IF

    ! sanity check of the input
    IF (radius <= 0) &
      CALL errore(__FILE__, "radius in the complex plane must be > 0", 1)
    IF (deg_num < 0) &
      CALL errore(__FILE__, "degree of numerator must be positive", deg_num)
    IF (deg_den < 0) &
      CALL errore(__FILE__, "degree of denominator must be positive", deg_den)

  END SUBROUTINE pade_robust

  !> determine the derivatives of the function
  !!
  !! we use a FFT to evaluate the derivative of the function
  SUBROUTINE pade_derivative(radius, func, num_deriv, deriv)

    USE constants,  ONLY: eps14
    USE fft_scalar, ONLY: cft_1z
    USE kinds,      ONLY: dp

    !> The radius of the circle in the complex plane.
    REAL(dp),    INTENT(IN) :: radius

    !> The values of the function evaluated on a circle in the complex plane.
    COMPLEX(dp), INTENT(IN) :: func(:)

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

    num_row = SIZE(row)
    num_col = SIZE(col)
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
        IF (irow > icol) THEN
          ! use row for upper triangle
          matrix(irow, icol) = row(irow - icol + 1)
        ELSE
          ! use col for lower triangle and diagonal
          matrix(irow, icol) = col(icol - irow + 1)
        END IF
        !
      END DO ! irow
    END DO ! icol

  END SUBROUTINE toeplitz_nonsym

END MODULE pade_module

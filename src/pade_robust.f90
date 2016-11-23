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
    IF (deg_num < 0) &
      CALL errore(__FILE__, "degree of numerator must be positive", deg_num)
    IF (deg_den < 0) &
      CALL errore(__FILE__, "degree of denominator must be positive", deg_den)

  END SUBROUTINE pade_robust

END MODULE pade_module

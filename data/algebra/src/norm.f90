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
!> This module provides an interface to the LAPACK norm routines.
MODULE norm_module

  IMPLICIT NONE

  !> evaluate the norm of a real or complex vector
  INTERFACE norm
    MODULE PROCEDURE norm_real, norm_complex
  END INTERFACE norm

CONTAINS

  !> determine norm of a real vector
  FUNCTION norm_real(vector, infinity)

    USE kinds, ONLY: dp

    !> the vector of which the norm is determined
    REAL(dp), INTENT(IN) :: vector(:)

    !> optionally the infinity norm is calculated (default 2-norm)
    LOGICAL,  INTENT(IN), OPTIONAL :: infinity

    !> the norm of the vector
    REAL(dp) norm_real

    !> character that identifies which norm is used
    CHARACTER(1) cnorm

    !> work for infinity norm
    REAL(dp) work

    !> declare interface to LAPACK's DLANGE
    REAL(dp), EXTERNAL :: DLANGE

    ! default to 2-norm
    cnorm = 'F'

    ! if optional flag is present and set switch to infinity norm
    IF (PRESENT(infinity)) THEN
      IF (infinity) cnorm = '1'
    END IF

    ! evaluate infinity norm of the vector
    ! 1 = vector, i.e. one row
    norm_real = DLANGE(cnorm, 1, SIZE(vector), vector, 1, work)

  END FUNCTION norm_real
  
  !> determine norm of a complex vector
  FUNCTION norm_complex(vector, infinity)

    USE kinds, ONLY: dp

    !> the vector of which the norm is determined
    COMPLEX(dp), INTENT(IN) :: vector(:)

    !> optionally the infinity norm is calculated (default 2-norm)
    LOGICAL,     INTENT(IN), OPTIONAL :: infinity

    !> the norm of the vector
    REAL(dp) norm_complex

    !> character that identifies which norm is used
    CHARACTER(1) cnorm

    !> work for infinity norm
    REAL(dp) work

    !> declare interface to LAPACK's ZLANGE
    REAL(dp), EXTERNAL :: ZLANGE

    ! default to 2-norm
    cnorm = 'F'

    ! if optional flag is present and set switch to infinity norm
    IF (PRESENT(infinity)) THEN
      IF (infinity) cnorm = '1'
    END IF

    ! evaluate infinity norm of the vector
    ! 1 = vector, i.e. one row
    norm_complex = ZLANGE(cnorm, 1, SIZE(vector), vector, 1, work)

  END FUNCTION norm_complex

END MODULE norm_module

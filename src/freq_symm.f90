!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2017
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
!> Provides helper routines to extend a reduced mesh of frequencies to a denser
!! set of points.
MODULE freq_symm_module

  IMPLICIT NONE

CONTAINS

  !> Expand an array from the symmetry reduced to the full frequency mesh.
  !!
  !! We assume an array with the following property
  !! \f{equation}{
  !!   A(\omega) = A^\ast(-\omega)~.
  !! \f}
  !! Then, we can extend an given input array for which only the first half
  !! of the frequencies was calculated to the full mesh by adding the complex
  !! conjugate to the end.
  SUBROUTINE freq_symm(freq, array)

    USE kinds, ONLY: dp

    !> *on input*: the first half of the frequency array
    !! *on output*: the full frequency array
    COMPLEX(dp), INTENT(INOUT) :: freq(:)

    !> *on input*: The first half of the array \f$A\f$ <br>
    !! *on output*: The full array \f$A\f$ (extended by its complex conjugate).
    COMPLEX(dp), INTENT(INOUT) :: array(:)

    ! mid point in the array
    INTEGER middle

    !
    ! sanity test
    !
    ! frequency and array should habe same dimension
    IF (SIZE(freq) /= SIZE(array)) THEN
      CALL errore(__FILE__, "array and frequency mesh inconsistent", 1)
    END IF
    ! total number of frequencies should be even
    IF (MOD(SIZE(array), 2) /= 0) THEN
      CALL errore(__FILE__, "expected an array with even number of frequencies", SIZE(array))
    END IF

    !
    ! copy first half to second half and complex conjugate
    !
    middle = SIZE(array) / 2
    freq(middle + 1:) = -freq(:middle)
    array(middle + 1:) = CONJG(array(:middle))

  END SUBROUTINE freq_symm

END MODULE freq_symm_module

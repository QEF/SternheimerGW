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
  SUBROUTINE freq_symm(freq_in, freq_out, array)

    USE freqbins_module, ONLY: freqbins_type
    USE kinds,           ONLY: dp

    !> the definition of the input frequency mesh
    TYPE(freqbins_type), INTENT(IN)         :: freq_in

    !> the full frequency array
    COMPLEX(dp), INTENT(INOUT), ALLOCATABLE :: freq_out(:)

    !> *on input*: The first half of the array \f$A\f$ <br>
    !! *on output*: The full array \f$A\f$ (extended by its complex conjugate).
    COMPLEX(dp), INTENT(INOUT), OPTIONAL    :: array(:)

    ! mid point in the array
    INTEGER middle

    ! create frequency mesh
    IF (.NOT.ALLOCATED(freq_out)) ALLOCATE(freq_out(freq_in%num_freq()))

    ! trivial case - symmetry not used => return same frequency used for solver
    IF (.NOT. freq_in%use_symmetry) THEN
      freq_out = freq_in%solver
      RETURN
    END IF

    !
    ! sanity test
    !
    IF (PRESENT(array)) THEN

      ! frequency and array should habe same dimension
      IF (SIZE(freq_out) /= SIZE(array)) THEN
        CALL errore(__FILE__, "array and frequency mesh inconsistent", 1)
      END IF

    END IF

    !
    ! copy first half to second half and complex conjugate
    !
    middle = SIZE(freq_in%solver)
    freq_out(:middle)     =  freq_in%solver
    freq_out(middle + 1:) = -freq_in%solver
    IF (PRESENT(array)) array(middle + 1:) = CONJG(array(:middle))

  END SUBROUTINE freq_symm

END MODULE freq_symm_module

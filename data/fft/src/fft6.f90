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
!> Provide routines for 6 dimensional FFT.
!!
!! This module allows to do a Fourier transform of quantities that have to
!! spacial coordinates into reciprocal space and the reverse
!! \f{equation}{
!!   f(r, r') \underset{\text{fwfft6}}{\longrightarrow} f(G, G')
!!            \underset{\text{invfft6}}{\longrightarrow} f(r, r')~.
!! \f}
MODULE fft6_module

  IMPLICIT NONE

  !> interface for the fwfft6 routines
  INTERFACE fwfft6
    MODULE PROCEDURE fwfft6_same, fwfft6_diff
  END INTERFACE fwfft6

  !> interface for the invfft6 routines
  INTERFACE invfft6
    MODULE PROCEDURE invfft6_same, invfft6_diff
  END INTERFACE invfft6

CONTAINS

  !> Transform an input array \f$f(r, r')\f$ from real to reciprocal space
  !! \f$f(G, G')\f$.
  !!
  !! This calls internally the fftw6_diff routine.
  SUBROUTINE fwfft6_same(grid_type, f, dfft, omega)

    USE kinds,          ONLY: dp
    USE fft_interfaces, ONLY: fwfft
    USE fft_types,      ONLY: fft_type_descriptor

    !> grid type used for the Fourier transform, this is passed to fwfft;
    !! check the definition of fwfft for a list of options
    CHARACTER(*), INTENT(IN)    :: grid_type

    !> *on input*  the array in real space \f$f(r, r')\f$ <br>
    !! *on output* the array in reciprocal space \f$f(G, G')\f$
    COMPLEX(dp),  INTENT(INOUT) :: f(:,:)

    !> FFT descriptor - this defines the way the Fourier transform is
    !! executed, must be consistent with grid_type (used for both FFT)
    !! @note the code does not check explicitly if dfft and grid_type are compatible
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft

    !> volume of the unit cell
    REAL(dp),     INTENT(IN)    :: omega

    ! do the FFT via the routine which allows for different arguments
    CALL fwfft6_diff(grid_type, f, dfft, dfft, omega)

  END SUBROUTINE fwfft6_same

  !> Transform an input array \f$f(r, r')\f$ from real to reciprocal space
  !! \f$f(G, G')\f$.
  !!
  !! This routine allows to have different FFT grids for the coordinates.
  !!
  !! This is done in the following steps:
  SUBROUTINE fwfft6_diff(grid_type, f, dfft, dfft_p, omega)

    USE kinds,          ONLY: dp
    USE fft_interfaces, ONLY: fwfft
    USE fft_types,      ONLY: fft_type_descriptor
    USE timing_module,  ONLY: time_fwfft6

    !> grid type used for the Fourier transform, this is passed to fwfft;
    !! check the definition of fwfft for a list of options
    CHARACTER(*), INTENT(IN)    :: grid_type

    !> *on input*  the array in real space \f$f(r, r')\f$ <br>
    !! *on output* the array in reciprocal space \f$f(G, G')\f$
    COMPLEX(dp),  INTENT(INOUT) :: f(:,:)

    !> FFT descriptor - this defines the way the Fourier transform is
    !! executed, must be consistent with grid_type (used for r/G FFT)
    !! @note the code does not check explicitly if dfft and grid_type are compatible
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft

    !> FFT descriptor - this defines the way the Fourier transform is
    !! executed, must be consistent with grid_type (used for r'/G' FFT)
    !! @note the code does not check explicitly if dfft and grid_type are compatible
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft_p

    !> volume of the unit cell
    REAL(dp),     INTENT(IN)    :: omega

    !> work array for the second index (which is not contigous in memory
    COMPLEX(dp), ALLOCATABLE    :: work(:)

    !> number of points in reciprocal space
    INTEGER num_g, num_g_p

    !> counter on reciprocal space points
    INTEGER ig

    !> number of points in real space
    INTEGER num_r, num_r_p

    !> counter on real space points
    INTEGER ir

    !> check for error in allocation
    INTEGER ierr

    CALL start_clock(time_fwfft6)

    ! initialize helper variables
    num_g   = SIZE(dfft%nl)
    num_g_p = SIZE(dfft_p%nl)
    num_r   = dfft%nnr
    num_r_p = dfft_p%nnr

    !
    ! check that f has size compatible with dfft
    IF (SIZE(f, 1) /= num_r) &
      CALL errore(__FILE__, "size of array not compatible with 1st Fourier transform", 1)
    IF (SIZE(f, 2) /= num_r_p) &
      CALL errore(__FILE__, "size of array not compatible with 2nd Fourier transform", 1)

    ! create work array
    ALLOCATE(work(num_r), STAT = ierr)
    IF (ierr /= 0) &
      CALL errore(__FILE__, "error allocating the work array for 1st FFT", ierr)

    ! loop over second index
    DO ir = 1, num_r_p
      !!
      !! 1. we store the conjugate of a row in the work array \f$w(r) = f^\ast(r, r')\f$
      !!
      work = CONJG(f(:, ir))
      !!
      !! 2. we transform the work array \f$w(r) \rightarrow w(G)\f$
      !!
      CALL fwfft(grid_type, work, dfft)
      !!
      !! 3. we extract the G vectors inside of the sphere \f$f(G, r') = w^\ast(G)\f$
      !!
      f(:num_g, ir) = CONJG(work(dfft%nl))
      !!
    END DO ! ir

    ! create work array
    DEALLOCATE(work)
    ALLOCATE(work(num_r_p), STAT = ierr)
    IF (ierr /= 0) &
      CALL errore(__FILE__, "error allocating the work array for 2nd FFT", ierr)

    ! loop over first index (now in reciprocal space)
    DO ig = 1, num_g
      !!
      !! 4. we store a column in the work array \f$w(r') = f(G, r')\f$
      !!
      work = f(ig, :)
      !!
      !! 5. we transform the work array \f$w(r') \rightarrow w(G')\f$
      !!
      CALL fwfft(grid_type, work, dfft_p)
      !!
      !! 6. extract the G vectors within the sphere \f$f(G, G') = w(G')\f$
      !!
      f(ig, :num_g_p) = work(dfft_p%nl) * omega
      !!
    END DO ! ig

    CALL stop_clock(time_fwfft6)

  END SUBROUTINE fwfft6_diff

  !> Transform an input array \f$f(G, G')\f$ from reciprocal to real space
  !! \f$f(r, r')\f$.
  !!
  !! This call internally the invfft6_diff routine.
  SUBROUTINE invfft6_same(grid_type, f, dfft, omega)

    USE kinds,          ONLY: dp
    USE fft_interfaces, ONLY: invfft
    USE fft_types,      ONLY: fft_type_descriptor

    !> grid type used for the Fourier transform, this is passed to invfft;
    !! check the definition of invfft for a list of options
    CHARACTER(*), INTENT(IN)    :: grid_type

    !> *on input*  the array in reciprocal space \f$f(G, G')\f$ <br>
    !! *on output* the array in real space \f$f(r, r')\f$
    COMPLEX(dp),  INTENT(INOUT) :: f(:,:)

    !> FFT descriptor - this defines the way the Fourier transform is
    !! executed, must be consistent with grid_type
    !! @note the code does not check explicitly if dfft and grid_type are compatible
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft

    !> volume of the unit cell
    REAL(dp),     INTENT(IN)    :: omega

    ! do the FFT via the routine which allows for different arguments
    CALL invfft6_diff(grid_type, f, dfft, dfft, omega)
    
  END SUBROUTINE invfft6_same

  !> Transform an input array \f$f(G, G')\f$ from reciprocal to real space
  !! \f$f(r, r')\f$.
  !!
  !! This routine allows to have different FFT grids for the coordinates.
  !!
  !! This is done in the following steps:
  SUBROUTINE invfft6_diff(grid_type, f, dfft, dfft_p, omega)

    USE kinds,          ONLY: dp
    USE fft_interfaces, ONLY: invfft
    USE fft_types,      ONLY: fft_type_descriptor
    USE timing_module,  ONLY: time_invfft6

    !> grid type used for the Fourier transform, this is passed to invfft;
    !! check the definition of invfft for a list of options
    CHARACTER(*), INTENT(IN)    :: grid_type

    !> *on input*  the array in reciprocal space \f$f(G, G')\f$ <br>
    !! *on output* the array in real space \f$f(r, r')\f$
    COMPLEX(dp),  INTENT(INOUT) :: f(:,:)

    !> FFT descriptor - this defines the way the Fourier transform is
    !! executed, must be consistent with grid_type (used for the r/G FFT)
    !! @note the code does not check explicitly if dfft and grid_type are compatible
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft

    !> FFT descriptor - this defines the way the Fourier transform is
    !! executed, must be consistent with grid_type (used for the r'/G' FFT)
    !! @note the code does not check explicitly if dfft and grid_type are compatible
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft_p

    !> volume of the unit cell
    REAL(dp),     INTENT(IN)    :: omega

    !> complex constant of zero
    COMPLEX(dp),  PARAMETER     :: zero = CMPLX(0.0_dp, 0.0_dp, KIND = dp)

    !> work array for the second index (which is not contigous in memory
    COMPLEX(dp), ALLOCATABLE    :: work(:)

    !> number of points in reciprocal space
    INTEGER num_g, num_g_p

    !> counter on reciprocal space points
    INTEGER ig

    !> number of points in real space
    INTEGER num_r, num_r_p

    !> counter on real space points
    INTEGER ir

    !> check for error in allocation
    INTEGER ierr

    CALL start_clock(time_invfft6)

    ! initialize helper variables
    num_g   = SIZE(dfft%nl)
    num_g_p = SIZE(dfft_p%nl)
    num_r   = dfft%nnr
    num_r_p = dfft_p%nnr

    ! check that f has size compatible with dfft
    IF (SIZE(f, 1) /= num_r) &
      CALL errore(__FILE__, "size of array not compatible with 1st Fourier transform", 1)
    IF (SIZE(f, 2) /= num_r_p) &
      CALL errore(__FILE__, "size of array not compatible with 2nd Fourier transform", 1)

    ! create work array for r'/G' FFT
    ALLOCATE(work(num_r_p), STAT = ierr)
    IF (ierr /= 0) &
      CALL errore(__FILE__, "error allocating the work array", ierr)

    ! loop over first index
    DO ig = 1, num_g
      !!
      !! 1. we store the conjugate of a column in the work array \f$w(G') = f(G, G')\f$
      !!
      ! note - we initialize work, because map only contains a subset of all entries
      work = zero
      work(dfft_p%nl) = f(ig, :num_g_p) / omega
      !!
      !! 2. we transform the work array \f$w(G') \rightarrow w(r')\f$
      !!
      CALL invfft(grid_type, work, dfft_p)
      !!
      !! 3. we extract the values from the work array \f$f(G, r') = w(r')\f$
      !!
      f(ig, :) = work
      !!
    END DO ! ig

    ! create work array for r/G FFT
    DEALLOCATE(work)
    ALLOCATE(work(num_r), STAT = ierr)
    IF (ierr /= 0) &
      CALL errore(__FILE__, "error allocating the work array", ierr)

    ! loop over second index (now in real space)
    DO ir = 1, num_r_p
      !!
      !! 4. we store a row in the work array \f$w(G) = f^\ast(G, r')\f$
      !!
      ! note - we initialize work, because map only contains a subset of all entries
      work = zero
      work(dfft%nl) = CONJG(f(:num_g, ir))
      !!
      !! 5. we transform the work array \f$w(G) \rightarrow w(r')\f$
      !!
      CALL invfft(grid_type, work, dfft)
      !!
      !! 6. extract the values from the work array \f$f(r, r') = w^\ast(G')\f$
      !!
      f(:, ir) = CONJG(work)
      !!
    END DO ! ir

    CALL stop_clock(time_invfft6)

  END SUBROUTINE invfft6_diff

  !> generate mapping from local to global G indices
  SUBROUTINE fft_map_generate(dfft, gvec, fft_map)

    USE fft_types, ONLY: fft_type_descriptor, fft_stick_index

    !> the FFT for which the mapping is generated
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft

    !> the G vectors in crystal coordinates
    INTEGER, INTENT(IN) :: gvec(:,:)

    !> the resulting FFT map
    INTEGER, ALLOCATABLE :: fft_map(:)

    !> local and global index on G vectors
    INTEGER local, globl

    ALLOCATE(fft_map(dfft%ngm))

    globl = 0
    DO local = 1, dfft%ngm
      DO
        globl = globl + 1
        IF (fft_stick_index(dfft, gvec(1,globl), gvec(2,globl)) /= 0) EXIT
      END DO
      fft_map(local) = globl
    END DO

  END SUBROUTINE fft_map_generate

END MODULE fft6_module

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
!> Provide wrapper routines around commonly used LAPACK functions.
MODULE lapack_wrapper

  IMPLICIT NONE

  !> LAPACK flag to determine work size
  INTEGER, PARAMETER :: determine = -1

CONTAINS

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

END MODULE lapack_wrapper

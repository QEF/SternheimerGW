!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
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
!> This module provides the routines to set up the linear problem from input.
MODULE linear_problem_module

  IMPLICIT NONE

CONTAINS

  !> Read the shifted linear problem written by the debug routine
  !!
  !! The linear problem is defined by
  !! \f{equation}{
  !!   (A + \sigma I) x = b
  !! \f}
  !! with the idea to solve the shifted systems in the same step as
  !! the unshifted one.
  SUBROUTINE linear_problem_read(filename, AA, sigma, bb, xx)

    USE kinds,       ONLY: dp
    USE iotk_module, ONLY: iotk_free_unit, iotk_scan_dat, &
                           iotk_open_read, iotk_close_read

    !> the filename from which the data is read
    CHARACTER(*), INTENT(IN) :: filename

    !> the unit of the file
    INTEGER iunit

    !> dimension of the linear problem
    INTEGER vec_size

    !> number of shifts in the file
    INTEGER num_shift

    !> the linear operator
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: AA(:,:)

    !> the shifts of the linear problem
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: sigma(:)

    !> the right hand side of the equation
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: bb(:)

    !> the incorrect solution of the problem we want to test
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: xx(:,:)

    !
    ! read the data from the given file
    !
    CALL iotk_free_unit(iunit)
    CALL iotk_open_read(iunit, filename, binary = .TRUE.)
    !
    ! read dimension
    CALL iotk_scan_dat(iunit, 'DIMENSION', vec_size)
    CALL iotk_scan_dat(iunit, 'NUMBER_SHIFT', num_shift)
    !
    ! create arrays to contain the data of the file
    ALLOCATE(sigma(num_shift))
    ALLOCATE(AA(vec_size, vec_size))
    ALLOCATE(bb(vec_size))
    ALLOCATE(xx(vec_size, num_shift))
    !
    ! read the linear problem (A + sigma I) x = b
    CALL iotk_scan_dat(iunit, 'LIST_SHIFT', sigma)
    CALL iotk_scan_dat(iunit, 'LINEAR_OPERATOR', AA)
    CALL iotk_scan_dat(iunit, 'RIGHT_HAND_SIDE', bb)
    CALL iotk_scan_dat(iunit, 'INCORRECT_SOLUTION', xx)
    !
    CALL iotk_close_read(iunit)

  END SUBROUTINE linear_problem_read

END MODULE linear_problem_module

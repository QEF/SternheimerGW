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
!> This module provides a routine, so that the test output looks homogenous
MODULE test_log_module

  IMPLICIT NONE

  !> unit to which the output is written
  INTEGER, PARAMETER :: test_log = 0

CONTAINS

  !> print that a test was started
  SUBROUTINE test_log_start(filename, line_number, message, rank, size)

    !> name of the file starting the test
    CHARACTER(*), INTENT(IN) :: filename

    !> line number in which the test starts
    INTEGER,      INTENT(IN) :: line_number

    !> message documenting the start of a test
    CHARACTER(*), INTENT(IN) :: message

    !> id of the current process (default to 0)
    INTEGER, INTENT(IN), OPTIONAL :: rank

    !> size of the process grid (default to 1)
    INTEGER, INTENT(IN), OPTIONAL :: size

    !> string to indicate the process grid
    CHARACTER(:), ALLOCATABLE :: proc_info

    ! only the head node prints
    IF (PRESENT(rank)) THEN
      IF (rank > 0) RETURN
    END IF

    ! determine process grid info
    IF (PRESENT(size)) THEN
      ALLOCATE(CHARACTER(17) :: proc_info)
      WRITE(proc_info, '(a6,i1,a10)') ' with ', size, ' processes'
    ELSE
      proc_info = ''
    END IF

    ! empty line to separate from previous output
    WRITE(test_log, *)

    ! write the message to document start of the routine
    WRITE(test_log, '(a)') message

    ! write the filename which starts the test
    WRITE(test_log, '(a, i0, a)') 'In file ' // filename // ', line ', line_number, proc_info

  END SUBROUTINE test_log_start

END MODULE test_log_module

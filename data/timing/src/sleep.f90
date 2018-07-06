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
!> this module provides a wrapper to the C sleep routine
MODULE sleep_module

  IMPLICIT NONE

  !> number of seconds in one minute
  INTEGER, PARAMETER :: one_min = 60

  !> number of seconds in two minutes
  INTEGER, PARAMETER :: two_min = 120

  !> number of second in five minutes
  INTEGER, PARAMETER :: five_min = 300

  !> number of seconds in 10 minutes
  INTEGER, PARAMETER :: ten_min = 600

  !> number of seconds in half an hour
  INTEGER, PARAMETER :: half_hour = 1800

  !> number of seconds in one hour
  INTEGER, PARAMETER :: one_hour = 3600

  INTERFACE

    !> sleep for specified amount of seconds
    FUNCTION sleep(second) BIND(C, NAME = "sleep") RESULT (remain)

      USE, INTRINSIC :: iso_c_binding, ONLY: C_int

      !> the number of seconds to sleep
      INTEGER(C_int), INTENT(IN), VALUE :: second

      !> the time remaining in the sleep interval (if interrupted)
      INTEGER(C_int) remain

    END FUNCTION sleep

  END INTERFACE

END MODULE sleep_module

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
!> Provides the input/output routines to read and write sigma to file.
!!
!! To access Sigma from different parts of the code consistently, we write
!! some metadata into the file that can be used to access the correct element.
!!
MODULE sigma_io_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

CONTAINS

  !> Write Sigma to disk
  SUBROUTINE sigma_io_write
  END SUBROUTINE sigma_io_write

  !> Read Sigma from disk
  SUBROUTINE sigma_io_read
  END SUBROUTINE sigma_io_read

END MODULE sigma_io_module

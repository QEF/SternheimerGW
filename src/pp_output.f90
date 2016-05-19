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
!> Provides routines that generate output for QE's plotband program
MODULE pp_output_mod

  USE kinds, ONLY : dp

  IMPLICIT NONE

  !> Print the an array for QE's plotband program
  !!
  !! The array has the shape (a1, ..., an, b), where a1 to an are
  !! multiplied to form the nbnd flag for plotband. b is the number
  !! of k-points of the system.
  !! \param filename file to which the data is printed
  !! \param kpt array containing the k-points
  !! \param data array containing the data
  INTERFACE pp_output
    MODULE PROCEDURE pp_output_2d, pp_output_3d
  END INTERFACE pp_output

  PRIVATE pp_output_2d, pp_output_3d

CONTAINS

  !> specialization of the interface for 2d data
  SUBROUTINE pp_output_2d(filename, kpt, data)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(dp),         INTENT(IN) :: kpt(:,:)
    REAL(dp),         INTENT(IN) :: data(:,:)

  END SUBROUTINE pp_output_2d

  !> specialization of the interface for 3d data
  SUBROUTINE pp_output_3d(filename, kpt, data)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(dp),         INTENT(IN) :: kpt(:,:)
    REAL(dp),         INTENT(IN) :: data(:,:,:)

  END SUBROUTINE pp_output_3d

END MODULE pp_output_mod

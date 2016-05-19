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

  USE kinds,       ONLY : dp
  USE gw_type_mod, ONLY : pp_output_type

  IMPLICIT NONE

  !> Print the an array for QE's plotband program
  !!
  !! The array has the shape (a1, ..., an), where a1 to an are
  !! multiplied to form the nbnd flag for plotband. It is assumed that
  !! this routine is called nks times, where nks is the number of k-points
  !! specified when opening the file.
  !!
  !! \param output file (as pp_output_type) to which the data is printed
  !! \param kpt vector containing the current k-point
  !! \param data array containing the data
  INTERFACE pp_output
    MODULE PROCEDURE pp_output_1d, pp_output_2d
  END INTERFACE pp_output

  PRIVATE pp_output_1d, pp_output_2d

CONTAINS

  !> specialization of the interface for 1d data
  SUBROUTINE pp_output_1d(output, kpt, data)

    TYPE(pp_output_type), INTENT(IN) :: output
    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: data(:)

    LOGICAL opnd

    !
    ! sanity test of the input
    !
    CALL errore(__FILE__, 'data array size inconsistent', output%num_band - SIZE(data))
    INQUIRE(UNIT = output%iunit, OPENED = opnd)
    IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)

    !
    ! write the data to the file
    !
    WRITE(output%iunit, '(5x,3f10.6)') kpt
    WRITE(output%iunit, '(10f10.5)') data
    ! add an empty line at the end of one data set
    WRITE(output%iunit,*)

  END SUBROUTINE pp_output_1d

  !> specialization of the interface for 2d data
  SUBROUTINE pp_output_2d(output, kpt, data)

    TYPE(pp_output_type), INTENT(IN) :: output
    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: data(:,:)

    INTEGER ii
    LOGICAL opnd

    !
    ! sanity test of the input
    !
    CALL errore(__FILE__, 'data array size inconsistent', output%num_band - SIZE(data))
    INQUIRE(UNIT = output%iunit, OPENED = opnd)
    IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)

    !
    ! write the data to the file
    !
    WRITE(output%iunit, '(5x,3f10.6)') kpt
    DO ii = 1, SIZE(data,2)
      WRITE(output%iunit, '(10f10.5)') data(:,ii)
      ! add an empty line if data would fill line completely
      IF (MOD(SIZE(data,1), 10) == 0) WRITE(output%iunit,*)
    END DO ! ii

    ! add an empty line at the end of one data set
    WRITE(output%iunit,*)

  END SUBROUTINE pp_output_2d

END MODULE pp_output_mod

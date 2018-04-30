!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
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
SUBROUTINE close_gwq(flag)
  !----------------------------------------------------------------------------
  !
  ! ... Close all files.
  ! ... Called at the end of the run with flag=.TRUE. (removes 'recover')
  ! ... or during execution with flag=.FALSE. (does not remove 'recover')
  !
  USE buffers,         ONLY : close_buffer
  USE control_flags,   ONLY : twfcollect
  USE control_gw,      ONLY : do_coulomb
  USE units_gw,        ONLY : iuwfc, iubar

  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: flag

  IF (twfcollect) THEN
    CALL close_buffer(iuwfc, 'delete')
  ELSE
    CALL close_buffer(iuwfc, 'keep')
  END IF

  IF (flag .AND. do_coulomb) THEN
    CALL close_buffer(iubar, 'keep')
  ELSE
    CALL close_buffer(iubar, 'keep')
  END IF

END SUBROUTINE close_gwq

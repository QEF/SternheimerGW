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
SUBROUTINE stop_gw(flag)
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  USE environment,   ONLY : environment_end
  USE mp_global,     ONLY : mp_global_end
  USE timing_module, ONLY : timing_print_clock
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  !
  !
  CALL timing_print_clock() 
  !
  CALL environment_end('SternheimerGW')
  !
  CALL mp_global_end()
  !
  IF ( flag ) THEN
     !
     STOP
     !
  ELSE
     !
     STOP 1
     !
  ENDIF
  !
END SUBROUTINE stop_gw

SUBROUTINE stop_smoothly_gw(flag)

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: flag

  CALL close_gwq(.FALSE.)

  CALL stop_gw(flag)

END SUBROUTINE stop_smoothly_gw

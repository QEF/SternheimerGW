!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_gw( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  USE kinds, ONLY : DP
  USE mp_global, ONLY : mp_global_end
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  !
  !
  CALL print_clock_gw()
  !
  CALL mp_global_end()
  !
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

CALL collect_grid_files()

CALL close_gwq(.FALSE.)

CALL stop_gw(flag)

END SUBROUTINE stop_smoothly_gw

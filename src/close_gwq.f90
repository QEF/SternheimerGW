!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_gwq( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files.
  ! ... Called at the end of the run with flag=.TRUE. (removes 'recover')
  ! ... or during execution with flag=.FALSE. (does not remove 'recover')
  !
  USE io_files,      ONLY : iunigk
  USE control_flags, ONLY : twfcollect
  USE paw_variables, ONLY : okpaw
  USE io_global,     ONLY : ionode, stdout
  USE buffers,       ONLY : close_buffer
  USE uspp,          ONLY : okvan
  USE units_gw,      ONLY : iuwfc, iudwf, iudwfp, iudwfm, iubar, iudrhous, iuebar, iudrho, &
                            iudvscf
  USE output,        ONLY : fildrho, fildvscf
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: exst, opnd
  !
  !
  IF ( twfcollect ) THEN
     !
     CALL close_buffer(iuwfc,'delete')
     !
  ELSE
     !
     CALL close_buffer(iuwfc,'keep')
     !
  END IF
  !
  IF (flag) THEN
     CALL close_buffer(iudwf,  'delete')
     CALL close_buffer(iudwfp, 'delete')
     CALL close_buffer(iudwfm, 'delete')
     CALL close_buffer(iubar,  'delete')
     !
     IF ( okvan ) CALL close_buffer(iudrhous,'delete')
     !
  ELSE
     CALL close_buffer(iudwf,  'keep')
     CALL close_buffer(iudwfp, 'keep')
     CALL close_buffer(iudwfm, 'keep')
     CALL close_buffer(iubar,  'keep')
     !
     IF ( okvan ) CALL close_buffer(iudrhous,'keep')
     !
  ENDIF
  !
  IF ( ionode .AND. fildrho /= ' ') THEN
     INQUIRE( UNIT=iudrho, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iudrho, STATUS = 'KEEP' )
  ENDIF
  !
  !
  IF ( fildvscf /= ' ' ) THEN
     INQUIRE( UNIT=iudvscf, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iudvscf, STATUS = 'KEEP' )
  ENDIF
  !
  !
  INQUIRE( UNIT=iunigk, OPENED=opnd ) 
  IF (opnd) CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  !
  RETURN
  !
END SUBROUTINE close_gwq

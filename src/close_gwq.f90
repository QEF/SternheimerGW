!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
SUBROUTINE close_gwq( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files.
  ! ... Called at the end of the run with flag=.TRUE. (removes 'recover')
  ! ... or during execution with flag=.FALSE. (does not remove 'recover')
  !
  USE buffers,         ONLY : close_buffer
  USE control_flags,   ONLY : twfcollect
  USE control_gw,      ONLY : do_coulomb
  USE io_global,       ONLY : ionode
  USE output_mod,      ONLY : fildrho, fildvscf
  USE units_gw,        ONLY : iuwfc, iudwf, iudwfp, iudwfm, iubar, iudrhous, iudrho, &
                              iudvscf
  USE uspp,            ONLY : okvan

  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: opnd
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
  IF (flag) THEN
    IF (do_coulomb) THEN
     CALL close_buffer(iudwf,  'keep')
     CALL close_buffer(iudwfp, 'keep')
     CALL close_buffer(iudwfm, 'keep')
     CALL close_buffer(iubar,  'keep')
    ENDIF
    CONTINUE
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
  RETURN
  !
END SUBROUTINE close_gwq

!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE clean_pw_gw(iq)
  !-----------------------------------------------------------------------
  !
  ! This routine deallocates all the variables of pwscf and of the
  ! GW code and resets the same variables as after reading input in
  ! gwq_readin, so that it is possible to start a calculation at
  ! a new q.
  !
  USE kinds,           ONLY : DP
  USE control_flags,   ONLY : twfcollect
  USE modes,           ONLY : nirr
  USE partial,         ONLY : done_irr
  USE disp,            ONLY : done_iq
  USE control_gw,      ONLY : done_bands, rec_code_read
  USE save_gw,         ONLY : restore_gw_input_variables
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !
  INTEGER :: irr
  !
  done_bands=.FALSE.
  done_iq(iq)=1
  DO irr=1,nirr
     IF (done_irr(irr)==0) done_iq(iq)=0
  ENDDO

  twfcollect=.FALSE. 

  CALL clean_pw( .FALSE. )
! CALL clean_pw( .TRUE. )

  CALL deallocate_gwq()
  rec_code_read=-1000
  !
  ! ... Close the files
  !
  CALL close_gwq( .TRUE. )
  !
  !HL not sure this is necessary
  !all stems from partial module which should probably be trashed... 
  !CALL restore_gw_input_variables()
  !
RETURN
END SUBROUTINE clean_pw_gw

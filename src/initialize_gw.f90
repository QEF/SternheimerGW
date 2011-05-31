!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE initialize_gw()
  !-----------------------------------------------------------------------
  !
  ! This is a driver to the GW initialization routine. Again ever so slightly adapted from phonon initialization
  !
  USE klist,  ONLY : nks
  USE qpoint, ONLY : nksq, ikks, ikqs
  USE control_gw, ONLY : lgamma
  USE gw_restart, ONLY : gw_writefile
  !

  IMPLICIT NONE
  INTEGER :: ik

  !
  ! ... nksq is the number of k-points, NOT including k+q points
  !

     nksq = nks / 2
     ALLOCATE(ikks(nksq), ikqs(nksq))

     DO ik=1,nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  !
  !  Save again the status of the run because now the bands have been
  !  calculated
  !
  CALL gw_writefile('init',0)
  !
  !Allocate the gw variables
  !
  CALL allocate_gwq()
  !
  !Set the main control variable of the gw code
  !
  CALL gwq_setup()
  !
  !Recover the status if available
  !
  CALL gwq_recover()
  !
  !Output summary of the main variables of the gw code
  !
  !HL want to cutdown on ouput for now  
  !CALL gwq_summary()
  
  !  Open Relevant GW files
  !  Moving this to beginning of program

  CALL openfilq()
  
  !  Initialize all quantities which do not depend on the 
  !  linear response to the perturbation
  !  All necessary quantities to describe the local and nonlocal 
  !  pseudopotential in the GW-q program.
  !
  CALL gwq_init()
  !
  CALL print_clock( 'GW' )
  !
  RETURN
  
END SUBROUTINE initialize_gw

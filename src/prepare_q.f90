!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE prepare_q(do_band, do_iq, setup_pw, iq, minq)
  !-----------------------------------------------------------------------
  !
  !  This routine prepares a few variables that are needed to control
  !  the GW run after the q point has been decided, but before
  !  doing the band calculation. 
  !  In particular if ldisp=true it sets:
  !  xq : the q point for the q calculation
  !  current_iq : the current q point
  !  do_iq : if .true. q point has to be calculated
  !  setup_pw : if .true. the pw_setup has to be run
  !  do_band : if .true. the bands need to be calculated before phonon
  
  USE control_flags,   ONLY : modenum
  USE io_global,       ONLY : stdout, ionode
  USE klist,           ONLY : lgauss
  USE qpoint,          ONLY : xq
  USE disp,            ONLY : x_q, comp_iq, xk_kpoints
  USE control_gw,      ONLY : ldisp, lgamma, current_iq, done_bands, do_epsil
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  LOGICAL, INTENT(INOUT) :: minq
  LOGICAL, INTENT(OUT) :: do_band, do_iq, setup_pw
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER :: irr
  !
  !
  current_iq = iq
  !
  IF ( ldisp ) THEN
     ! ... set the name for the output file
     ! ... set the q point
     !For convolution we calculate  W_{-\q}
     !  xq(1:3)  = -x_q(1:3,iq)
        xq(1:3)  = x_q(1:3,iq)
        if(.not.lgauss) then
          !if ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 ) xq(1) = -0.001d0
          if ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 ) xq(1) = 0.001d0
        endif
     !In case we want to calulate eps(q) where q is given in the input file:
        if (do_epsil) xq(:) = xk_kpoints(:, iq)
        lgamma = (xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0)
  ENDIF
  !
  !
  ! ... In the case of q != 0, we make first a non selfconsistent run
  !
  RETURN
END SUBROUTINE prepare_q

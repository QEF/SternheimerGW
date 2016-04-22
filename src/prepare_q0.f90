!------------------------------------------------------------------------------ 
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
SUBROUTINE prepare_q0(do_band, do_iq, setup_pw, iq, minq)
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
        xq(1:3)  = x_q(1:3,iq)
        if ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 ) xq(1) = 0.001d0
     !In case we want to calulate eps(q) where q is given in the input file:
        if (do_epsil) xq(:) = xk_kpoints(:, iq)
        lgamma = (xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0)
  ENDIF
  do_band=.true.
  setup_pw=.true.
  do_iq=.true.
  !
  !
  ! ... In the case of q != 0, we make first a non selfconsistent run
  !
  RETURN
END SUBROUTINE prepare_q0

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
SUBROUTINE initialize_gw(coulomb)
  !-----------------------------------------------------------------------
  !
  ! This is a driver to the GW initialization routine. Again ever so slightly adapted from phonon initialization
  !
  USE klist,  ONLY : nks
  USE qpoint, ONLY : nksq, ikks, ikqs
  USE control_lr, ONLY : lgamma

  IMPLICIT NONE

  !> are we in the Coulomb or the self-energy step?
  LOGICAL, INTENT(IN) :: coulomb

  INTEGER :: ik

  !
  ! ... nksq is the number of k-points, NOT including k+q points
  !
  IF (coulomb) THEN
    !
    ! for the screened Coulomb interaction, k and k + q alternate for
    ! all points except Gamma, where k and k + q are identical
    !
    IF ( lgamma ) THEN
      !
      nksq = nks
      ALLOCATE(ikks(nksq), ikqs(nksq))
      DO ik=1,nksq
        ikks(ik) = ik
        ikqs(ik) = ik
      ENDDO
     !
    ELSE
     !
      nksq = nks / 2
      ALLOCATE(ikks(nksq), ikqs(nksq))
      DO ik=1,nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
      ENDDO
      !
    END IF
    !
  ELSE ! not coulomb
    !
    ! nksq, ikks, and ikqs is already set in setup_nscf_green
    !
  END IF
  !
  !Allocate the gw variables
  !
  CALL allocate_gwq()
  !
  !Set the main control variable of the gw code
  !
  CALL gwq_setup()
  !
  !Open Relevant GW files
  CALL openfilq()
  !Initialize all quantities which do not depend on the 
  !linear response to the perturbation
  !All necessary quantities to describe the local and nonlocal 
  !pseudopotential in the GW program.
  !
  CALL gwq_init(coulomb)
  !
  CALL print_clock( 'GW' )
  !
  RETURN
END SUBROUTINE initialize_gw

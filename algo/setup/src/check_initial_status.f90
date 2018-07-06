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
SUBROUTINE check_initial_status()
  !-----------------------------------------------------------------------
  ! This routine checks the initial status of the phonon run and sets
  ! the variables that control the run, dealing with the image
  ! and GRID parallelization features of the phonon code.
  ! The size of the grid is determined by the following variables:
  ! nqs : the number of q points
  ! x_q : the coordinates of the q points
  ! nfs : the number of imaginary frequencies
  ! fiu : which frequencies 
  ! The flags that control which tensors to calculate
  !
  USE control_gw,      ONLY : tmp_dir_gw, do_epsil
  USE io_files,        ONLY : tmp_dir
  USE io_rho_xml,      ONLY : write_scf
  USE lsda_mod,        ONLY : nspin
  USE mp,              ONLY : mp_bcast
  USE mp_global,       ONLY : mp_global_end
  USE scf,             ONLY : rho
  !
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  !
  tmp_dir=tmp_dir_gw
  !
  ! If this not a recover run, we generate the q mesh. Otherwise at this
  ! point the code has read the q mesh from the files contained in 
  ! prefix.phsave
  IF (.NOT.do_epsil) CALL q_points()
  !
  ! this is the standard treatment
  CALL write_scf( rho, nspin )
  !
  RETURN
  END SUBROUTINE check_initial_status

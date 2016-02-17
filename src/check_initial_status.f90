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
SUBROUTINE check_initial_status(auxdyn)
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
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : modenum
  USE ions_base,       ONLY : nat
  USE io_files,        ONLY : tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE scf,             ONLY : rho
  USE disp,            ONLY : nqs, x_q, nq1, nq2, nq3
  USE qpoint,          ONLY : xq
  USE output,          ONLY : fildyn
  USE control_gw,      ONLY : ldisp, recover, done_bands,  &
                              start_q, last_q, current_iq, tmp_dir_gw, lgamma, &
                              ext_recover, ext_restart
  USE save_gw,         ONLY : tmp_dir_save
  USE units_gw,        ONLY : iudyn
  USE io_rho_xml,      ONLY : write_rho
  USE mp_images,       ONLY : nimage, intra_image_comm
  USE io_global,       ONLY : ionode, ionode_id
  USE io_files,        ONLY : prefix
  USE mp,              ONLY : mp_bcast
  USE xml_io_base,     ONLY : create_directory
  USE mp_global,       ONLY : mp_global_end
  !
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: auxdyn, filename
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER :: iq, iq_start, ierr
  !
  tmp_dir=tmp_dir_gw
  !
  ! If this not a recover run, we generate the q mesh. Otherwise at this
  ! point the code has read the q mesh from the files contained in 
  ! prefix.phsave
  call q_points()
  !
  ! this is the standard treatment
  CALL write_rho( rho, nspin )
  !
  RETURN
  END SUBROUTINE check_initial_status

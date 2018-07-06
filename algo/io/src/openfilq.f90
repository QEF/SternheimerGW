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
SUBROUTINE openfilq()
!----------------------------------------------------------------------------
! ... This subroutine opens all the files necessary for the GW q
! ... calculation.
 
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level
  USE control_gw,       ONLY : tmp_dir_gw, do_coulomb
  USE control_lr,       ONLY : lgamma
  USE io_files,         ONLY : tmp_dir, diropn, seqopn, prefix
  USE noncollin_module, ONLY : npol
  USE save_gw,          ONLY : tmp_dir_save
  USE units_gw,         ONLY : iuwfc, lrwfc, iubar, lrbar
  USE wvfct,            ONLY : nbnd, npwx
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst, exst_mem
  ! logical variable to check file existe
  !
  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfilq', 'wrong prefix', 1)
  !     There are six direct access files to be opened in the tmp area
  !     The file with the wavefunctions. In the lgamma case reads those
  !     written by pw.x. In the other cases those calculated by gw.x
  tmp_dir=tmp_dir_gw
  IF (lgamma) tmp_dir=tmp_dir_save

  iuwfc = 20
  lrwfc = nbnd * npwx * npol
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir_gw)
  IF (.NOT.exst.AND..NOT.exst_mem) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF

  !
  ! From now on all files are written with the _gw prefix
  !
  tmp_dir = tmp_dir_gw
  !
  !    The file with deltaV_{bare} * psi
  !
  IF (do_coulomb) THEN
    iubar = 21
    lrbar = nbnd * npwx * npol
    CALL open_buffer(iubar, 'bar', lrbar, io_level, exst_mem, exst, tmp_dir)
  END IF

END SUBROUTINE openfilq

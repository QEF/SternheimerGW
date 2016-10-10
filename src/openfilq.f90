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
SUBROUTINE openfilq()
!----------------------------------------------------------------------------
! ... This subroutine opens all the files necessary for the GW q
! ... calculation.
 
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level
  USE control_gw,       ONLY : ext_recover, tmp_dir_gw, lgamma, do_coulomb
  USE fft_base,         ONLY : dfftp
  USE io_files,         ONLY : tmp_dir, diropn, seqopn, prefix
  USE noncollin_module, ONLY : npol, nspin_mag
  USE save_gw,          ONLY : tmp_dir_save
  USE units_gw,         ONLY : iuwfc, iudwf, iubar, iudrhous, iudrho, lrwfc, &
                               lrdwf, lrbar, lrdrhous, lrdrho, iudwfm, iudwfp
  USE uspp,             ONLY : okvan
  USE wvfct,            ONLY : nbnd, npwx
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst, exst_mem
  ! logical variable to check file existe
  !
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
  tmp_dir=tmp_dir_gw
  !
  !    The file with deltaV_{bare} * psi
  !
  if(do_coulomb) then
    iubar = 21
    lrbar = nbnd * npwx * npol
    CALL open_buffer (iubar, 'bar', lrbar, io_level, exst_mem, exst, tmp_dir)
    IF (ext_recover.AND..NOT.exst) &
       CALL errore ('openfilq','file '//trim(prefix)//'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !
    iudwf = 22
    lrdwf =  nbnd * npwx * npol
    CALL open_buffer (iudwf, 'dwf', lrdwf, io_level, exst_mem, exst, tmp_dir)
    IF (ext_recover.AND..NOT.exst) &
    CALL errore ('openfilq','file '//trim(prefix)//'.dwf not found', 1)

    iudwfm = 29 
    iudwfp = 30
    CALL open_buffer (iudwfp, 'dwfp', lrwfc, io_level, exst_mem, exst, tmp_dir)
    CALL open_buffer (iudwfm, 'dwfm', lrwfc, io_level, exst_mem, exst, tmp_dir)
  endif
  !
  !   open a file with the static change of the charge
  !
    IF (okvan) THEN
     iudrhous = 25
     lrdrhous =  dfftp%nnr * nspin_mag
     CALL open_buffer (iudrhous, 'prd', lrdrhous, io_level, exst_mem, exst,tmp_dir)
     IF (ext_recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.prd not found', 1)
  ENDIF
  !
  !  Optional file(s) containing Delta\rho (opened and written in solve_e
  !  and solve_linter). Used for third-order calculations.
  !
  iudrho = 23
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
!HL write files for \Delta\psi^{\pm}
  RETURN
END SUBROUTINE openfilq

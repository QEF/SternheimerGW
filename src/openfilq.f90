!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.
!
!----------------------------------------------------------------------------
SUBROUTINE openfilq()
!----------------------------------------------------------------------------
! ... This subroutine opens all the files necessary for the GW q
! ... calculation.
 
  USE kinds,          ONLY : DP
  USE control_flags, ONLY : modenum
  USE units_gw,       ONLY : iuwfc, iudwf, iubar, iucom, iudvkb3, &
                             iudrhous, iuebar, iudrho, iudyn, iudvscf, &
                             lrwfc, lrdwf, lrbar, lrcom, lrdvkb3, lrcoul, &
                             lrdrhous, lrebar, lrdrho, iudwfm, iudwfp, iuncoul, lrgrn, iungreen
  USE io_files,       ONLY : tmp_dir
  USE freq_gw,        ONLY : nfs, nwgreen
  USE gsmooth,        ONLY : ngms

!HL  USE control_ph,     ONLY : epsil, zue, ext_recover, trans, elph, lgamma, &
!                               tmp_dir_ph, start_irr, last_irr

  USE control_gw,     ONLY : ext_recover, trans, tmp_dir_gw
  USE save_gw,        ONLY : tmp_dir_save
  USE qpoint,         ONLY : nksq
  USE output,         ONLY : fildyn, fildvscf
  USE wvfct,          ONLY : nbnd, npwx
  USE gvect,          ONLY : nrx1, nrx2, nrx3, nrxx
  USE lsda_mod,       ONLY : nspin
  USE uspp,           ONLY : nkb, okvan
  USE io_files,       ONLY : prefix, iunigk
  USE noncollin_module, ONLY : npol, nspin_mag
  USE control_flags,  ONLY : twfcollect
  USE mp_global,      ONLY : me_pool
  USE io_global,      ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint
  ! the name of the file
  LOGICAL :: exst
  ! logical variable to check file existe
  !
  REAL(DP) :: edum(1,1), wdum(1,1)
  INTEGER :: ndr, ierr
  !
  !
  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfilq', 'wrong prefix', 1)
  !     There are six direct access files to be opened in the tmp area
  !     The file with the wavefunctions. In the lgamma case reads those
  !     written by pw.x. In the other cases those calculated by gw.x
 
  tmp_dir=tmp_dir_gw

  !IF (lgamma.AND.modenum==0) tmp_dir=tmp_dir_save

  iuwfc = 20
  lrwfc = 2 * nbnd * npwx * npol

  CALL diropn (iuwfc, 'wfc', lrwfc, exst)

  IF (.NOT.exst) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  !
  ! From now on all files are written with the _gw prefix
  !
  tmp_dir=tmp_dir_gw
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx * npol
  CALL diropn (iubar, 'bar', lrbar, exst)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx * npol
  CALL diropn (iudwf, 'dwf', lrdwf, exst)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
     iudrhous = 25
     lrdrhous = 2 * nrxx * nspin_mag
     CALL diropn (iudrhous, 'prd', lrdrhous, exst)
     IF (ext_recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.prd not found', 1)
  ENDIF
  !
  !  Optional file(s) containing Delta\rho (opened and written in solve_e
  !  and solve_linter). Used for third-order calculations.
  !
  iudrho = 23
  lrdrho = 2 * nrx1 * nrx2 * nrx3 * nspin_mag
  !
  !
  !   Here the sequential files
  !
  !   The igk at a given k (and k+q if q!=0)
  !
  iunigk = 24
  IF (nksq > 1) CALL seqopn (iunigk, 'igk', 'unformatted', exst)

!HL write files for \Delta\psi^{\pm}
  iudwfm = 29 
  iudwfp = 30

  CALL diropn (iudwfp, 'dwfp', lrdwf, exst)
  CALL diropn (iudwfm, 'dwfm', lrdwf, exst)

400 IF (fildvscf.NE.' ') THEN
     iudvscf = 27
     IF ( me_pool == 0 ) THEN
        CALL diropn (iudvscf, fildvscf, lrdrho, exst)
     END IF
  END IF
  !
  RETURN
  !
END SUBROUTINE openfilq

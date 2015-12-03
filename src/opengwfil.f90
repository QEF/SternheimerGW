  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE opengwfil()
  USE units_gw,        ONLY : iuncoul, iungreen, lrgrn, lrcoul, iunsigma, &
                              lrsigma, lrsex, iunsex, iunresid, lrresid, & 
                              iunalphabeta, lralphabeta, iunsigext, lrsigext
  USE io_files,        ONLY : tmp_dir, diropn, seqopn
  USE freq_gw,         ONLY : nfs, nwsigma
  USE gwsigma,         ONLY : sigma_x_st, sigma_c_st, gcutcorr
  USE control_gw,      ONLY : multishift, do_serial, do_sigma_exxG
  USE wvfct,           ONLY : nbnd,npwx
  USE io_global,       ONLY : meta_ionode

IMPLICIT  NONE
  LOGICAL :: exst

    iuncoul = 28
   !lrcoul = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt * nfs
    lrcoul = 2 * gcutcorr * gcutcorr * nfs
    if(meta_ionode) CALL diropn (iuncoul, 'coul', lrcoul, exst)
    iungreen = 31
    !lrgrn  = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt
    lrgrn  = 2 * gcutcorr * gcutcorr
!    CALL diropn (iungreen, 'green', lrgrn, exst)
!\Sigma^{c}(\G,\G';\omega)
    iunsigma = 32
    if(do_serial) then
       !lrsigma = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt
       lrsigma = 2 * gcutcorr * gcutcorr
    else  
       !lrsigma = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt * nwsigma
       lrsigma = 2 * gcutcorr * gcutcorr*nwsigma
    endif
    if(meta_ionode) CALL diropn(iunsigma, 'sigma', lrsigma, exst)
!\Sigma^{x}(\G,\G';\omega)
    if(.not. do_sigma_exxG) then
       iunsex = 33
       lrsex = 2 * sigma_x_st%ngmt * sigma_x_st%ngmt
       if(meta_ionode) CALL diropn(iunsex, 'sigma_ex', lrsex, exst)
    endif
    if(multishift) then
       iunresid = 34
!!!!!This only needs to be the same as ngmsco!!!
!      lrresid  = 2*npwx
! Factor of 2 for complex and second factor so we don't
! throw away g sphere
       !lrresid  = 2*2*sigma_c_st%ngmt
       lrresid  = 2*2*gcutcorr
       iunalphabeta = 35
       lralphabeta  = 4
    endif
END SUBROUTINE opengwfil

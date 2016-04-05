!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
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
  USE output_mod,      ONLY : filcoul, filsigx, filsigc

IMPLICIT  NONE
  LOGICAL :: exst

    iuncoul = 28
   !lrcoul = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt * nfs
    lrcoul = 2 * gcutcorr * gcutcorr * nfs
    if(meta_ionode) CALL diropn (iuncoul, filcoul, lrcoul, exst)
    iungreen = 31
    !lrgrn  = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt
    lrgrn  = 2 * gcutcorr * gcutcorr
!    CALL diropn (iungreen, 'green', lrgrn, exst)
!\Sigma^{c}(\G,\G';\omega)
    iunsigma = 32
    !lrsigma = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt * nwsigma
    lrsigma = 2 * gcutcorr * gcutcorr*nwsigma
    if(meta_ionode) CALL diropn(iunsigma, filsigc, lrsigma, exst)
!\Sigma^{x}(\G,\G';\omega)
    if(.not. do_sigma_exxG) then
       iunsex = 33
       lrsex = 2 * sigma_x_st%ngmt * sigma_x_st%ngmt
       if(meta_ionode) CALL diropn(iunsex, filsigx, lrsex, exst)
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

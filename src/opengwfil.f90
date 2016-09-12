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

  USE control_gw,      ONLY : multishift, output
  USE disp,            ONLY : xk_kpoints, num_k_pts
  USE freq_gw,         ONLY : nfs, nwsigma
  USE gwsigma,         ONLY : sigma_x_st, gcutcorr
  USE io_files,        ONLY : diropn, seqopn
  USE io_global,       ONLY : meta_ionode
  USE output_mod,      ONLY : filcoul, filsigx, filsigc
  USE sigma_io_module, ONLY : sigma_io_open_write
  USE units_gw,        ONLY : iuncoul, lrcoul, iunsigma, &
                              lrsigma, lrsex, iunsex, iunresid, lrresid, & 
                              iunalphabeta, lralphabeta

IMPLICIT  NONE
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit

  ! open file for coulomb 
  iuncoul = 28
  lrcoul = 2 * gcutcorr * gcutcorr * nfs
  IF (meta_ionode) THEN 
    CALL diropn (iuncoul, filcoul, lrcoul, exst)
  END IF

  ! file for \Sigma^{c}(\G,\G';\omega)
  lrsigma = 2 * gcutcorr * gcutcorr*nwsigma
  IF (meta_ionode) THEN
    iunsigma = find_free_unit()
    CALL diropn(iunsigma, filsigc, lrsigma, exst)
  END IF

  ! file for \Sigma^{x}(\G,\G';\omega)
  lrsex = 2 * sigma_x_st%ngmt * sigma_x_st%ngmt
  IF (meta_ionode) THEN
    iunsex = find_free_unit()
    CALL diropn(iunsex, filsigx, lrsex, exst)
  END IF

  ! file for output of sigma
  IF (meta_ionode) THEN
    CALL sigma_io_open_write(output%file_sigma, xk_kpoints(:,:num_k_pts), &
                     sigma_x_st%ngmt, gcutcorr, nwsigma, output%unit_sigma)
  END IF

  IF (multishift) THEN
    ! This only needs to be the same as ngmsco
    ! Factor of 2 for complex and second factor so we don't
    ! throw away g sphere
    iunresid = 34
    lrresid  = 2*2*gcutcorr
    iunalphabeta = 35
    lralphabeta  = 4
  END IF

END SUBROUTINE opengwfil

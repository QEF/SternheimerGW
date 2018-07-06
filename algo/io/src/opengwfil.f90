!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
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
SUBROUTINE opengwfil(grid)

  USE control_gw,        ONLY : output
  USE disp,              ONLY : xk_kpoints, num_k_pts
  USE freq_gw,           ONLY : nfs, nwsigma
  USE io_files,          ONLY : diropn, seqopn
  USE io_global,         ONLY : meta_ionode
  USE output_mod,        ONLY : filcoul, filsigx, filsigc
  USE sigma_grid_module, ONLY : sigma_grid_type
  USE sigma_io_module,   ONLY : sigma_io_open_write
  USE units_gw,          ONLY : iuncoul, lrcoul, iunsigma, lrsigma, lrsex, iunsex

IMPLICIT  NONE

  !> the grid used for the FFT
  TYPE(sigma_grid_type), INTENT(IN) :: grid

  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit

  !> number of G vectors in exchange grid
  INTEGER num_g_exch

  !> number of G vectors in correlation grid
  INTEGER num_g_corr

  ! initialize helper variables
  num_g_exch = grid%exch_fft%ngm
  num_g_corr = grid%corr_fft%ngm

  ! open file for coulomb 
  iuncoul = 28
  lrcoul = 2 * num_g_corr * num_g_corr * nfs
  IF (meta_ionode) THEN 
    CALL diropn (iuncoul, filcoul, lrcoul, exst)
  END IF

  ! file for \Sigma^{c}(\G,\G';\omega)
  lrsigma = 2 * num_g_corr * num_g_corr
  IF (meta_ionode) THEN
    iunsigma = find_free_unit()
    CALL diropn(iunsigma, filsigc, lrsigma, exst)
  END IF

  ! file for \Sigma^{x}(\G,\G';\omega)
  lrsex = 2 * num_g_exch * num_g_exch
  IF (meta_ionode) THEN
    iunsex = find_free_unit()
    CALL diropn(iunsex, filsigx, lrsex, exst)
  END IF

  ! file for output of sigma
  IF (meta_ionode) THEN
    CALL sigma_io_open_write(output%file_sigma, xk_kpoints(:,:num_k_pts), &
                     num_g_exch, num_g_corr, nwsigma, output%unit_sigma)
  END IF

END SUBROUTINE opengwfil

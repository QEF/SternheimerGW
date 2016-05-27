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
!> Provides the input/output routines to read and write sigma to file.
!!
!! To access Sigma from different parts of the code consistently, we write
!! some metadata into the file that can be used to access the correct element.
!!
MODULE sigma_io_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  !> xml tag used as root tag
  CHARACTER(*), PARAMETER :: tag_root = "SELF_ENERGY"

  !> xml tag for the number of k-points
  CHARACTER(*), PARAMETER :: tag_num_kpoint = "NUM_KPOINT"

  !> xml tag for the number of G-vectors for exchange part of Sigma
  CHARACTER(*), PARAMETER :: tag_num_exchange = "NUM_EXCHANGE"

  !> xml tag for the number of G-vectors for correlation part of Sigma
  CHARACTER(*), PARAMETER :: tag_num_correlation = "NUM_CORRELATION"

  !> xml tag for the k-point
  CHARACTER(*), PARAMETER :: tag_kpoint = "KPOINT"

  !> xml tag for the sigma matrix
  CHARACTER(*), PARAMETER :: tag_sigma = "SIGMA"

  !> xml tag for the exchange part
  CHARACTER(*), PARAMETER :: tag_exchange = "EXCHANGE"

  !> xml tag for the correlation part
  CHARACTER(*), PARAMETER :: tag_correlation = "CORRELATION"

CONTAINS

  !> Open Sigma file to write.
  !!
  !! Open Sigma and write metadata about it to the header.
  !!
  !! \param[in] filename Name of the file in which the data is written.
  !! \param[in] kpt List of k-points for which Sigma is generated
  !! \param[in] ngm_x Number of G vectors for exchange.
  !! \param[in] ngm_c Number of G vectors for correlation.
  !! \param[out] iunit Unit to access the file.
  SUBROUTINE sigma_io_open(filename, kpt, ngm_x, ngm_c, iunit)

    USE iotk_module, ONLY: iotk_free_unit, iotk_open_write, &
                           iotk_write_dat

    CHARACTER(*), INTENT(IN)  :: filename
    REAL(dp),     INTENT(IN)  :: kpt(:,:)
    INTEGER,      INTENT(IN)  :: ngm_x, ngm_c
    INTEGER,      INTENT(OUT) :: iunit

    ! find a free unit
    CALL iotk_free_unit(iunit)

    ! open the file
    CALL iotk_open_write(iunit, filename, binary=.TRUE., root=tag_root)

    !
    ! write the meta data
    !
    ! number of G-vectors in Sigma_x
    CALL iotk_write_dat(iunit, tag_num_exchange, ngm_x)

    ! number of G-vectors in Sigma_c
    CALL iotk_write_dat(iunit, tag_num_correlation, ngm_c)

    ! number of k-points
    CALL iotk_write_dat(iunit, tag_num_kpoint, SIZE(kpt, 2))

    ! write k-points
    CALL iotk_write_dat(iunit, tag_kpoint, kpt)

  END SUBROUTINE sigma_io_open

  !> Write Sigma to disk.
  !!
  !! Write the self-energy for exchange and correlation separately.
  !!
  !! \param[in] iunit Unit to which the Sigma is written, must be opened with iotk_module.
  !! \param[in] ikpt Index of the k-point.
  !! \param[in] sigma_x Exchange part of Sigma.
  !! \param[in] sigma_c Correlation part of Sigma (frequency dependent).
  !!
  SUBROUTINE sigma_io_write(iunit, ikpt, sigma_x, sigma_c)

    USE iotk_module, ONLY: iotk_index, iotk_namlenx, &
                           iotk_write_begin, iotk_write_end, iotk_write_dat

    INTEGER,     INTENT(IN) :: iunit
    INTEGER,     INTENT(IN) :: ikpt
    COMPLEX(dp), INTENT(IN) :: sigma_x(:,:)
    COMPLEX(dp), INTENT(IN) :: sigma_c(:,:,:)

    CHARACTER(LEN=iotk_namlenx) tag_sigma_loc

    tag_sigma_loc = tag_sigma // iotk_index(iunit)

    CALL iotk_write_begin(iunit, tag_sigma_loc)

    ! write the exchange part
    CALL iotk_write_dat(iunit, tag_exchange, sigma_x)

    ! write the correlation part
    CALL iotk_write_dat(iunit, tag_correlation, sigma_c)

    CALL iotk_write_end(iunit, tag_sigma_loc)

  END SUBROUTINE sigma_io_write

  !> Read Sigma from disk
  SUBROUTINE sigma_io_read
  END SUBROUTINE sigma_io_read

END MODULE sigma_io_module

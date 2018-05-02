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

  !> xml tag for the number of G-vectors for exchange part of Sigma
  CHARACTER(*), PARAMETER :: tag_num_exchange = "NUM_EXCHANGE"

  !> xml tag for the number of G-vectors for correlation part of Sigma
  CHARACTER(*), PARAMETER :: tag_num_correlation = "NUM_CORRELATION"

  !> xml tag for the number of frequencies
  CHARACTER(*), PARAMETER :: tag_frequency = "NUM_FREQUENCY"

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
  !! \param[in] kpt List of k-points for which Sigma is generated.
  !! \param[in] ngm_x Number of G vectors for exchange.
  !! \param[in] ngm_c Number of G vectors for correlation.
  !! \param[in] num_freq Number of frequency points.
  !! \param[out] iunit Unit to access the file.
  SUBROUTINE sigma_io_open_write(filename, kpt, ngm_x, ngm_c, num_freq, iunit)

    USE iotk_module, ONLY: iotk_free_unit, iotk_open_write, &
                           iotk_write_dat

    CHARACTER(*), INTENT(IN)  :: filename
    REAL(dp),     INTENT(IN)  :: kpt(:,:)
    INTEGER,      INTENT(IN)  :: ngm_x, ngm_c
    INTEGER,      INTENT(IN)  :: num_freq
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

    ! number of frequency points
    CALL iotk_write_dat(iunit, tag_frequency, num_freq)

    ! write k-points
    CALL iotk_write_dat(iunit, tag_kpoint, kpt)

  END SUBROUTINE sigma_io_open_write

  !> Open Sigma file to read.
  !!
  !! Open Sigma and read metadata about it from the header.
  !!
  !! \param[in] filename Name of the file from which the data is read.
  !! \param[out] kpt List of k-points for which Sigma is generated.
  !! \param[out] ngm_x Number of G vectors for exchange.
  !! \param[out] ngm_c Number of G vectors for correlation.
  !! \param[out] num_freq Number of frequency points.
  !! \param[out] iunit Unit to access the file.
  SUBROUTINE sigma_io_open_read(filename, kpt, ngm_x, ngm_c, num_freq, iunit)

    USE iotk_module, ONLY: iotk_free_unit, iotk_open_read, &
                           iotk_scan_dat

    CHARACTER(*), INTENT(IN)  :: filename
    REAL(dp),     INTENT(OUT) :: kpt(:,:)
    INTEGER,      INTENT(OUT) :: ngm_x, ngm_c
    INTEGER,      INTENT(OUT) :: num_freq
    INTEGER,      INTENT(OUT) :: iunit

    ! find a free unit
    CALL iotk_free_unit(iunit)

    ! open the file
    CALL iotk_open_read(iunit, filename, binary=.TRUE.)

    !
    ! read the meta data
    !
    ! number of G-vectors in Sigma_x
    CALL iotk_scan_dat(iunit, tag_num_exchange, ngm_x)

    ! number of G-vectors in Sigma_c
    CALL iotk_scan_dat(iunit, tag_num_correlation, ngm_c)

    ! number of frequency points
    CALL iotk_scan_dat(iunit, tag_frequency, num_freq)

    ! write k-points
    CALL iotk_scan_dat(iunit, tag_kpoint, kpt)

  END SUBROUTINE sigma_io_open_read

  !> Write correlation part of Sigma to disk.
  !!
  !! Write the correlation part of the self-energy to disk. Use the separate routine
  !! sigma_io_write_x to write the exchange part.
  !!
  !! \param[in] iunit Unit to which the Sigma is written, must be opened with iotk_module.
  !! \param[in] ikpt Index of the k-point.
  !! \param[in] sigma_c Correlation part of Sigma (frequency dependent).
  !!
  !! \note The format of the output file expects that sigma_c is written first.
  !!
  SUBROUTINE sigma_io_write_c(iunit, ikpt, sigma_c)

    USE iotk_module, ONLY: iotk_index, iotk_namlenx, &
                           iotk_write_begin, iotk_write_end, iotk_write_dat

    INTEGER,     INTENT(IN) :: iunit
    INTEGER,     INTENT(IN) :: ikpt
    COMPLEX(dp), INTENT(IN) :: sigma_c(:,:,:)

    CHARACTER(LEN=iotk_namlenx) tag_sigma_loc

    tag_sigma_loc = tag_sigma // iotk_index(ikpt)

    CALL iotk_write_begin(iunit, tag_sigma_loc)

    ! write the correlation part
    CALL iotk_write_dat(iunit, tag_correlation, sigma_c)

  END SUBROUTINE sigma_io_write_c

  !> Write exchange part of Sigma to disk.
  !!
  !! Write the exchange part of the self-energy to disk. Use the separate routine
  !! sigma_io_write_c to write the correlation part.
  !!
  !! \param[in] iunit Unit to which the Sigma is written, must be opened with iotk_module.
  !! \param[in] ikpt Index of the k-point.
  !! \param[in] sigma_x Exchange part of Sigma.
  !!
  !! \note The format of the output file expects that sigma_c is written first.
  !!
  SUBROUTINE sigma_io_write_x(iunit, ikpt, sigma_x)

    USE iotk_module, ONLY: iotk_index, iotk_namlenx, &
                           iotk_write_begin, iotk_write_end, iotk_write_dat

    INTEGER,     INTENT(IN) :: iunit
    INTEGER,     INTENT(IN) :: ikpt
    COMPLEX(dp), INTENT(IN) :: sigma_x(:,:)

    CHARACTER(LEN=iotk_namlenx) tag_sigma_loc

    tag_sigma_loc = tag_sigma // iotk_index(ikpt)

    ! write the exchange part
    CALL iotk_write_dat(iunit, tag_exchange, sigma_x)

    CALL iotk_write_end(iunit, tag_sigma_loc)

  END SUBROUTINE sigma_io_write_x

  !> Read Sigma from disk
  !!
  !! Read the self-energy for exchange and correlation separately.
  !!
  !! \param[in] iunit Unit to which the Sigma is written, must be opened with iotk_module.
  !! \param[in] ikpt Index of the k-point.
  !! \param[out] sigma_x Exchange part of Sigma.
  !! \param[out] sigma_c Correlation part of Sigma (frequency dependent).
  !!
  SUBROUTINE sigma_io_read(iunit, ikpt, sigma_x, sigma_c)

    USE iotk_module, ONLY: iotk_index, iotk_namlenx, &
                           iotk_scan_begin, iotk_scan_end, iotk_scan_dat

    INTEGER,     INTENT(IN)  :: iunit
    INTEGER,     INTENT(IN)  :: ikpt
    COMPLEX(dp), INTENT(OUT) :: sigma_x(:,:)
    COMPLEX(dp), INTENT(OUT) :: sigma_c(:,:,:)

    CHARACTER(LEN=iotk_namlenx) tag_sigma_loc

    tag_sigma_loc = tag_sigma // iotk_index(ikpt)

    CALL iotk_scan_begin(iunit, tag_sigma_loc)

    ! write the correlation part
    CALL iotk_scan_dat(iunit, tag_correlation, sigma_c)

    ! write the exchange part
    CALL iotk_scan_dat(iunit, tag_exchange, sigma_x)

    CALL iotk_scan_end(iunit, tag_sigma_loc)

  END SUBROUTINE sigma_io_read

  !> Close a file opened for writing.
  !! \param[in] iunit Unit of the file to close.
  SUBROUTINE sigma_io_close_write(iunit)

    USE iotk_module, ONLY: iotk_close_write

    INTEGER, INTENT(IN) :: iunit

    CALL iotk_close_write(iunit)

  END SUBROUTINE sigma_io_close_write

  !> Close a file opened for reading.
  !! \param[in] iunit Unit of the file to close.
  SUBROUTINE sigma_io_close_read(iunit)

    USE iotk_module, ONLY: iotk_close_read

    INTEGER, INTENT(IN) :: iunit

    CALL iotk_close_read(iunit)

  END SUBROUTINE sigma_io_close_read

END MODULE sigma_io_module

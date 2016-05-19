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
!> Provides routines that generate output for QE's plotband program
MODULE pp_output_mod

  USE kinds,       ONLY : dp
  USE gw_type_mod, ONLY : pp_output_type, output_type

  IMPLICIT NONE

  !> Print the an array for QE's plotband program
  !!
  !! The array has the shape (a1, ..., an), where a1 to an are
  !! multiplied to form the nbnd flag for plotband. It is assumed that
  !! this routine is called nks times, where nks is the number of k-points
  !! specified when opening the file.
  !!
  !! \param output file (as pp_output_type) to which the data is printed
  !! \param kpt vector containing the current k-point
  !! \param data array containing the data
  INTERFACE pp_output
    MODULE PROCEDURE pp_output_1d, pp_output_2d
  END INTERFACE pp_output

  !> tag used for the number of k-point
  CHARACTER(*), PARAMETER :: tag_num_kpoint = 'NUM_KPOINT'
  !> tag used for every individual k-point
  CHARACTER(*), PARAMETER :: tag_kpoint = 'KPOINT'
  !> tag used for the k-point vector
  CHARACTER(*), PARAMETER :: tag_kpoint_vector = 'VECTOR'
  !> tag used for the number of frequencies
  CHARACTER(*), PARAMETER :: tag_num_freq = 'NUM_FREQ'
  !> tag used for every individual frequency
  CHARACTER(*), PARAMETER :: tag_freq = 'FREQUENCY'
  !> tag used for the value of the frequency
  CHARACTER(*), PARAMETER :: tag_freq_value = 'VALUE'
  !> tag used for the number of bands
  CHARACTER(*), PARAMETER :: tag_num_band = 'NUM_BAND'
  !> tag used for the band data
  CHARACTER(*), PARAMETER :: tag_band = 'BAND'

  PRIVATE pp_output_1d, pp_output_2d

CONTAINS

  !> Open all files that are use for PP printing
  !!
  !! \param nks number of k-points the data will contain
  !! \param nbnd number of bands for GW quasiparticle energies
  !! \param nw_re number of frequency points on real axis
  !! \param nw_im number of frequency point on imaginary axis
  !! \param output type that contains all files for PP
  !!
  SUBROUTINE pp_output_open_all(nks, nbnd, nw_re, nw_im, output)

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(IN) :: nw_re
    INTEGER, INTENT(IN) :: nw_im
    TYPE(output_type), INTENT(INOUT) :: output

    INTEGER dim_re, dim_im

    ! bands for band structures
    CALL pp_output_open(nks, nbnd, output%pp_dft)
    CALL pp_output_open(nks, nbnd, output%pp_gw)
    CALL pp_output_open(nks, nbnd, output%pp_vxc)
    CALL pp_output_open(nks, nbnd, output%pp_exchange)
    CALL pp_output_open(nks, nbnd, output%pp_renorm)

    ! bands * frequencies on real frequency axis
    CALL pp_output_open_xml(nks, nbnd, nw_re, output%pp_re_corr)
    CALL pp_output_open_xml(nks, nbnd, nw_re, output%pp_im_corr)
    CALL pp_output_open_xml(nks, nbnd, nw_re, output%pp_spec)

    ! bands * frequencies on imaginary axis
    CALL pp_output_open_xml(nks, nbnd, nw_im, output%pp_re_corr_iw)
    CALL pp_output_open_xml(nks, nbnd, nw_im, output%pp_im_corr_iw)
    CALL pp_output_open_xml(nks, nbnd, nw_im, output%pp_spec_iw)

  END SUBROUTINE pp_output_open_all

  !> Open a file to print the data for QE's plotband program
  !!
  !! If the filename is not set, this routine will just clear the to_file
  !! flag, so that later parts of the code can test whether data is meant
  !! to be written. If the filename is present, the file is opened and the
  !! unit is stored in the type. Then the header for the data is written
  !! into the file.
  !!
  !! \param nks number of k-points the data will contain
  !! \param nbnd number of bands the data will have
  !! \param output type that contains the filename on input and the unit
  !! and some metadata after the return of the function
  !!
  SUBROUTINE pp_output_open(nks, nbnd, output)

    USE io_files, ONLY : seqopn

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    TYPE(pp_output_type), INTENT(INOUT) :: output

    INTEGER, EXTERNAL :: find_free_unit
    LOGICAL exst

    NAMELIST /plot/ nks, nbnd

    ! if no filename is present, clear to_file flag and exit
    output%to_file = (output%filename /= '')
    IF (.NOT.output%to_file) RETURN

    ! set metadata
    output%num_band = nbnd
    output%num_freq = 1
    output%num_kpoint = nks
    output%xml_format = .FALSE.

    ! open the file
    output%iunit = find_free_unit()
    CALL seqopn(output%iunit, output%filename, "FORMATTED", exst)

    ! write namelist to file
    WRITE(output%iunit, NML=plot)

  END SUBROUTINE pp_output_open

  !> Open a file to print the data in xml format.
  !!
  !! If the filename is not set, this routine will just clear the to_file
  !! flag, so that later parts of the code can test whether data is meant
  !! to be written. If the filename is present, the file is opened and the
  !! unit is stored in the type. Then the header for the data is written
  !! into the file.
  !!
  !! \param nks number of k-points the data will contain
  !! \param nbnd number of bands the data will have
  !! \param nfreq number of frequency points the data will have
  !! \param output type that contains the filename on input and the unit
  !! and some metadata after the return of the function
  !!
  SUBROUTINE pp_output_open_xml(nks, nbnd, nfreq, output)

    USE iotk_module, ONLY: iotk_free_unit, iotk_open_write, &
                           iotk_write_begin, iotk_write_dat

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(IN) :: nfreq
    TYPE(pp_output_type), INTENT(INOUT) :: output

    ! if no filename is present, clear to_file flag and exit
    output%to_file = (output%filename /= '')
    IF (.NOT.output%to_file) RETURN

    ! set metadata
    output%num_band = nbnd
    output%num_freq = nfreq
    output%num_kpoint = nks
    output%xml_format = .FALSE.

    ! open the file
    CALL iotk_free_unit(output%iunit)
    CALL iotk_open_write(output%iunit, output%filename)

    ! write header with metadata
    CALL iotk_write_dat(output%iunit, tag_num_kpoint, nks)
    CALL iotk_write_dat(output%iunit, tag_num_band, nbnd)
    CALL iotk_write_dat(output%iunit, tag_num_freq, nfreq)

  END SUBROUTINE pp_output_open_xml

  !> Close a file associated to a particular output type
  !!
  !! Depending on the format of the file either the xml or the
  !! regular file closing routine are used.
  !!
  !! \param output file (as pp_output_type) which is closed
  !!
  SUBROUTINE pp_output_close(output)

    USE iotk_module, ONLY: iotk_write_end, iotk_close_write

    TYPE(pp_output_type), INTENT(IN) :: output

    LOGICAL opnd

    ! don't do anything if file is not open
    IF (.NOT.output%to_file) RETURN

    ! close the file (XML version)
    IF (output%xml_format) THEN
      CALL iotk_close_write(output%iunit)

    ! close the file (regular version)
    ELSE
      INQUIRE(UNIT = output%iunit, OPENED = opnd)
      IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)
      CLOSE(output%iunit)

    END IF

  END SUBROUTINE pp_output_close

  !> specialization of the interface for 1d data
  SUBROUTINE pp_output_1d(output, kpt, data)

    TYPE(pp_output_type), INTENT(IN) :: output
    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: data(:)

    LOGICAL opnd

    ! don't do anything if file is not required
    IF (.NOT.output%to_file) RETURN

    !
    ! sanity test of the input
    !
    CALL errore(__FILE__, 'data array size inconsistent (band)', output%num_band - SIZE(data))
    CALL errore(__FILE__, 'data array size inconsistent (freq)', output%num_freq - 1)
    INQUIRE(UNIT = output%iunit, OPENED = opnd)
    IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)

    !
    ! write the data to the file
    !
    WRITE(output%iunit, '(5x,3f10.6)') kpt
    WRITE(output%iunit, '(10f15.8)') data
    ! add an empty line at the end of one data set
    WRITE(output%iunit,*)

  END SUBROUTINE pp_output_1d

  !> specialization of the interface for 2d data
  SUBROUTINE pp_output_2d(output, kpt, data)

    USE iotk_module, ONLY: iotk_write_begin, iotk_write_end, iotk_write_dat

    TYPE(pp_output_type), INTENT(IN) :: output
    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: data(:,:)

    INTEGER ifreq

    ! don't do anything if file is not required
    IF (.NOT.output%to_file) RETURN

    !
    ! sanity test of the input
    !
    CALL errore(__FILE__, 'data array size inconsistent (band)', output%num_band - SIZE(data,2))
    CALL errore(__FILE__, 'data array size inconsistent (freq)', output%num_freq - SIZE(data,1))

    !
    ! write the data to the file
    !
    CALL iotk_write_begin(output%iunit, tag_kpoint)
    CALL iotk_write_dat(output%iunit, tag_kpoint_vector, kpt)
    DO ifreq = 1, SIZE(data,1)
      CALL iotk_write_begin(output%iunit, tag_freq)
!      CALL iotk_write_dat(output%iunit, tag_freq_value, 0)
      CALL iotk_write_dat(output%iunit, tag_band, data(ifreq, :))
      CALL iotk_write_end(output%iunit, tag_freq)
    END DO ! ifreq
    CALL iotk_write_end(output%iunit, tag_kpoint)

  END SUBROUTINE pp_output_2d

END MODULE pp_output_mod

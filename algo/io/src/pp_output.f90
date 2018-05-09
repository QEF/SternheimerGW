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
!> Provides routines that generate output for QE's plotband program or in xml format.
MODULE pp_output_mod

  USE kinds,       ONLY : dp
  USE gw_type_mod, ONLY : pp_output_type

  IMPLICIT NONE

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
  !> xml suffix
  CHARACTER(*), PARAMETER :: xml = '.xml'

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

    USE gw_type_mod, ONLY : name_length, output_type

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(IN) :: nw_re
    INTEGER, INTENT(IN) :: nw_im
    TYPE(output_type), INTENT(INOUT) :: output

    CHARACTER(LEN=name_length) prefix

    ! bands for band structures
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_dft)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_gw)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_vxc)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_exchange)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_renorm)

    ! prefix contains directory and file prefix
    prefix = TRIM(output%directory)//"/"//TRIM(output%prefix)//"."

    ! bands * frequencies on real frequency axis
    CALL pp_output_open_xml(nks, nbnd, nw_re, prefix, output%pp_re_corr)
    CALL pp_output_open_xml(nks, nbnd, nw_re, prefix, output%pp_im_corr)
    CALL pp_output_open_xml(nks, nbnd, nw_re, prefix, output%pp_spec)

    ! bands * frequencies on imaginary axis
    CALL pp_output_open_xml(nks, nbnd, nw_im, prefix, output%pp_re_corr_iw)
    CALL pp_output_open_xml(nks, nbnd, nw_im, prefix, output%pp_im_corr_iw)
    CALL pp_output_open_xml(nks, nbnd, nw_im, prefix, output%pp_spec_iw)

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
  !! \param directory directory in which files are created
  !! \param output type that contains the filename on input and the unit
  !! and some metadata after the return of the function
  !!
  SUBROUTINE pp_output_open(nks, nbnd, directory, output)

    USE io_files, ONLY : seqopn

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    CHARACTER(*), INTENT(IN) :: directory
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
    CALL seqopn(output%iunit, output%filename, "FORMATTED", exst, TRIM(directory)//"/")

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
  !! \param prefix directory and/or prefix added to the front of all files
  !! \param output type that contains the filename on input and the unit
  !! and some metadata after the return of the function
  !!
  SUBROUTINE pp_output_open_xml(nks, nbnd, nfreq, prefix, output)

    USE iotk_module, ONLY: iotk_free_unit, iotk_open_write, &
                           iotk_write_begin, iotk_write_dat

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(IN) :: nfreq
    CHARACTER(*), INTENT(IN) :: prefix
    TYPE(pp_output_type), INTENT(INOUT) :: output

    ! if no filename is present, clear to_file flag and exit
    output%to_file = (output%filename /= '')
    IF (.NOT.output%to_file) RETURN

    ! set metadata
    output%num_band = nbnd
    output%num_freq = nfreq
    output%num_kpoint = nks
    output%xml_format = .TRUE.

    ! open the file
    CALL iotk_free_unit(output%iunit)
    CALL iotk_open_write(output%iunit, TRIM(prefix)//TRIM(output%filename)//xml)

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

  !> Print an 1D array for QE's plotband program
  !!
  !! The array has the dimension of the number of bands. It is assumed that
  !! this routine is called nks times, where nks is the number of k-points
  !! specified when opening the file.
  !!
  !! \param output file (as pp_output_type) to which the data is printed
  !! \param kpt vector containing the current k-point
  !! \param data array containing the data
  SUBROUTINE pp_output(output, kpt, data)

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
    IF (output%xml_format) CALL errore(__FILE__, 'must not be called with xml file', 1)
    INQUIRE(UNIT = output%iunit, OPENED = opnd)
    IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)

    !
    ! write the data to the file
    !
    WRITE(output%iunit, '(5x,3f10.6)') kpt
    WRITE(output%iunit, '(10f15.8)') data
    ! add an empty line at the end of one data set
    WRITE(output%iunit,*)

  END SUBROUTINE pp_output

  !> Print a 2D array in xml format
  !!
  !! The array has the shape (nfreq, nbnd), where nfreq and nbnd are the
  !! number of frequencies and bands, respectively. It is assumed that
  !! this routine is called nks times, where nks is the number of k-points
  !! specified when opening the file.
  !!
  !! \param output file (as pp_output_type) to which the data is printed
  !! \param ikq index of the k-point
  !! \param kpt vector containing the current k-point
  !! \param omega frequency values
  !! \param data array containing the data
  SUBROUTINE pp_output_xml(output, ikq, kpt, omega, data)

    USE iotk_module, ONLY: iotk_index, iotk_write_begin, iotk_write_end, iotk_write_dat

    TYPE(pp_output_type), INTENT(IN) :: output
    INTEGER,  INTENT(IN) :: ikq
    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: omega(:)
    REAL(dp), INTENT(IN) :: data(:,:)

    INTEGER ifreq

    ! don't do anything if file is not required
    IF (.NOT.output%to_file) RETURN

    !
    ! sanity test of the input
    !
    CALL errore(__FILE__, 'data array size inconsistent (band)', output%num_band - SIZE(data,2))
    CALL errore(__FILE__, 'data array size inconsistent (freq)', output%num_freq - SIZE(data,1))
    CALL errore(__FILE__, 'frequency and data array inconsistent', SIZE(data,1) - SIZE(omega))
    IF (.NOT.output%xml_format) CALL errore(__FILE__, 'must be called with xml file', 1)

    !
    ! write the data to the file
    !
    CALL iotk_write_begin(output%iunit, tag_kpoint//TRIM(iotk_index(ikq)))
    CALL iotk_write_dat(output%iunit, tag_kpoint_vector, kpt)
    DO ifreq = 1, SIZE(data,1)
      CALL iotk_write_begin(output%iunit, tag_freq//TRIM(iotk_index(ifreq)))
      CALL iotk_write_dat(output%iunit, tag_freq_value, omega(ifreq))
      CALL iotk_write_dat(output%iunit, tag_band, data(ifreq, :))
      CALL iotk_write_end(output%iunit, tag_freq//TRIM(iotk_index(ifreq)))
    END DO ! ifreq
    CALL iotk_write_end(output%iunit, tag_kpoint//TRIM(iotk_index(ikq)))

  END SUBROUTINE pp_output_xml

END MODULE pp_output_mod

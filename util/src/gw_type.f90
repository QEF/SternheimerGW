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
!> Defines types to contain data used in GW calculations.
MODULE gw_type_mod

  !> maximal length of file name
  INTEGER, PARAMETER :: name_length = 256

  !> Type that contains the information about a file linked for PP output
  !!
  !! This type facilitates the pp_output_mod module. The type stores the
  !! filename of the file to open and once the file is opened the associated
  !! unit. To check for consistency of the data the number of k-points and
  !! the number of bands are stored, too.
  !!
  !! \param filename   name of the file to which the data is written
  !! \param to_file    flag indicating if data is print to file    
  !! \param xml_format flag indicating if the file is written in XML format
  !! \param iunit      unit access to file
  !! \param num_kpoint number of k-points
  !! \param num_band   number of bands
  !! \param num_freq   number of frequency points (if any)
  !! \see pp_output_mod
  !!
  TYPE pp_output_type
    CHARACTER(LEN=name_length) filename
    LOGICAL to_file
    LOGICAL xml_format
    INTEGER iunit
    INTEGER num_kpoint
    INTEGER num_band
    INTEGER num_freq
  END TYPE pp_output_type

  !> Type defining the output format.
  !!
  !! Defines the filename of different output that can be produced by SternheimerGW.
  !! All filenames used for PP are stored in a pp_output_type that contains
  !! some additional metadata to check for consistency of the data written to
  !! file. Filenames not set will lead to the corresponding output not being
  !! written to file.
  !!
  !! \param directory      directory in which output is written
  !! \param prefix         prefix added to all file names
  !! \param pp_dft         output of DFT Kohn-Sham eigenvalues
  !! \param pp_gw          output of GW quasi-particle eigenvalues
  !! \param pp_vxc         output of expectation values of V\f$_\text{xc}\f$
  !! \param pp_exchange    output of the exchange part of \f$\Sigma\f$
  !! \param pp_renorm      output the quasiparticle renormalization
  !! \param pp_re_corr     output of real part of correlation on real frequency axis
  !! \param pp_re_corr_iw  output of real part of correlation on imaginary frequency axis
  !! \param pp_im_corr     output of imaginary part of correlation on real frequency axis
  !! \param pp_im_corr_iw  output of imaginary part of correlation on imaginary frequency axis
  !! \param pp_spec        output of spectral function on real frequency axis
  !! \param pp_spec_iw     output of spectral function on imaginary frequency axis
  !! \param unit_sigma     unit under which file_sigma is opened
  !! \param file_sigma     filename in which Sigma_x and Sigma_c are written
  !!
  TYPE output_type
    CHARACTER(LEN=name_length) directory
    CHARACTER(LEN=name_length) prefix
    TYPE(pp_output_type) pp_dft
    TYPE(pp_output_type) pp_gw
    TYPE(pp_output_type) pp_vxc
    TYPE(pp_output_type) pp_exchange
    TYPE(pp_output_type) pp_renorm
    TYPE(pp_output_type) pp_re_corr
    TYPE(pp_output_type) pp_re_corr_iw
    TYPE(pp_output_type) pp_im_corr
    TYPE(pp_output_type) pp_im_corr_iw
    TYPE(pp_output_type) pp_spec
    TYPE(pp_output_type) pp_spec_iw

    ! file for exchange and correlation part of sigma
    INTEGER unit_sigma
    CHARACTER(LEN=name_length) file_sigma

  END TYPE output_type

END MODULE gw_type_mod

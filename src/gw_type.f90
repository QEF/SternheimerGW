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
!> Defines types to contain data used in GW calculations.
MODULE gw_type_mod

  !> maximal length of file name
  INTEGER, PARAMETER :: name_length = 256

  !> Type defining the output format.
  !!
  !! Defines the filename of different output that can be produced by
  !! SGW. Filenames not set will lead to the corresponding output not being
  !! written to file.
  !!
  !! \param file_dft         output of DFT Kohn-Sham eigenvalues
  !! \param print_dft        if DFT Kohn-Sham eigenvalues are printed
  !! \param file_gw          output of GW quasi-particle eigenvalues
  !! \param print_gw         if GW quasi-particle eigenvalues are printed
  !! \param file_vxc         output of expectation values of V\f$_\text{xc}\f$
  !! \param print_vxc        if expectation values of V\f$_\text{xc}\f$ are printed
  !! \param file_exchange    output of the exchange part of \f$\Sigma\f$
  !! \param print_exchange   if exchange part of \f$\Sigma\f$ is printed
  !! \param file_renorm      output the quasiparticle renormalization
  !! \param print_renorm     if quasiparticle renormalization is printed
  !! \param file_re_corr     output of real part of correlation on real frequency axis
  !! \param print_re_corr    if real part of correlation on real frequency axis is printed
  !! \param file_re_corr_iw  output of real part of correlation on imaginary frequency axis
  !! \param print_re_corr_iw if real part of correlation on imaginary frequency axis is printed
  !! \param file_im_corr     output of imaginary part of correlation on real frequency axis
  !! \param print_im_corr    if imaginary part of correlation on real frequency axis is printed
  !! \param file_im_corr_iw  output of imaginary part of correlation on imaginary frequency axis
  !! \param print_im_corr_iw if imaginary part of correlation on imaginary frequency axis is printed
  !! \param file_spec        output of spectral function on real frequency axis
  !! \param print_spec       if spectral function on real frequency axis is printed
  !! \param file_spec_iw     output of spectral function on imaginary frequency axis
  !! \param print_spec_iw    if spectral function on imaginary frequency axis is printed
  !!
  TYPE output_type
    CHARACTER(LEN=name_length) file_dft
    LOGICAL                    print_dft
    CHARACTER(LEN=name_length) file_gw
    LOGICAL                    print_gw
    CHARACTER(LEN=name_length) file_vxc
    LOGICAL                    print_vxc
    CHARACTER(LEN=name_length) file_exchange
    LOGICAL                    print_exchange
    CHARACTER(LEN=name_length) file_renorm
    LOGICAL                    print_renorm
    CHARACTER(LEN=name_length) file_re_corr
    LOGICAL                    print_re_corr
    CHARACTER(LEN=name_length) file_re_corr_iw
    LOGICAL                    print_re_corr_iw
    CHARACTER(LEN=name_length) file_im_corr
    LOGICAL                    print_im_corr
    CHARACTER(LEN=name_length) file_im_corr_iw
    LOGICAL                    print_im_corr_iw
    CHARACTER(LEN=name_length) file_spec
    LOGICAL                    print_spec
    CHARACTER(LEN=name_length) file_spec_iw
    LOGICAL                    print_spec_iw
  END TYPE output_type

END MODULE gw_type_mod

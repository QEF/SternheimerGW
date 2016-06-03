!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
subroutine print_clock_gw
  !-----------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  USE uspp, only: okvan
  USE control_gw, ONLY : trans, zue, epsil
  USE nlcc_gw, ONLY: nlcc_any
  implicit none
  !
  WRITE( stdout, * )
  call print_clock ('GW')
  WRITE( stdout,  * ) '    INITIALIZATION: '
  call print_clock ('gwq_setup')
  call print_clock ('gwq_init')
  WRITE( stdout, * )
  call print_clock ('gwq_init')
  if (nlcc_any) call print_clock ('set_drhoc')
  call print_clock ('init_vloc')
  call print_clock ('init_us_1')
  call print_clock ('init_us_2')
  call print_clock ('dvanqq')
  call print_clock ('drho')

  WRITE( stdout,  * ) ' G, W and sigma: '
  call print_clock ('coulomb')
  call print_clock ('green')
  call print_clock ('sigma_c')
  call print_clock ('sigma_x')
  WRITE( stdout, * )
  WRITE( stdout, * )
  call print_clock ('gwqscf')
  call print_clock ('solve_linter')
  WRITE( stdout, * )
  call print_clock ('dvqpsi_us')
  call print_clock ('ortho')
  call print_clock ('incdrhoscf')
  call print_clock ('vpsifft')
  call print_clock ('dv_of_drho')
  call print_clock ('mix_pot')
  call print_clock ('ef_shift')
  call print_clock ('localdos')
#ifdef __PARA
!  call print_clock ('psymdvscf')
#else
!  call print_clock ('symdvscf')
#endif
  WRITE( stdout, * )
  WRITE( stdout, * )
  call print_clock ('cgsolve')
  call print_clock ('cbcgsolve')
  WRITE( stdout, * )
  call print_clock ('ch_psi')
  call print_clock ('first')
  call print_clock ('h_psiq')

  call print_clock ('last')
  WRITE( stdout, * )
  call print_clock ('firstfft')
  call print_clock ('product')
  call print_clock ('secondfft')

  call print_clock ('add_vuspsi')
  WRITE( stdout, * )
  call print_clock ('incdrhoscf')

  call print_clock ('addusdbec')
  WRITE( stdout, * )
  call print_clock ('drhodvus')
  WRITE( stdout, * )
  WRITE( stdout,  * ) '     General routines'
  call print_clock ('calbec')
  call print_clock ('cft3')
  call print_clock ('cft3s')
  call print_clock ('cinterpolate')
  call print_clock ('davcio')
  call print_clock ('write_rec')
  WRITE( stdout, * )
  return
end subroutine print_clock_gw

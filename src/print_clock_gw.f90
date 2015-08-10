!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
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
  call print_clock ('greenlinsys')
  call print_clock ('sigmac')
  call print_clock ('sigma_exch')
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

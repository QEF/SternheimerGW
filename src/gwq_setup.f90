!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
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
  subroutine gwq_setup
!-----------------------------------------------------------------------
  !
  ! MODIFIED Henry Lambert
  !  This subroutine prepares several variables which are needed in the
  !  GW !  program:
  !  1) computes the total local potential (external+scf) on the smooth
  !     grid to be used in h_psi and similiar
  !  2) computes dmuxc 3) with GC if needed
  !  4) set the inverse of every matrix invs
  !  5) for metals sets the occupied bands
  !  6) computes alpha_pv
  !  9) set the variables needed to deal with nlcc
  !
  ! LEGACY FROM THE PHONON STUFF (SHOULD BE REMOVED)
  !
  !  7) computes the variables needed to pass to the pattern representation
  !     gi     the G associated to each symmetry operation
  !     gimq   the G of the q -> -q+G symmetry
  !     irgq   the small group indices
  !     nsymq  the order of the small group of q
  !     minus_q true if there is a symmetry sending q -> -q+G
  !  8) for testing purposes it sets ubar
  !  10) set the variables needed for the partial computation
  !       of the dynamical matrix
  !
  !  IMPORTANT NOTE ABOUT SYMMETRIES:
  !  nrot  is the number of sym.ops. of the Bravais lattice
  !        read from data file, only used in set_default_pw
  !  nsym  is the number of sym.ops. of the crystal symmetry group
  !        read from data file, should never be changed
  !  nsymq is the number of sym.ops. of the small group of q
  !        it is calculated in set_defaults_pw for each q
  !  The matrices "s" of sym.ops are ordered as follows:
  !   first the nsymq sym.ops. of the small group of q
  !   (the ordering is done in subroutine copy_sym in set_defaults_pw),
  !   followed by the remaining nsym-nsymq sym.ops. of the crystal group,
  !   followed by the remaining nrot-nsym sym.ops. of the Bravais  group
  !

  USE constants,          ONLY : eps14
  USE control_flags,      ONLY : noinv
  USE control_gw,         ONLY : niter_gw, alpha_mix, set_alpha_pv
  USE fft_base,           ONLY : dfftp
  USE funct,              ONLY : dmxc, dmxc_spin, dmxc_nc, dft_is_gradient, get_icorr
  USE gvecs,              ONLY : doublegrid
  USE ions_base,          ONLY : ntyp => nsp
  USE lsda_mod,           ONLY : nspin
  USE mp,                 ONLY : mp_max, mp_min
  USE nlcc_gw,            ONLY : nlcc_any
  USE noncollin_module,   ONLY : noncolin
  USE scf,                ONLY : v, vrs, vltot, kedtau
  USE spin_orb,           ONLY : domag
  USE symm_base,          ONLY : time_reversal, inverse_s
  USE uspp_param,         ONLY : upf

  implicit none

  integer :: it
  ! counters

  logical :: magnetic_sym

  call start_clock ('gwq_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Set non linear core correction variables.
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  !  3) If necessary calculate the local magnetization. This information is
  !      needed in find_sym
  !
  !IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  !IF (noncolin.and.domag) THEN
  !   DO na = 1, nat
  !      !
  !      m_loc(1,na) = starting_magnetization(ityp(na)) * &
  !                    SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
  !      m_loc(2,na) = starting_magnetization(ityp(na)) * &
  !                    SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
  !      m_loc(3,na) = starting_magnetization(ityp(na)) * &
  !                    COS( angle1(ityp(na)) )
  !   END DO
  !   ux=0.0_DP
  !   if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
  !ENDIF
  !
  ! 3) Computes the derivative of the xc potential
  !
  call setup_dmuxc()
  !
  ! Setup all gradient correction stuff
  !
  call setup_dgc()
  !
  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !
  call inverse_s()
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  call setup_nbnd_occ()
  !
  ! 6) Computes alpha_pv
  !
  if (set_alpha_pv) call setup_alpha_pv()
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym

  !set the alpha_mix parameter
  !HL probably want restart functionality here.
  do it = 2, niter_gw
     if (alpha_mix (it) < eps14) alpha_mix (it) = alpha_mix (it - 1)
  enddo

  CALL stop_clock ('gwq_setup')
  RETURN
END SUBROUTINE gwq_setup

!FXN to only include correlation energy:
!-----------------------------------------------------------------------
function dmxc_corr (rho)
  !-----------------------------------------------------------------------
  !
  !  derivative of the xc potential with respect to the local density
  !
  !
  USE kinds,         ONLY : DP
  USE funct,         ONLY : get_icorr, xc
  implicit none
  !
  real(DP), intent(in) :: rho
  ! input: the charge density ( positive )
  real(DP) :: dmxc_corr
  ! output: the derivative of the xc potential
  !
  ! local variables
  !
  real(DP) :: dr, vxp, vcp, vxm, vcm, ex, ec, rs
  real(DP), external :: dpz
  integer :: iflg
  !
  real(DP), parameter :: small = 1.E-30_DP, e2 = 2.0_DP, &
       pi34 = 0.75_DP / 3.141592653589793_DP, third = 1.0_DP /3.0_DP
  !
  dmxc_corr = 0.0_DP
  if (rho < small) then
     return
  endif
  !
  !    first case: analytical derivatives available
  !
  if (get_icorr() == 1) then
     rs = (pi34 / rho)**third
     !..exchange
     !call slater (rs, ex, vx)
     !dmxc = vx / (3.0_DP * rho)
     !..correlation
     iflg = 2
     if (rs < 1.0_DP) iflg = 1
     dmxc_corr = dpz (rs, iflg)
  else
     !
     !     second case: numerical derivatives
     !
     dr = min (1.E-6_DP, 1.E-4_DP * rho)
     call xc (rho + dr, ex, ec, vxp, vcp)
     call xc (rho - dr, ex, ec, vxm, vcm)
     !dmxc = (vxp + vcp - vxm - vcm) / (2.0_DP * dr)
     dmxc_corr = (vcp - vcm) / (2.0_DP * dr)
  endif
  !
  ! bring to rydberg units
  !
  dmxc_corr = e2 * dmxc_corr
  return
  !
end function dmxc_corr
!
!-----------------------------------------------------------------------

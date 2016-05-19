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

  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg  
  USE io_global,     ONLY : stdout
  USE ener,          ONLY : Ef
  USE klist,         ONLY : xk, lgauss, degauss, ngauss, nks, nelec
  USE ktetra,        ONLY : ltetra, tetra
  USE lsda_mod,      ONLY : nspin, lsda, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm
  USE gvecs,       ONLY : doublegrid
  USE symm_base,     ONLY : nrot, nsym, s, ftau, irt, t_rev, time_reversal, &
                            sr, invs, inverse_s
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants,     ONLY : degspin, pi
  USE noncollin_module,   ONLY : noncolin, m_loc, angle1, angle2, ux, nspin_mag
  USE wvfct,              ONLY : nbnd, et
  USE nlcc_gw,       ONLY : drc, nlcc_any
  USE eqv,           ONLY : dmuxc
  USE control_gw,    ONLY : rec_code, lgamma_gamma, search_sym, start_irr, &
                            last_irr, niter_gw, alpha_mix, all_done, &
                            epsil, lgamma, recover, where_rec, alpha_pv, &
                            nbnd_occ, flmixdpot, reduce_io, rec_code_read, &
                            done_epsil, zeu, done_zeu, current_iq, just_corr
  USE output_mod,    ONLY : fildrho
  USE lr_symm_base,  ONLY : gi, gimq, irgq, minus_q, nsymq, rtau
  USE efield_mod,    ONLY : epsilon, zstareu
  USE qpoint,        ONLY : xq
  USE partial,       ONLY : comp_irr, atomo, nat_todo, list, nrapp, all_comp, &
                            done_irr
  USE gamma_gamma,   ONLY : has_equivalent, asr, nasr, n_diff_sites, &
                            equiv_atoms, n_equiv_atoms, with_symmetry
  USE control_flags, ONLY : iverbosity, modenum, noinv
  USE disp,          ONLY : comp_irr_iq
  USE funct,         ONLY : dmxc, dmxc_spin, dmxc_nc, dft_is_gradient, get_icorr
  USE mp,            ONLY : mp_max, mp_min
  USE mp_pools,      ONLY : inter_pool_comm, npool


  implicit none

  real(DP) :: rhotot, rhoup, rhodw, target, small, fac, xmax, emin, emax
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy

  real(DP) :: sr_is(3,3,48)

  integer :: ir, isym, jsym, irot, ik, ibnd, ipol, &
       mu, nu, imode0, irr, ipert, na, it, nt, is, js, nsym_is, last_irr_eff
  ! counters

  real(DP) :: auxdmuxc(4,4)
  real(DP), allocatable :: w2(:)

  logical :: sym (48), is_symmorgwic, magnetic_sym
  integer, allocatable :: ifat(:)
  integer :: ierr

  real(DP) :: dmxc_corr

  call start_clock ('gwq_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Set non linear core correction variables.
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  if (nlcc_any) allocate (drc( ngm, ntyp))    
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
  dmuxc(:,:,:) = 0.d0
  if (lsda) then
     do ir = 1, dfftp%nnr
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        call dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     enddo
  else
     IF (noncolin.and.domag) THEN
        do ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           call dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
           DO is=1,nspin_mag
              DO js=1,nspin_mag
                 dmuxc(ir,is,js)=auxdmuxc(is,js)
              END DO
           END DO
        enddo
     ELSE
        do ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           IF(.not.just_corr) THEN
             if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
             if (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
           ELSE
           !Only correlation energy is included in the RPA + Vxc.
             if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc_corr      (rhotot)
             if (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc_corr (- rhotot)
           ENDIF
        enddo
     END IF
  endif

  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !HL-
   call inverse_s ( )
  ! 5) Computes the number of occupied bands for each k point
  if (lgauss) then
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint:
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     if (ngauss.eq. - 99) then
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     endif
     target = ef + xmax * degauss
     do ik = 1, nks
        do ibnd = 1, nbnd
           if (et (ibnd, ik) .lt.target) nbnd_occ (ik) = ibnd
        enddo
        if (nbnd_occ (ik) .eq. nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     enddo
  else if (ltetra) then
     call errore('gwq_setup','phonon + tetrahedra not implemented', 1)
  else
     if (lsda) call infomsg('gwq_setup','occupation numbers probably wrong')
     if (noncolin) then
        nbnd_occ = nint (nelec) 
     else
        do ik = 1, nks
           nbnd_occ (ik) = nint (nelec) / degspin
        enddo
     endif
  endif
  !
  ! 6) Computes alpha_pv
  !
  emin = et (1, 1)
  do ik = 1, nks
     do ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     enddo
  enddo
#ifdef __PARA
  ! find the minimum across pools
  call mp_min( emin, inter_pool_comm )
#endif
  if (lgauss) then
     emax = target
     alpha_pv = emax - emin
  else
     emax = et (1, 1)
     do ik = 1, nks
        do ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        enddo
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  endif
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym

  !set the alpha_mix parameter
  !HL probably want restart functionality here.
  do it = 2, niter_gw
     if (alpha_mix (it) .eq.0.d0) alpha_mix (it) = alpha_mix (it - 1)
  enddo

  flmixdpot = 'mixd'


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
  real(DP) :: dr, vxp, vcp, vxm, vcm, vx, ex, ec, rs
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

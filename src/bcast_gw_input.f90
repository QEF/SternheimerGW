!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_gw_input ( )
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the GW input to all
  !     the other processors
  !
  !
#ifdef __PARA
  use mp,          ONLY : mp_bcast
  use mp_world,    ONLY : world_comm
  USE io_global,   ONLY : meta_ionode_id
  USE control_gw,  ONLY : start_irr, last_irr, start_q, last_q, nmix_gw, &
                         niter_gw, lnoloc, alpha_mix, tr2_gw, lrpa, recover, &
                         ldisp, elgw, reduce_io, zue, zeu, epsil, trans, &
                         lgamma, eta, modielec, do_coulomb, do_sigma_c,& 
                         do_sigma_exx, do_sigma_exxG, do_green, do_sigma_matel, tr2_green,&
                         do_q0_only, maxter_green, godbyneeds, cohsex, padecont,&
                         multishift, do_sigma_extra, solve_direct, w_green_start, tinvert,&
                         coul_multishift, trunc_2d, do_epsil, do_serial, do_diag_g, do_diag_w,&
                         do_imag, do_pade_coul, newgrid, high_io, freq_gl, prec_direct, tmp_dir_coul, prec_shift, just_corr
  USE disp,     ONLY : iq1, iq2, iq3, nq1, nq2, nq3, kpoints, w_of_q_start
  USE partial,  ONLY : nat_todo, nrapp
  USE freq_gw,  ONLY : fpol, wsigmamin, wsigmamax, wcoulmax, deltaw, plasmon, greenzero
  USE output,   ONLY : fildvscf, fildyn, fildrho
  use io_files, ONLY : tmp_dir, prefix
  USE control_flags,    ONLY: iverbosity, modenum
  USE input_parameters, ONLY: max_seconds
  USE units_gw,         ONLY : iuncoul, iungreen, lrgrn, lrcoul, iunsigma, lrsigma, lrsex, iunsex
  USE ions_base,        ONLY : amass
  USE run_info, ONLY : title
  USE gwsigma,       ONLY : nbnd_sig, ecutsex, ecutsco, ecutprec
  USE gwsymm,        ONLY : use_symm
 
  implicit none
  call mp_bcast (trans, meta_ionode_id, world_comm )
  call mp_bcast (reduce_io, meta_ionode_id, world_comm )
  call mp_bcast (fpol, meta_ionode_id, world_comm )
  call mp_bcast (ldisp, meta_ionode_id, world_comm )
  call mp_bcast (recover, meta_ionode_id, world_comm )
  call mp_bcast (lrpa, meta_ionode_id, world_comm )
  call mp_bcast (lnoloc, meta_ionode_id, world_comm )
  !
  ! integers
  !
  call mp_bcast (start_irr, meta_ionode_id, world_comm )
  call mp_bcast (last_irr, meta_ionode_id, world_comm )
  call mp_bcast (start_q, meta_ionode_id, world_comm )
  call mp_bcast (last_q, meta_ionode_id, world_comm )
  call mp_bcast (niter_gw, meta_ionode_id, world_comm )
  call mp_bcast (nmix_gw, meta_ionode_id, world_comm )
  call mp_bcast (iverbosity, meta_ionode_id, world_comm )
  call mp_bcast (modenum, meta_ionode_id, world_comm )
  call mp_bcast (nat_todo, meta_ionode_id, world_comm )
  call mp_bcast (nrapp, meta_ionode_id, world_comm )
  CALL mp_bcast( nq1, meta_ionode_id, world_comm )
  CALL mp_bcast( nq2, meta_ionode_id, world_comm )
  CALL mp_bcast( nq3, meta_ionode_id, world_comm )
  CALL mp_bcast( iq1, meta_ionode_id, world_comm )
  CALL mp_bcast( iq2, meta_ionode_id, world_comm )
  CALL mp_bcast( iq3, meta_ionode_id, world_comm )
!
! real*8
!
  call mp_bcast (tr2_gw, meta_ionode_id, world_comm )
  call mp_bcast (tr2_green, meta_ionode_id, world_comm )
  call mp_bcast (amass, meta_ionode_id, world_comm )
  call mp_bcast (alpha_mix, meta_ionode_id, world_comm )
  call mp_bcast (max_seconds, meta_ionode_id, world_comm )
! characters
  call mp_bcast (title, meta_ionode_id, world_comm )
  call mp_bcast (fildyn, meta_ionode_id, world_comm )
  call mp_bcast (fildvscf, meta_ionode_id, world_comm )
  call mp_bcast (fildrho, meta_ionode_id, world_comm )
  call mp_bcast (tmp_dir, meta_ionode_id, world_comm )
  call mp_bcast (tmp_dir_coul, meta_ionode_id, world_comm )
  call mp_bcast (prefix, meta_ionode_id, world_comm )

 !SGW cutoffs and control
  call mp_bcast (ecutsex, meta_ionode_id, world_comm)
  call mp_bcast (ecutsco, meta_ionode_id, world_comm)
  call mp_bcast (ecutprec, meta_ionode_id, world_comm)
  call mp_bcast (nbnd_sig, meta_ionode_id, world_comm)
  call mp_bcast (modielec, meta_ionode_id, world_comm)
  call mp_bcast (godbyneeds, meta_ionode_id, world_comm)
  call mp_bcast (padecont, meta_ionode_id, world_comm)
  call mp_bcast (multishift, meta_ionode_id, world_comm)
  call mp_bcast (cohsex, meta_ionode_id, world_comm)
  call mp_bcast (eta, meta_ionode_id, world_comm)
  call mp_bcast (kpoints, meta_ionode_id, world_comm)
  call mp_bcast (do_coulomb, meta_ionode_id, world_comm)
  call mp_bcast (do_sigma_c, meta_ionode_id, world_comm)
  call mp_bcast (do_sigma_exx, meta_ionode_id, world_comm)
  call mp_bcast (do_sigma_exxG, meta_ionode_id, world_comm)
  call mp_bcast (do_green, meta_ionode_id, world_comm)
  call mp_bcast (do_sigma_matel, meta_ionode_id, world_comm)
  call mp_bcast (do_q0_only, meta_ionode_id, world_comm)
  call mp_bcast (do_sigma_extra, meta_ionode_id, world_comm)
  call mp_bcast (plasmon,meta_ionode_id, world_comm)
  call mp_bcast (greenzero,meta_ionode_id, world_comm)
  call mp_bcast (solve_direct,meta_ionode_id, world_comm)
  call mp_bcast (tinvert,meta_ionode_id, world_comm)
  call mp_bcast (coul_multishift,meta_ionode_id, world_comm)
  call mp_bcast (trunc_2d,meta_ionode_id, world_comm)
  call mp_bcast (do_epsil,meta_ionode_id, world_comm)
  call mp_bcast (do_serial, meta_ionode_id, world_comm)
  call mp_bcast (do_diag_w, meta_ionode_id, world_comm)
  call mp_bcast (do_diag_g, meta_ionode_id, world_comm)
  call mp_bcast (do_imag, meta_ionode_id, world_comm)
  call mp_bcast (do_pade_coul, meta_ionode_id, world_comm)

!Frequency grid
  call mp_bcast (wsigmamin, meta_ionode_id, world_comm)
  call mp_bcast (wsigmamax, meta_ionode_id, world_comm)
  call mp_bcast (wcoulmax,  meta_ionode_id, world_comm)
  call mp_bcast (deltaw,    meta_ionode_id, world_comm)
  call mp_bcast (just_corr,    meta_ionode_id, world_comm)

!units information
!  call mp_bcast (iuncoul,    meta_ionode_id, world_comm)
!  call mp_bcast (iungreen,    meta_ionode_id, world_comm)
!  call mp_bcast (iunsigma,    meta_ionode_id, world_comm)
  call mp_bcast (lrcoul,    meta_ionode_id, world_comm)
  call mp_bcast (lrgrn,    meta_ionode_id, world_comm)
  call mp_bcast (high_io,    meta_ionode_id, world_comm)
  call mp_bcast (freq_gl,    meta_ionode_id, world_comm)
  call mp_bcast (prec_direct,    meta_ionode_id, world_comm)
  call mp_bcast (prec_shift,    meta_ionode_id, world_comm)

  call mp_bcast (use_symm, meta_ionode_id, world_comm)
  call mp_bcast (w_of_q_start, meta_ionode_id, world_comm)
  call mp_bcast (w_green_start, meta_ionode_id, world_comm)
  call mp_bcast (maxter_green, meta_ionode_id, world_comm)
  call mp_bcast (newgrid, meta_ionode_id, world_comm)
#endif
  return
end subroutine bcast_gw_input

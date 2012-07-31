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

  use mp, only: mp_bcast
  USE mp_global, only : intra_image_comm
  USE control_gw, ONLY : start_irr, last_irr, start_q, last_q, nmix_gw, &
                         niter_gw, lnoloc, alpha_mix, tr2_gw, lrpa, recover, &
                         ldisp, elgw, reduce_io, zue, zeu, epsil, trans, &
                         lgamma, eta, modielec, do_coulomb, do_sigma_c,& 
                         do_sigma_exx, do_green, do_sigma_matel, tr2_green,&
                         do_q0_only, maxter_green, godbyneeds, cohsex, padecont,&
                         multishift
                         !do_sigma_extra
 
  USE gamma_gamma, ONLY : asr
  USE disp, ONLY : iq1, iq2, iq3, nq1, nq2, nq3, kpoints, w_of_q_start
  USE partial, ONLY : nat_todo, nrapp
 !HL new freq variable
  USE freq_gw, ONLY : fpol, wsigmamin, wsigmamax, wcoulmax, deltaw
  USE output, ONLY : fildvscf, fildyn, fildrho
  use io_files, ONLY : tmp_dir, prefix
  USE control_flags, only: iverbosity, modenum
  USE input_parameters, ONLY: max_seconds
  USE units_gw,         ONLY : iuncoul, iungreen, lrgrn, lrcoul, iunsigma, lrsigma, lrsex, iunsex
  USE ions_base,     ONLY : amass
  USE io_global, ONLY : ionode_id
  USE printout_base, ONLY : title
  USE gwsigma,       ONLY : nbnd_sig, ecutsig, ecutsex, ecutsco, ecutgrn, ecutpol
  USE gwsymm,        ONLY : use_symm
 
  implicit none
 !logicals
 !HL lgamma, epsi
 !call mp_bcast (lgamma, ionode_id )
 !call mp_bcast (epsil, ionode_id )
  call mp_bcast (trans, ionode_id )
  call mp_bcast (reduce_io, ionode_id )
  call mp_bcast (fpol, ionode_id )
  call mp_bcast (ldisp, ionode_id )
  call mp_bcast (recover, ionode_id )
  call mp_bcast (asr, ionode_id )
  call mp_bcast (lrpa, ionode_id )
  call mp_bcast (lnoloc, ionode_id )
  !
  ! integers
  !
  call mp_bcast (start_irr, ionode_id )
  call mp_bcast (last_irr, ionode_id )
  call mp_bcast (start_q, ionode_id )
  call mp_bcast (last_q, ionode_id )
  call mp_bcast (niter_gw, ionode_id )
  call mp_bcast (nmix_gw, ionode_id )
  call mp_bcast (iverbosity, ionode_id )
  call mp_bcast (modenum, ionode_id )
  call mp_bcast (nat_todo, ionode_id )
  call mp_bcast (nrapp, ionode_id )
  CALL mp_bcast( nq1, ionode_id )
  CALL mp_bcast( nq2, ionode_id )
  CALL mp_bcast( nq3, ionode_id )
  CALL mp_bcast( iq1, ionode_id )
  CALL mp_bcast( iq2, ionode_id )
  CALL mp_bcast( iq3, ionode_id )
!
! real*8
!
  call mp_bcast (tr2_gw, ionode_id )
  call mp_bcast (tr2_green, ionode_id )
  call mp_bcast (amass, ionode_id )
  call mp_bcast (alpha_mix, ionode_id )
  call mp_bcast (max_seconds, ionode_id )
! characters
  call mp_bcast (title, ionode_id )
  call mp_bcast (fildyn, ionode_id )
  call mp_bcast (fildvscf, ionode_id )
  call mp_bcast (fildrho, ionode_id )
  call mp_bcast (tmp_dir, ionode_id )
  call mp_bcast (prefix, ionode_id )

 !SGW cutoffs and control
  call mp_bcast (ecutsex, ionode_id)
  call mp_bcast (ecutsco, ionode_id)
  call mp_bcast (ecutsig, ionode_id)
  call mp_bcast (ecutpol, ionode_id)
  call mp_bcast (ecutgrn, ionode_id)
  call mp_bcast (nbnd_sig, ionode_id)
  call mp_bcast (modielec, ionode_id)
  call mp_bcast (godbyneeds, ionode_id)
  call mp_bcast (padecont, ionode_id)
  call mp_bcast (multishift, ionode_id)
  call mp_bcast (cohsex, ionode_id)
  call mp_bcast (eta, ionode_id)
  call mp_bcast (kpoints, ionode_id)
  call mp_bcast (do_coulomb, ionode_id)
  call mp_bcast (do_sigma_c, ionode_id)
  call mp_bcast (do_sigma_exx, ionode_id)
  call mp_bcast (do_green, ionode_id)
  call mp_bcast (do_sigma_matel, ionode_id)
  call mp_bcast (do_q0_only, ionode_id)
!HLS  call mp_bcast (do_sigma_extra, ionode_id)

!Frequency grid
  call mp_bcast (wsigmamin, ionode_id)
  call mp_bcast (wsigmamax, ionode_id)
  call mp_bcast (wcoulmax,  ionode_id)
  call mp_bcast (deltaw,    ionode_id)

!units information
  call mp_bcast (iuncoul,    ionode_id)
  call mp_bcast (lrcoul,    ionode_id)
  call mp_bcast (lrgrn,    ionode_id)
  call mp_bcast (iungreen,    ionode_id)

  call mp_bcast (use_symm, ionode_id)
  call mp_bcast (w_of_q_start, ionode_id)
  call mp_bcast (maxter_green, ionode_id)
#endif
  return
end subroutine bcast_gw_input

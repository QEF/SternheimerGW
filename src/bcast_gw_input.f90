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
                         lgamma, eta
 
   !HL elph out of control_ph

  USE gamma_gamma, ONLY : asr
  USE disp, ONLY : iq1, iq2, iq3, nq1, nq2, nq3
  USE partial, ONLY : nat_todo, nrapp
 !HL new freq variable
  USE freq_gw, ONLY : fpol
  USE output, ONLY : fildvscf, fildyn, fildrho
  use io_files, ONLY : tmp_dir, prefix
  USE control_flags, only: iverbosity, modenum
  USE input_parameters, ONLY: max_seconds
  USE ions_base,     ONLY : amass
  USE io_global, ONLY : ionode_id
  USE printout_base, ONLY : title
  USE gwsigma,       ONLY : nbnd_sig, ecutsig

  implicit none
  !
  ! logicals
  !
  !HL lgamma, epsi
  !  call mp_bcast (lgamma, ionode_id )
  !  call mp_bcast (epsil, ionode_id )
  call mp_bcast (trans, ionode_id )

  !HL 
  !call mp_bcast (zue, ionode_id )
  !call mp_bcast (zeu, ionode_id )
  !call mp_bcast (elph, ionode_id )
  !call mp_bcast (lraman, ionode_id )
  !call mp_bcast (elop, ionode_id )
  
  call mp_bcast (reduce_io, ionode_id )
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
  call mp_bcast (amass, ionode_id )
  call mp_bcast (alpha_mix, ionode_id )
  call mp_bcast (max_seconds, ionode_id )

  ! HL
  ! call mp_bcast (eth_rps, ionode_id )
  ! call mp_bcast (eth_ns, ionode_id )
  ! call mp_bcast (dek, ionode_id )
  !

  ! characters
  !
  call mp_bcast (title, ionode_id )
  call mp_bcast (fildyn, ionode_id )
  call mp_bcast (fildvscf, ionode_id )
  call mp_bcast (fildrho, ionode_id )
  call mp_bcast (tmp_dir, ionode_id )
  call mp_bcast (prefix, ionode_id )

 !SGW specific
  call mp_bcast (ecutsig, ionode_id)
  call mp_bcast (nbnd_sig, ionode_id)
  call mp_bcast (modielec, ionode_id)
  call mp_bcast (eta, ionode_id)



#endif
  return
end subroutine bcast_gw_input

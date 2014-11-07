!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------
subroutine deallocate_gwq
!----------========-----------------------
!
!  deallocates the variables allocated by allocate_gwq
!
  USE noncollin_module,      ONLY: m_loc
  USE becmod,                ONLY: bec_type, becp, deallocate_bec_type
  USE wavefunctions_module,  ONLY: evc
  USE qpoint,                ONLY: eigqts, igkq, ikks, ikqs, nksq
  USE gc_gw,                 ONLY: grho, gmag, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s, &
                                   vsgga, segni
  USE gamma_gamma,           ONLY : with_symmetry, has_equivalent, equiv_atoms, &
                                    n_equiv_atoms
  USE eqv,                   ONLY : dmuxc, vlocq, dpsi, dvpsi, evq, eprec
  USE nlcc_gw,               ONLY : drc
  USE units_gw,              ONLY : this_dvkb3_is_on_file, this_pcxpsi_is_on_file
  USE control_gw,            ONLY : lgamma
  USE gwsigma,               ONLY: scrcoul, green
  USE eqv,                   ONLY: dpsim, dpsip, dvbare

  IMPLICIT NONE
  INTEGER :: ik, ipol

  if(allocated(dvpsi)) deallocate (dvpsi)    
  if(allocated(dpsi)) deallocate ( dpsi)    
  if(allocated(vlocq)) deallocate (vlocq)
  if(allocated(dmuxc)) deallocate (dmuxc)
  if(allocated(eprec)) deallocate (eprec)
  if(allocated(ikks)) deallocate (ikks)
  if(allocated(ikqs)) deallocate (ikqs)
  if(allocated(eigqts)) deallocate (eigqts)
  !HL \psi\pm W,G, dvbare 
  if(allocated(dpsim)) deallocate(dpsim)
  if(allocated(dpsip)) deallocate(dpsip)
  if(allocated(dvbare)) deallocate(dvbare)
  if(allocated(this_dvkb3_is_on_file)) deallocate (this_dvkb3_is_on_file)    
  if(allocated(this_pcxpsi_is_on_file)) deallocate (this_pcxpsi_is_on_file)
  return
end subroutine deallocate_gwq

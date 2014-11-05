!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_gwq
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !

  USE kinds, ONLY : DP
  USE klist, ONLY : nks, nkstot
  USE wvfct, ONLY : nbnd, igk, npwx
  USE gvect, ONLY : ngm
  USE lsda_mod, ONLY : nspin
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wavefunctions_module,  ONLY: evc
  USE spin_orb,       ONLY : lspinorb
  USE becmod,         ONLY: bec_type, becp, allocate_bec_type
  USE uspp,           ONLY: okvan, nkb
  USE paw_variables,  ONLY : okpaw
  USE uspp_param,     ONLY: nhm
  USE freq_gw,        ONLY : fpol, fiu, nfs, nfsmax
  USE gwsigma,        ONLY: scrcoul, green, sigma
  USE ions_base,      ONLY : nat, ntyp => nsp

  USE qpoint,        ONLY : nksq, eigqts, igkq
  USE gwus,          ONLY : becsumort, dpqq, dpqq_so, alphasum, alphasum_nc,& 
                            becsum_nc, becp1, alphap
  USE eqv,        ONLY : dpsi, evq, vlocq, dmuxc, dvpsi, eprec, dvbare, dpsim, dpsip
  USE units_gw,   ONLY : this_pcxpsi_is_on_file, this_dvkb3_is_on_file
  USE control_gw, ONLY : lgamma  
  !USE gsmooth,    ONLY : nrxxs 
  USE fft_base,    ONLY : dfftp
  USE disp,       ONLY : gmap, eval_occ


  implicit none
  INTEGER :: ik, ipol
  !
  !   FOR LGAMMA
  if (lgamma) then
     !
     !  q=0  : evq and igkq are pointers to evc and igk
     !
     evq  => evc
     igkq => igk
  else
 !
 !q!=0 : evq, igkq are allocated and calculated at point k+q
 !
    allocate (evq ( npwx*npol , nbnd))    
    allocate (igkq ( npwx))    
  endif

  allocate (dvpsi ( npwx*npol , nbnd))    
  allocate (vlocq ( ngm , ntyp))    
  allocate (eprec ( nbnd, nksq) )
  allocate (eigqts ( nat))
!  allocate (dmuxc (nrxx , nspin_mag , nspin_mag))    
!  allocate (dvbare(nrxxs))    
!HL
  allocate (dmuxc (dfftp%nnr , nspin_mag , nspin_mag))    
  allocate (dvbare(dfftp%nnr))    

   if (okvan) then
     allocate (alphasum ( nhm * (nhm + 1)/2 , 3 , nat , nspin_mag)) 
     allocate (this_pcxpsi_is_on_file(nksq,3))
     this_pcxpsi_is_on_file(:,:)=.false.
   endif
  ALLOCATE (becp1(nksq))

  DO ik=1,nksq
     call allocate_bec_type ( nkb, nbnd, becp1(ik) )
  END DO
  CALL allocate_bec_type ( nkb, nbnd, becp )
  return
end subroutine allocate_gwq

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

  USE kinds, only : DP
  USE klist, only : nks, nkstot
  USE wvfct, ONLY : nbnd, igk, npwx
  USE gvect, ONLY : nrxx, ngm
  USE lsda_mod, ONLY : nspin
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wavefunctions_module,  ONLY: evc
  USE spin_orb, ONLY : lspinorb
  USE becmod, ONLY: bec_type, becp, allocate_bec_type
  USE uspp, ONLY: okvan, nkb
  USE paw_variables, ONLY : okpaw
  USE uspp_param, ONLY: nhm
  USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax

 
  USE gwsigma, ONLY: scrcoul, green, sigma
  

  USE ions_base, ONLY : nat, ntyp => nsp

  USE qpoint, ONLY : nksq, eigqts, igkq
  USE gwus, ONLY : int1, int1_nc, int2, int2_so, int3, int3_nc, int3_paw, &
                   int4, int4_nc, int5, int5_so, becsumort, dpqq, &
                   dpqq_so, alphasum, alphasum_nc, becsum_nc, &
                   becp1, alphap
  USE efield_mod, ONLY : zstareu, zstareu0, zstarue0, zstarue0_rec, zstarue
  USE eqv, ONLY : dpsi, evq, vlocq, dmuxc, dvpsi, eprec, dvbare, dpsim, dpsip
  USE units_gw, ONLY : this_pcxpsi_is_on_file, this_dvkb3_is_on_file
  USE dynmat, ONLY : dyn00, dyn, dyn_rec, w2
  USE modes, ONLY : u, ubar, rtau, npert, name_rap_mode
  USE control_gw, ONLY : lgamma  
  USE gsmooth, ONLY : nrxxs 
  USE disp,    ONLY : gmap, eval_occ


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
  !
  allocate (dvpsi ( npwx*npol , nbnd))    

!  if(.not.coul_multi) then
!  allocate (dpsi  ( npwx*npol , nbnd))    
!  allocate (dpsim ( npwx*npol , nbnd))    
!  allocate (dpsip ( npwx*npol , nbnd))    
!  endif

  !
  allocate (vlocq ( ngm , ntyp))    
  allocate (eprec ( nbnd, nksq) )
  allocate (eigqts ( nat))
  allocate (dmuxc (nrxx , nspin_mag , nspin_mag))    
  allocate (dvbare(nrxxs))    

   if (okvan) then
  !int1, int2, etc are all defined in PRB 64 235118
     allocate (int1 ( nhm, nhm, 3, nat, nspin_mag))
     allocate (int2 ( nhm , nhm , 3 , nat , nat))
     if (okpaw) then
        allocate (becsumort ( nhm*(nhm+1)/2 , nat , nspin, 3*nat))
     endif

     allocate (int4 ( nhm * (nhm + 1)/2,  3 , 3 , nat, nspin_mag))
     allocate (int5 ( nhm * (nhm + 1)/2 , 3 , 3 , nat , nat))
     allocate (dpqq( nhm, nhm, 3, ntyp))
     IF (noncolin) THEN
        ALLOCATE(int1_nc( nhm, nhm, 3, nat, nspin))
        ALLOCATE(int4_nc( nhm, nhm, 3, 3, nat, nspin))
        ALLOCATE(becsum_nc( nhm*(nhm+1)/2, nat, npol, npol))
        ALLOCATE(alphasum_nc( nhm*(nhm+1)/2, 3, nat, npol, npol))
        IF (lspinorb) THEN
           ALLOCATE(int2_so( nhm, nhm, 3, nat , nat, nspin))
           ALLOCATE(int5_so( nhm, nhm, 3, 3, nat , nat, nspin))
           allocate(dpqq_so( nhm, nhm, nspin, 3, ntyp))
        END IF
     END IF
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

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
subroutine allocate_gwq
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
  USE becmod,                ONLY : bec_type, becp, allocate_bec_type
  USE control_lr,            ONLY : lgamma
  USE eqv,                   ONLY : evq, vlocq, dmuxc, dvpsi
  USE eqv_gw,                ONLY : dvbare
  USE fft_base,              ONLY : dfftp, dffts
  USE gvect,                 ONLY : ngm
  USE ions_base,             ONLY : nat, ntyp => nsp
  USE lrus,                  ONLY : becp1
  USE noncollin_module,      ONLY : npol, nspin_mag
  USE qpoint,                ONLY : nksq, eigqts, igkq
  USE uspp,                  ONLY : nkb
  USE wvfct,                 ONLY : nbnd, npwx
  USE wavefunctions_module,  ONLY : evc

  IMPLICIT NONE

  INTEGER :: ik
  !
  !   FOR LGAMMA
  if (lgamma) then
    !
    ! q=0  : evq is a pointer to evc
    !
    evq  => evc
  else
    !
    ! q!=0 : evq is and calculated at point k+q
    !
    allocate (evq ( npwx*npol , nbnd))    
    allocate (igkq ( npwx))    
  endif

  allocate (dvpsi ( npwx*npol , nbnd))    
  allocate (vlocq ( ngm , ntyp))    
  allocate (eigqts ( nat))
  allocate (dmuxc (dfftp%nnr , nspin_mag , nspin_mag))    
  allocate (dvbare(dffts%nnr))    

  ALLOCATE (becp1(nksq))
  DO ik=1,nksq
     call allocate_bec_type ( nkb, nbnd, becp1(ik) )
  END DO
  CALL allocate_bec_type ( nkb, nbnd, becp )

  return
end subroutine allocate_gwq

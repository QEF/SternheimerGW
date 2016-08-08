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
subroutine incdrhoscf_w (drhoscf, weight, ik, dbecsum, dpsi)
!-----------------------------------------------------------------------
!
!     This routine computes the change of the charge density due to the
!     perturbation. It is called at the end of the computation of the
!     change of the wavefunction for a given k point.
!

  USE kinds, only : DP
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat
  USE wvfct,     ONLY : npwx, nbnd
  USE uspp_param,ONLY: nhm
  USE wavefunctions_module,  ONLY: evc
  USE klist,     ONLY : ngk, igk_k
  USE qpoint,    ONLY : ikks, ikqs
  USE control_gw, ONLY : nbnd_occ
  !USE gsmooth,   ONLY : dffts%nnr, nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  USE gvecs,     ONLY : nls
  USE fft_base,  ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : invfft
  USE mp_bands,   ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp, ONLY : mp_sum


  implicit none

  integer :: ik
  ! input: the k point

  real(DP) :: weight
  ! input: the weight of the k point
  complex(DP) :: drhoscf (dffts%nnr), dbecsum (nhm*(nhm+1)/2,nat)
  complex(DP) :: dpsi(npwx, nbnd) 
  ! output: the change of the charge densit
  ! inp/out: the accumulated dbec
  ! here the local variable

  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable  :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: npw, npwq
  integer :: ibnd, ikk, ikq, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  allocate (dpsic(  dffts%nnr))    
  allocate (psi  (  dffts%nnr))    
  wgt = 2.d0 * weight / omega
 !HLTIL
 !wgt = weight / omega

  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  ! dpsi contains the perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_occ (ikk)
     psi (:) = (0.d0, 0.d0)
     do ig = 1, npw
        psi (nls (igk_k (ig, ikk) ) ) = evc (ig, ibnd)
     enddo

     CALL invfft('Wave', psi, dffts)

     dpsic(:) = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsic (nls (igk_k (ig, ikq) ) ) = dpsi (ig, ibnd)
     enddo

     CALL invfft('Wave', dpsic, dffts)  
  
     do ir = 1, dffts%nnr
        drhoscf (ir) = drhoscf (ir) + wgt * (CONJG(psi (ir) ) * dpsic (ir))
     enddo
  enddo
! HL adds B15 of DalCorso. REQUIRED FOR ULTRASOFT
!  call addusdbec (ik, weight, dpsi, dbecsum)   
  deallocate (psi)
  deallocate (dpsic)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf_w


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
subroutine dvqpsi_us (dvbarein, ik, addnlcc)
!----------------------------------------------------------------------

  !
  !Initializes the electric field potential the perturbing potential generated in Coulomb.
  !

  USE eqv,                  ONLY: dvpsi
  USE fft_base,             ONLY: dfftp, dffts
  USE fft_interfaces,       ONLY: invfft, fwfft
  USE gvecs,                ONLY: doublegrid
  USE kinds,                ONLY: dp
  USE klist,                ONLY: igk_k
  USE nlcc_gw,              ONLY: nlcc_any
  USE noncollin_module,     ONLY: npol
  USE qpoint,               ONLY: npwq, ikks, ikqs
  USE wavefunctions_module, ONLY: evc
  USE wvfct,                ONLY: nbnd, npw, npwx

  implicit none
  !
  !   The dummy variables
  integer :: ik
  ! input: the k point
  ! counter for G-vector for plane wave perturbation  
  logical :: addnlcc
  !
  !   And the local variables
  !
  integer :: ikk, ikq, ig, ibnd, ir, ip
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh
  !
  complex(DP) , allocatable, target :: aux (:)
  complex(DP) , allocatable :: aux1 (:), aux2 (:)
  complex(DP) , pointer :: auxs (:)

  complex(DP) dvbarein(dffts%nnr)

  ! work space

  call start_clock ('dvqpsi_us')
  if (nlcc_any.and.addnlcc) then
     allocate (aux( dfftp%nnr))    
     if (doublegrid) then
        allocate (auxs( dffts%nnr))    
     else
        auxs => aux
     endif
  endif

  allocate (aux1( dffts%nnr))    
  allocate (aux2( dffts%nnr))    
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  dvpsi(:,:) = (0.d0, 0.d0)
  aux1(:) = (0.d0, 0.d0)

  !HL INITIALIZE PERTURBATION dvbarein is in real space.

   aux1 (:) = dvbarein(:)

  !HL- Compute Delta V_q(r')*Psi_nk(r') 
  do ibnd = 1, nbnd
     do ip=1,npol
        aux2(:) = (0.d0, 0.d0)
        if (ip==1) then
           do ig = 1, npw
              aux2(dffts%nl(igk_k(ig, ikk))) = evc(ig, ibnd)
           enddo
        else
           do ig = 1, npw
              aux2(dffts%nl(igk_k(ig, ikk))) = evc(ig + npwx, ibnd)
           enddo
        end if
        !
        !  This wavefunction is transformed into real space
        !
        !call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
        CALL invfft('Rho', aux2, dffts)
        do ir = 1, dffts%nnr
           aux2 (ir) = aux2 (ir) * aux1 (ir)
        enddo
        !call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
        CALL fwfft('Rho', aux2, dffts)

        if (ip==1) then
           do ig = 1, npwq
              dvpsi(ig, ibnd) = aux2(dffts%nl(igk_k(ig, ikq)))
           enddo
        else
           do ig = 1, npwq
              dvpsi(ig+npwx, ibnd) = aux2(dffts%nl(igk_k(ig, ikq)))
           enddo
        end if
     enddo
  enddo

  deallocate (aux2)
  deallocate (aux1)
  if (nlcc_any.and.addnlcc) then
     deallocate (aux)
     if (doublegrid) deallocate (auxs)
  endif

  call stop_clock ('dvqpsi_us')
  return
end subroutine dvqpsi_us

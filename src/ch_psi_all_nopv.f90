!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ch_psi_all_nopv (n, h, ah, e, cw, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE kinds, only : DP
  USE becmod, ONLY : becp, calbec
  USE uspp, ONLY: nkb, vkb
  USE wvfct, ONLY : npwx, nbnd
  USE control_gw, ONLY : alpha_pv, nbnd_occ
  USE eqv,  ONLY : evq
  USE qpoint, ONLY : ikqs
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE mp_bands,             ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum

  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  complex(kind=DP) :: e (m), cw

  ! input: the eigenvalue + iw
  ! HL or just w + i eta for green linear system

  complex(kind=DP) :: h (npwx, m), ah (npwx, m)
  !  input: the vector
  !  output: the operator applied to the vector
  !  local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(kind=DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  call start_clock ('ch_psi')
  allocate (ps  ( nbnd , m))    
  allocate (hpsi( npwx , m))    
  allocate (spsi( npwx , m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  call h_psiq (npwx, n, m, h, hpsi, spsi)

  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah =(0.0d0, 0.0d0)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd)= hpsi(ig, ibnd) - (e(ibnd) + cw)*spsi(ig, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     do ibnd = 1, m
        do ig = 1, n
           ah (ig+npwx,ibnd)=hpsi(ig+npwx,ibnd) - (e(ibnd) + cw)*spsi(ig+npwx,ibnd)
        end do
     end do
  END IF
!  Here we compute the projector in the valence band
!  ikq = ikqs(ik)
!  ps (:,:) = (0.d0, 0.d0)
!  call zgemm ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evq, &
!       npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
!  HL need to remember the projector should be double complex
!  ps (:,:) = ps(:,:) * dcmplx(alpha_pv, 0.0d0)
!
!  hpsi (:,:) = (0.d0, 0.d0)
!  call zgemm ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
!       npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
!  spsi(:,:) = hpsi(:,:)
!  And apply S again
!  call calbec (n, vkb, hpsi, becp, m)
!  call s_psi (npwx, n, m, hpsi, spsi)
!  do ibnd = 1, m
!     do ig = 1, n
!        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
!     enddo
!  enddo
!  IF (noncolin) THEN
!     do ibnd = 1, m
!        do ig = 1, n
!           ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
!        enddo
!     enddo
!  ENDIF
  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all_nopv

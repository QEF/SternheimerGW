!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ch_psi_all_green (n, h, ah, e, cw, ik, m)
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

  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag


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
  
  logical :: use_projector

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
  !   then we compute the operator H-epsilon S
  !
  ah =(0.0d0, 0.0d0)
  !write(6,*) cw
  do ibnd = 1, m
     do ig = 1, n
       !ah (ig, ibnd)= hpsi(ig, ibnd)-cw*spsi(ig, ibnd)*h(ig,ibnd)
       !ah (ig, ibnd)= hpsi(ig, ibnd) - cw*(h(ig,ibnd))
        ah (ig, ibnd)= hpsi(ig, ibnd) - cw*(spsi(ig,ibnd))
     enddo
  enddo

  IF (noncolin) THEN
     do ibnd = 1, m
        do ig = 1, n
          ah (ig+npwx,ibnd)=hpsi(ig+npwx,ibnd)-e(ibnd)*spsi(ig+npwx,ibnd) - cw*(h(ig+npwx,ibnd)) 
        end do
     end do
  END IF

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all_green

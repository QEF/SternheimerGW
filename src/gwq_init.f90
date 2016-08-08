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
SUBROUTINE gwq_init()
!----------------------------------------------------------------------------
  !
  !     0) initializes the structure factors
  !     2) Computes the non-linear core correction in cases where that is
  !     required.
  !     3) Computes the preconditioner for the linear system solver.
  !
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : bg, tpiba, tpiba2, omega
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE becmod,               ONLY : calbec
  USE constants,            ONLY : eps8, tpi
  USE gvect,                ONLY : g, ngm
  USE klist,                ONLY : xk, nkstot, ngk
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE buffers,              ONLY : get_buffer
  USE io_global,            ONLY : stdout, meta_ionode_id
  USE atom,                 ONLY : msh, rgrid
  USE vlocal,               ONLY : strf
  USE spin_orb,             ONLY : lspinorb
  USE wvfct,                ONLY : g2kin, npwx, npw, nbnd
  USE gvecw,                ONLY : ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf
  USE eqv_gw,               ONLY : vlocq, evq, eprectot
  USE nlcc_gw,              ONLY : nlcc_any
  USE control_gw,           ONLY : nbnd_occ, lgamma
  USE units_gw,             ONLY : lrwfc, iuwfc
  USE qpoint,               ONLY : xq, npwq, nksq, eigqts, ikks, ikqs
  USE mp_bands,            ONLY : intra_bgrp_comm
  USE mp,                  ONLY : mp_sum
  USE mp_pools,            ONLY : inter_pool_comm, my_pool_id, npool, kunit
  USE mp_world,             ONLY : nproc, mpime
  !
  IMPLICIT NONE
  !
  ! ... local variables

  INTEGER :: nt, ik, ikq, ipol, ibnd, ikk, na, ig, irr, imode0
    ! counter on atom types
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on bands
    ! index for wavefunctions at k
    ! counter on atoms
    ! counter on G vectors
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)
    ! used to compute alphap
  COMPLEX(DP), EXTERNAL :: zdotc
  INTEGER :: nbase, nks, rest
  !
  !
  !
  CALL start_clock( 'gwq_init' )
  !
  ALLOCATE( aux1( npwx*npol, nbnd ) )    

!HL structure factors
  DO na = 1, nat
     arg = ( xq(1) * tau(1,na) + &
             xq(2) * tau(2,na) + &
             xq(3) * tau(3,na) ) * tpi
     eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) ,kind=DP)
  END DO
  !                 
  !
  ! ... a0) compute rhocore for each atomic-type if needed for nlcc
  !
  nks    = kunit * ( nkstot / kunit / npool )
  rest = ( nkstot - nks * npool ) / kunit
  IF ( ( my_pool_id + 1 ) <= rest ) nks = nks + kunit
  nbase = nks * my_pool_id
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit

  eprectot(:,:) = 0.0d0
  DO ik = 1, nksq
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npw  = ngk(ikk)
     npwq = ngk(ikq)
     !
     IF ( lsda ) current_spin = isk( ikk )
     !
     ! ... if there is only one k-point evc, evq, npw, igk stay in memory
     !
     IF ( .NOT. lgamma ) THEN
     !
        IF ( ABS( xq(1) - ( xk(1,ikq) - xk(1,ikk) ) ) > eps8 .OR. &
             ABS( xq(2) - ( xk(2,ikq) - xk(2,ikk) ) ) > eps8 .OR. &
             ABS( xq(3) - ( xk(3,ikq) - xk(3,ikk) ) ) > eps8 ) THEN
           WRITE( stdout,'(/,5x,"k points #",i6," and ", &
                  & i6,5x," total number ",i6)') ikk, ikq, nksq
           WRITE( stdout, '(  5x,"Expected q ",3f10.7)')(xq(ipol), ipol=1,3)
           WRITE( stdout, '(  5x,"Found      ",3f10.7)')((xk(ipol,ikq) &
                                                -xk(ipol,ikk)), ipol = 1, 3)
           CALL errore( 'gwq_init', 'wrong order of k points', 1 )
        END IF
        !
     END IF
     !
     ! ... read the wavefunctions at k+1
     !
     CALL get_buffer( evq, lrwfc, iuwfc, ikq )
     !
     !
     ! diagonal elements of the unperturbed Hamiltonian,
     ! needed for preconditioning
     !
     CALL g2_kin(ikq)
     !
     aux1 = (0.d0,0.d0)
     DO ig = 1, npwq
        aux1 (ig,1:nbnd_occ(ikk)) = g2kin (ig) * evq (ig, 1:nbnd_occ(ikk))
     END DO
     !
     !
     IF (noncolin) THEN
        DO ig = 1, npwq
           aux1 (ig+npwx,1:nbnd_occ(ikk)) = g2kin (ig)* &
                                  evq (ig+npwx, 1:nbnd_occ(ikk))
        END DO
     END IF
     !
     DO ibnd= 1, nbnd_occ(ikk)
        eprectot (ibnd, nbase+ik) = 1.35d0 * zdotc(npwx*npol,evq(1,ibnd),1,aux1(1,ibnd),1)
     END DO
     !
  END DO
  !
  CALL mp_sum   ( eprectot, inter_pool_comm )
  !
  DEALLOCATE( aux1 )
  !
  CALL stop_clock( 'gwq_init' )
  !
  RETURN
  !
END SUBROUTINE gwq_init

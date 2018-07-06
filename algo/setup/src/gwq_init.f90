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
SUBROUTINE gwq_init(coulomb)
!----------------------------------------------------------------------------
  !
  !     0) initializes the structure factors
  !     2) Computes the non-linear core correction in cases where that is
  !     required.
  !     3) Computes the preconditioner for the linear system solver.
  !
  !
  USE becmod,               ONLY : calbec
  USE buffers,              ONLY : get_buffer
  USE constants,            ONLY : eps8, tpi
  USE control_lr,           ONLY : lgamma
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat, tau
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, nkstot, ngk
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE mp,                   ONLY : mp_sum
  USE mp_pools,             ONLY : my_pool_id, npool, kunit
  USE noncollin_module,     ONLY : npol
  USE qpoint,               ONLY : xq, npwq, nksq, eigqts, ikks, ikqs
  USE wvfct,                ONLY : npwx, npw, nbnd
  !
  IMPLICIT NONE
  !
  !> are we in the Coulomb or the self-energy step?
  LOGICAL, INTENT(IN) :: coulomb
  ! ... local variables

  INTEGER :: ik, ikq, ipol, ikk, na
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

  DO ik = 1, nksq
     !
     IF (coulomb) THEN
       ikk  = ikks(ik)
       npw  = ngk(ikk)
     END IF
     !
     ikq  = ikqs(ik)
     npwq = ngk(ikq)
     !
     IF (lsda .AND. coulomb) current_spin = isk(ikk)
     !
     ! ... if there is only one k-point evc, evq, npw, igk stay in memory
     !
     IF ( .NOT. lgamma .AND. coulomb ) THEN
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
  END DO
  !
  CALL stop_clock( 'gwq_init' )
  !
  RETURN
  !
END SUBROUTINE gwq_init

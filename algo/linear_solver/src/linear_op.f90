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
!> Contains the routine to apply the shifted Hamiltonian to a vector.
!!
!! @note this reimplements ch_psi_all from within LR_Modules
MODULE linear_op_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC linear_op

CONTAINS

!> Apply the shifted Hamiltonian.
!!
!! We evaluate the following equation
!! \f{equation}{
!!   A_n \psi_n = (H + \omega_n S + \alpha_{\text{pv}} P_{\text{v}}) \psi_n
!! \f}
!! Here H is the Hamiltonian, S is the overlap matrix, \f$P_{\text{v}}\f$ is the
!! projector onto the valence bands.
!!  
SUBROUTINE linear_op(current_k, num_g, omega, alpha_pv, psi, A_psi)

  USE constants,        ONLY: eps14
  USE control_lr,       ONLY: nbnd_occ
  USE kinds,            ONLY: dp
  USE noncollin_module, ONLY: noncolin, npol
  USE timing_module,    ONLY: time_linear_op
  USE wvfct,            ONLY: npwx, current_k_ => current_k

  IMPLICIT NONE

  !> Active k-point: Note this routine will set the wvfct module current_k to the same value.
  INTEGER,     INTENT(IN)  :: current_k

  !> Number of G vectors at this k-point; only this part of \f$\psi\f$ is used.
  INTEGER,     INTENT(IN)  :: num_g

  !> The frequency \f$\omega\f$ by which the Hamiltonian is shifted.
  COMPLEX(dp), INTENT(IN)  :: omega(:)

  !> The prefactor \f$\alpha_{\text{pv}}\f$ in front of the projector.
  REAL(dp),    INTENT(IN)  :: alpha_pv

  !> The vector \f$\psi\f$ to which the shifted Hamiltonian is applied.
  COMPLEX(dp), INTENT(IN)  :: psi(:,:)

  !> On output, the vector \f$(H + \omega S + \alpha_{\text{pv}} P_{\text{v}}) \psi\f$
  COMPLEX(dp), INTENT(OUT) :: A_psi(:,:)

  !> total vector size (npwx * npol).
  INTEGER vec_size

  !> Number of bands.
  INTEGER num_band

  !> Overlap matrix applied to psi.
  COMPLEX(dp), ALLOCATABLE :: spsi(:,:)

  !> loop variables
  INTEGER iband

  !> complex constants
  COMPLEX(dp), PARAMETER :: one = 1

  CALL start_clock(time_linear_op)

  ! determine some helper variables
  vec_size = npwx * npol
  num_band = SIZE(omega)

  !
  ! sanity check of input
  !
  ! psi must have dimension npwx*npol, num_band
  IF (SIZE(psi, 1) /= vec_size) &
    CALL errore(__FILE__, "First dimension of vector psi should be npwx * npol", 1)
  IF (SIZE(psi, 2) /= num_band) &
    CALL errore(__FILE__, "Second dimension of vector psi should be num_band", 1)
  ! A_psi and psi must have same shape
  IF (ANY(SHAPE(psi) /= SHAPE(A_psi))) &
    CALL errore(__FILE__, "Input and output vector must have same shape", 1)

  !
  ! apply Hamiltonian and overlap matrix to psi 
  !
  current_k_ = current_k

  ALLOCATE(spsi(vec_size, num_band))

  ! note: h_psi and s_psi will initialize the arrays
  CALL h_psi (npwx, num_g, num_band, psi, A_psi)
  CALL s_psi (npwx, num_g, num_band, psi, spsi)

  !
  ! we compute (H + omega S) psi
  !
  DO iband = 1, num_band
    CALL ZAXPY(num_g, omega(iband), spsi(1, iband), 1, A_psi(1, iband), 1)
    IF (noncolin) &
      CALL ZAXPY(num_g, omega(iband), spsi(npwx + 1, iband), 1, A_psi(npwx + 1, iband), 1)
  END DO ! iband

  !
  ! we compute the projector in the valence band
  !
  IF (ABS(alpha_pv) > eps14) THEN
    CALL projector_psi(num_g, nbnd_occ(current_k), alpha_pv, spsi)
    DO iband = 1, num_band
      CALL ZAXPY(num_g, one, spsi(1, iband), 1, A_psi(1, iband), 1)  
      IF (noncolin) &
        CALL ZAXPY(num_g, one, spsi(1 + npwx, iband), 1, A_psi(1 + npwx, iband), 1)  
    END DO
  END IF

  DEALLOCATE(spsi)

  CALL stop_clock (time_linear_op)

END SUBROUTINE linear_op

!> Evaluate the term \f$\alpha_{\text{pv}} P_{\text{v}} \psi\f$.
SUBROUTINE projector_psi(num_g, num_band_occ, alpha_pv, spsi)

  USE becmod, ONLY: becp, calbec
  USE kinds,  ONLY: dp
  USE eqv,    ONLY: evq
  USE uspp,   ONLY: vkb
  USE wvfct,  ONLY: npwx, nbnd 
 
  !> Number of G vectors at this k-point; only this part of \f$\psi\f$ is used.
  INTEGER,     INTENT(IN)    :: num_g

  !> Number of occupied bands
  INTEGER,     INTENT(IN)    :: num_band_occ

  !> prefactor before projector
  REAL(dp),    INTENT(IN)    :: alpha_pv

  !> on input:  overlap matrix applied to psi
  !! on output: the projector applied to psi
  COMPLEX(dp), INTENT(INOUT) :: spsi(:,:)

  !> total vector size (npwx * npol).
  INTEGER vec_size

  !> Number of bands.
  INTEGER num_band

  !> projectors
  COMPLEX(dp), ALLOCATABLE :: proj(:,:)

  !> projector applied to psi
  COMPLEX(dp), ALLOCATABLE :: P_psi(:,:)

  !> complex value of alpha
  COMPLEX(dp) c_alpha_pv

  !> complex constants
  COMPLEX(dp), PARAMETER :: one = 1, zero = 0

  ! set helper variables
  vec_size = SIZE(spsi, 1)
  num_band = SIZE(spsi, 2)
  c_alpha_pv = alpha_pv

  ALLOCATE(proj(nbnd, num_band))    
  ALLOCATE(P_psi(vec_size, num_band))

  ! proj = alpha evq* (S psi)
  CALL ZGEMM('C', 'N', num_band_occ, num_band, vec_size, c_alpha_pv, evq, &
             vec_size, spsi, vec_size, zero, proj, nbnd)

  ! P psi = evq proj
  CALL ZGEMM('N', 'N', vec_size, num_band, num_band_occ, one, evq, &
             vec_size, proj, nbnd, zero, P_psi, vec_size)

  ! apply S again
  CALL calbec(num_g, vkb, P_psi, becp, num_band)
  CALL s_psi(npwx, num_g, num_band, P_psi, spsi)

END SUBROUTINE projector_psi
 
END MODULE linear_op_module

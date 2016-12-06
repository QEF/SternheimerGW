!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
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
!> This subroutine is the main driver of the COULOMB self consistent cycle
!! which calculates the dielectric matrix by generating the density response
!! to a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G).
!! The dielectric matrix is given by:
!! eps_{q}^{-1}(G,G',iw) = (\delta_{GG'} + drhoscfs^{scf}_{G,G',iw})
!!
!! @note: this is solver for the special case q + G = 0
!-----------------------------------------------------------------------
SUBROUTINE coulomb_q0G0(config, eps_m) 
!-----------------------------------------------------------------------
  USE constants,        ONLY : eps8
  USE control_gw,       ONLY : solve_direct, niter_gw
  USE eqv_gw,           ONLY : drhoscfs, dvbare
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE freq_gw,          ONLY : fiu, nfs
  USE gvecs,            ONLY : nls
  USE gvect,            ONLY : g
  USE io_global,        ONLY : stdout, ionode
  USE noncollin_module, ONLY : nspin_mag
  USE kinds,            ONLY : dp
  USE qpoint,           ONLY : xq
  USE select_solver_module, ONLY : select_solver_type
  USE solve_module,     ONLY : solve_linter

  IMPLICIT NONE

  !> stores the configuration of the linear solver for the screened Coulomb interaction
  TYPE(select_solver_type), INTENT(IN) :: config

  !> the screened coulomb interaction
  COMPLEX(dp), INTENT(OUT) :: eps_m(nfs)

  !> complex constant of zero
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND = dp)

  !> complex constant of one
  COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND = dp)

  !> square of the length |q + G|
  REAL(dp) qg2

  !> index of G (fixed to the G = 0)
  INTEGER, PARAMETER :: ig = 1

  !> index of G'
  INTEGER igp

  !> index of the frequency
  INTEGER iw

  !> number of iterations - set to 1 for direct solver
  INTEGER num_iter

  !> format used to print eps/inveps for direct/iterative solver
  CHARACTER(LEN=:), ALLOCATABLE :: format_str

  ! we use the frequency as equivalent of the perturbation in phonon
  ALLOCATE (drhoscfs(dfftp%nnr, nspin_mag, nfs))
  IF (nspin_mag /= 1) CALL errore(__FILE__, "magnetic calculation not implemented", 1)

  ! determine number of iterations and set format string
  IF (solve_direct) THEN
    num_iter = 1
    format_str = '(8x,"eps_{GG}(q,w) = ", 2f12.5)'
  ELSE IF (niter_gw > 1) THEN
    num_iter = niter_gw
    format_str = '(5x,"inveps_{GG}(q,w) = ", 2f12.5)'
  ELSE
    CALL errore(__FILE__, "for iterative solver, we need to use more iterations", 1)
    format_str = '(*)'
  END IF
  !
  ! square of length |q + G|
  qg2 = SUM((g(:, ig) + xq)**2)
  !
  ! sanity check - qg2 must not vanish
  IF (qg2 < eps8) THEN
    CALL errore(__FILE__, "the solver cannot converge for q + G = 0", 1)
  ELSE
    !
    ! initialize the potential for a single G to 1
    drhoscfs = zero
    dvbare   = zero
    dvbare(nls(ig)) = one
    !
    ! potential in real space
    CALL invfft('Smooth', dvbare, dffts)
    !
    ! solve for linear response due to this perturbation
    CALL solve_linter(config, num_iter, dvbare, fiu(:nfs), drhoscfs)
    !
    ! back to reciprocal space
    CALL fwfft('Smooth', dvbare, dffts)
    !
    ! loop over frequencies
    DO iw = 1, nfs
      !
      ! evaluate response in reciprocal space
      CALL fwfft ('Dense', drhoscfs(:, 1, iw), dffts)
      !
      ! copy to output array
      eps_m(iw) = drhoscfs(nls(ig), 1, iw)
      !
      ! for the direct solver at bare potential in diagonal
      IF (solve_direct) THEN
        eps_m(iw) = eps_m(iw) + dvbare(nls(ig))
      END IF
      !
      ! print the diagonal element
      IF (ionode) THEN
        igp = nls(ig)
        WRITE(stdout, format_str) drhoscfs(igp, 1, iw) + dvbare(igp)
      END IF
      !
    END DO ! iw
    !
  END IF

  DEALLOCATE(drhoscfs)

END SUBROUTINE coulomb_q0G0

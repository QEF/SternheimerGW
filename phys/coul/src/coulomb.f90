!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
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
!> This subroutine is the main driver of the COULOMB self consistent cycle
!! which calculates the dielectric matrix by generating the density response
!! to a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G).
!! The dielectric matrix is given by:
!! eps_{q}^{-1}(G,G',iw) = (\delta_{GG'} + drhoscfs^{scf}_{G,G',iw})
!-----------------------------------------------------------------------
SUBROUTINE coulomb(config, igstart, num_g_corr, num_task, scrcoul) 
!-----------------------------------------------------------------------
  USE constants,        ONLY : eps8
  USE control_gw,       ONLY : solve_direct, niter_gw
  USE eqv,              ONLY : drhoscfs
  USE eqv_gw,           ONLY : dvbare
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE freq_gw,          ONLY : fiu, nfs
  USE gvect,            ONLY : g
  USE gwsymm,           ONLY : ig_unique
  USE io_global,        ONLY : stdout, ionode
  USE noncollin_module, ONLY : nspin_mag
  USE kinds,            ONLY : dp
  USE qpoint,           ONLY : xq
  USE select_solver_module, ONLY : select_solver_type
  USE solve_module,     ONLY : solve_linter
  USE timing_module,    ONLY : time_coulomb

  IMPLICIT NONE

  !> stores the configuration of the linear solver for the screened Coulomb interaction
  TYPE(select_solver_type), INTENT(IN) :: config

  !> first index of the G vector evaluated on this process
  INTEGER, INTENT(IN) :: igstart

  !> number of G vectors on all processes
  INTEGER, INTENT(IN) :: num_g_corr

  !> number of G vector evaluated by this process
  INTEGER, INTENT(IN) :: num_task

  !> the screened coulomb interaction
  COMPLEX(dp), INTENT(OUT) :: scrcoul(num_g_corr, nfs, num_task)

  !> actual index inside the array
  INTEGER indx

  !> complex constant of zero
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND = dp)

  !> complex constant of one
  COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND = dp)

  !> square of the length |q + G|
  REAL(dp) qg2

  !> index of G and G'
  INTEGER ig, igp

  !> index of the frequency
  INTEGER iw

  !> number of iterations - set to 1 for direct solver
  INTEGER num_iter

  !> format used to print eps/inveps for direct/iterative solver
  CHARACTER(LEN=:), ALLOCATABLE :: format_str

  !> determine the time since starting the routine
  REAL(dp) get_clock

  !> time at the start of the routine
  REAL(dp) start_time

  start_time = get_clock(time_coulomb)

  ! we use the frequency as equivalent of the perturbation in phonon
  ALLOCATE (drhoscfs(dfftp%nnr, nspin_mag, nfs))
  IF (nspin_mag /= 1) CALL errore(__FILE__, "magnetic calculation not implemented", 1)
  scrcoul = zero

  ! determine number of iterations and set format string
  IF (solve_direct) THEN
    num_iter = 1
    format_str = '(8x,"eps_{GG}(q,w) = ", 2f12.5, f9.2, a)'
  ELSE IF (niter_gw > 1) THEN
    num_iter = niter_gw
    format_str = '(5x,"inveps_{GG}(q,w) = ", 2f12.5, f9.2, a)'
  ELSE
    CALL errore(__FILE__, "for iterative solver, we need to use more iterations", 1)
    format_str = '(*)'
  END IF

  !
  ! loop over tasks = G vectors done in this image
  !
  DO indx = 1, num_task
    !
    ! determine index of G vector
    ig = igstart + indx - 1
    !
    ! square of length |q + G|
    qg2 = SUM((g(:, ig_unique(ig)) + xq)**2)
    !
    ! the case q + G = 0 is treated by separate routine
    IF (qg2 < eps8) CYCLE
    !
    ! initialize the potential for a single G to 1
    drhoscfs = zero
    dvbare   = zero
    dvbare(dffts%nl(ig_unique(ig))) = one
    !
    ! potential in real space
    CALL invfft('Rho', dvbare, dffts)
    !
    ! solve for linear response due to this perturbation
    CALL solve_linter(config, num_iter, dvbare, fiu(:nfs), drhoscfs)
    !
    ! back to reciprocal space
    CALL fwfft('Rho', dvbare, dffts)
    !
    ! loop over frequencies
    DO iw = 1, nfs
      !
      ! evaluate response in reciprocal space
      CALL fwfft ('Rho', drhoscfs(:, 1, iw), dffts)
      !
      ! copy to output array
      DO igp = 1, num_g_corr
        scrcoul(igp, iw, indx) = drhoscfs(dffts%nl(igp), 1, iw)
      END DO ! igp
      !
      ! for the direct solver at bare potential in diagonal
      IF (solve_direct) THEN
        igp = ig_unique(ig)
        scrcoul(igp, iw, indx) = scrcoul(igp, iw, indx) + dvbare(dffts%nl(igp))
      END IF
      !
      ! print the diagonal element + timing at the last frequency
      IF (ionode) THEN
        igp = dffts%nl(ig_unique(ig))
        IF (iw == nfs) THEN
          WRITE(stdout, format_str) drhoscfs(igp, 1, iw) + dvbare(igp), &
            get_clock(time_coulomb) - start_time, "s"
        ELSE
          WRITE(stdout, format_str) drhoscfs(igp, 1, iw) + dvbare(igp)
        END IF
      END IF
      !
    END DO ! iw
    !
  END DO ! indx

  DEALLOCATE(drhoscfs)

END SUBROUTINE coulomb

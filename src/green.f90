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
!> Wraps the routines used to evaluate the Green's function of the system.
MODULE green_module

  IMPLICIT NONE

  PUBLIC green_function
  PRIVATE

CONTAINS

  !> Evaluate the Green's function of the system.
  !!
  !! We solve the linear equation in reciprocal space to obtain the Green's
  !! function \f$G\f$
  !! \f{equation}{
  !!   (H_k(G') - \omega) G_k(G,G',\omega) = -\delta_{G,G'}~,
  !! \f}
  !! where \f$H_k\f$ is the Hamiltonian at a certain k-point, \f$\omega\f$ is
  !! the frequency, and \f$\delta_{G,G'}\f$ is the Kronecker delta.
  !!
  SUBROUTINE green_function(comm, multishift, lmax, threshold, igk, num_g, omega, green)

    USE bicgstab_module, ONLY: bicgstab
    USE kinds,           ONLY: dp
    USE mp,              ONLY: mp_sum
    USE parallel_module, ONLY: parallel_task
    USE timing_module,   ONLY: time_green

    !> Parallelize the calculation over this communicator
    INTEGER,     INTENT(IN)  :: comm

    !> Use the multishift solver to determine the Green's function
    LOGICAL,     INTENT(IN)  :: multishift

    !> Depth of the GMRES part of the BiCGstab(l) algorithm.
    INTEGER,     INTENT(IN)  :: lmax

    !> Threshold for the convergence of the linear system.
    REAL(dp),    INTENT(IN)  :: threshold

    !> The list from global G vector order to current k-point.
    !! @note this should be reduced to the correlation cutoff
    INTEGER,     INTENT(IN)  :: igk(:)

    !> The number of G-vectors at the k-point.
    INTEGER,     INTENT(IN)  ::  num_g

    !> The list of frequencies for which the Green's function is evaluated.
    COMPLEX(dp), INTENT(IN)  :: omega(:)

    !> The Green's function of the system.
    COMPLEX(dp), INTENT(OUT) :: green(:,:,:)

    !> distribution of the tasks over the process grid
    INTEGER,     ALLOCATABLE :: num_task(:)

    !> The right hand side of the linear equation
    COMPLEX(dp), ALLOCATABLE :: bb(:)

    !> Helper array storing the result of the linear equation
    COMPLEX(dp), ALLOCATABLE :: green_part(:,:)

    !> The number of frequencies.
    INTEGER num_freq

    !> The number of G-vectors for the correlation
    INTEGER num_g_corr

    !> First and last G-vector done on this process
    INTEGER ig_start, ig_stop

    !> loop variable for G loop
    INTEGER ig

    !> loop over frequencies
    INTEGER ifreq

    !> check error in array allocation
    INTEGER ierr

    !> complex zero
    COMPLEX(dp), PARAMETER :: zero = 0.0_dp

    !> complex one
    COMPLEX(dp), PARAMETER :: one = 1.0_dp

    CALL start_clock(time_green)

    ! determine helper variables
    num_g_corr = SIZE(green, 1)
    num_freq = SIZE(omega)

    ! sanity test of the input
    IF (SIZE(green, 2) /= num_g_corr) &
      CALL errore(__FILE__, "Green's function should be num_g_c x num_g_c x num_freq", 1)
    IF (SIZE(green, 3) /= num_freq) &
      CALL errore(__FILE__, "Green's function and omega are inconsistent", 1)
    IF (SIZE(igk) /= num_g_corr) &
      CALL errore(__FILE__, "mismatch in size of igk and Green's function", 1)

    ! allocate array for the right hand side
    ALLOCATE(bb(num_g), STAT = ierr)
    CALL errore(__FILE__, "could not allocate bb", ierr)

    ! initialize array
    green = zero

    ! allocate array for result of the linear system
    IF (multishift) THEN
      ALLOCATE(green_part(num_g, num_freq), STAT = ierr)
    ELSE
      ALLOCATE(green_part(num_g, 1), STAT = ierr)
    END IF
    CALL errore(__FILE__, "could not allocate green_part", ierr)

    ! parallelize over communicator
    CALL parallel_task(comm, num_g_corr, ig_start, ig_stop, num_task)
    DEALLOCATE(num_task)

    ! loop over all G-vectors
    DO ig = ig_start, ig_stop

      ! set right-hand side
      bb = zero
      bb(ig) = -one

      ! if multishift is set, we solve all frequencies at once
      IF (multishift) THEN

        ! solve the linear system
        CALL bicgstab(lmax, threshold, green_operator, bb, omega, green_part)

        ! copy from temporary array to result
        green(igk, igk(ig), :) = green_part(:num_g_corr, :)

      ! without multishift, we solve every frequency separately
      ELSE

        DO ifreq = 1, num_freq

          ! solve the linear system
          CALL bicgstab(lmax, threshold, green_operator, bb, omega(ifreq:ifreq), green_part)

          ! copy from temporary array to result
          green(igk, igk(ig), ifreq) = green_part(:num_g_corr, 1)

        END DO ! ifreq

      END IF

    END DO ! ig

    ! free memory
    DEALLOCATE(bb)
    DEALLOCATE(green_part)

    ! ccollect on single process
    CALL mp_sum(green, comm)

    CALL stop_clock(time_green)

  END SUBROUTINE green_function

  !> Wrapper for the linear operator call.
  !!
  !! Sets the necessary additional variables.
  SUBROUTINE green_operator(omega, psi, A_psi)

    USE kinds,            ONLY: dp
    USE klist,            ONLY: nks
    USE linear_op_module, ONLY: linear_op
    USE wvfct,            ONLY: npwx

    !> The initial shift
    complex(dp), INTENT(IN)  :: omega

    !> The input vector.
    COMPLEX(dp), INTENT(IN)  :: psi(:)

    !> The operator applied to the vector.
    COMPLEX(dp), INTENT(OUT) :: A_psi(:)

    !> Local copy of omega, because linear_op require arrays.
    COMPLEX(dp) omega_(1)

    !> Local copies of psi and A_psi because linear_op requires arrays.
    COMPLEX(dp), ALLOCATABLE :: psi_(:,:), A_psi_(:,:)

    !> real value of zero
    REAL(dp), PARAMETER :: zero = 0.0_dp

    !> Error code.
    INTEGER ierr

    !> number of G vectors
    INTEGER num_g

    ! determine number of G vectors
    num_g = SIZE(psi)

    !> allocate temporary array
    ALLOCATE(psi_(npwx, 1), A_psi_(npwx, 1), STAT = ierr)
    CALL errore(__FILE__, "could not allocate temporary arrays psi_ and A_psi_", ierr)

    !
    ! initialize helper
    !

    ! negative omega so that (H - omega) psi is calculated
    omega_ = -omega

    ! zero the elements outside of the definition
    psi_(:num_g, 1) = psi
    psi_(num_g + 1:, 1) = zero

    !
    ! apply the linear operator (H - w) psi
    !
    ! for the Green's function
    ! current_k = nks + 1
    ! alpha_pv = 0
    CALL linear_op(nks + 1, num_g, omega_, zero, psi_, A_psi_)

    ! extract result
    A_psi = A_psi_(:num_g,1)

    ! free memory
    DEALLOCATE(psi_, A_psi_)

  END SUBROUTINE green_operator

END MODULE green_module

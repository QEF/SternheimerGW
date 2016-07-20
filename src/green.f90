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
  SUBROUTINE green_function(multishift, igk, omega, green)

    USE kinds,         ONLY: dp
    USE timing_module, ONLY: time_green

    !> Use the multishift solver to determine the Green's function
    LOGICAL,     INTENT(IN)  :: multishift

    !> The list from global G vector order to current k-point.
    INTEGER,     INTENT(IN)  :: igk(:)

    !> The list of frequencies for which the Green's function is evaluated.
    COMPLEX(dp), INTENT(IN)  :: omega(:)

    !> The Green's function of the system.
    COMPLEX(dp), INTENT(OUT) :: green(:,:,:)

    !> The number of frequencies.
    INTEGER num_freq

    !> The number of G-vectors.
    INTEGER num_g

    CALL start_clock(time_green)

    ! determine helper variables
    num_g = SIZE(green, 1)
    num_freq = SIZE(omega)

    ! sanity test of the input
    IF (SIZE(green, 2) /= num_g) &
      CALL errore(__FILE__, "Green's function should be num_g x num_g x num_freq", 1)
    IF (SIZE(green, 3) /= num_freq) &
      CALL errore(__FILE__, "Green's function and omega are inconsistent", 1)
    IF (SIZE(igk) /= num_g) &
      CALL errore(__FILE__, "mismatch in size of igk and Green's function", 1)

    CALL stop_clock(time_green)

  END SUBROUTINE green_function

END MODULE green_module

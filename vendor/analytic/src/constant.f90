! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

MODULE constant_module

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

  REAL(dp), PARAMETER :: eps6 = 1e-6_dp
  REAL(dp), PARAMETER :: eps8 = 1e-8_dp
  REAL(dp), PARAMETER :: eps10 = 1e-10_dp
  REAL(dp), PARAMETER :: eps12 = 1e-12_dp
  REAL(dp), PARAMETER :: eps14 = 1e-14_dp

  REAL(dp), PARAMETER    :: zero = 0.0_dp
  COMPLEX(dp), PARAMETER :: c_zero = CMPLX(zero, zero, KIND=dp)

  REAL(dp), PARAMETER    :: one = 1.0_dp
  COMPLEX(dp), PARAMETER :: c_one = CMPLX(one, zero, KIND=dp)

  COMPLEX(dp), PARAMETER :: imag = CMPLX(zero, one, KIND=dp)

  REAL(dp), PARAMETER :: two = 2.0_dp

  REAL(dp), PARAMETER :: quarter_pi = ATAN(1.0_dp)
  REAL(dp), PARAMETER :: half_pi = two * quarter_pi
  REAL(dp), PARAMETER :: pi = two * half_pi
  REAL(dp), PARAMETER :: two_pi = two * pi
  REAL(dp), PARAMETER :: four_pi = two * two_pi

  REAL(dp), PARAMETER :: sqrt_two = SQRT(two)
  REAL(dp), PARAMETER :: sqrt_half = one / sqrt_two

  COMPLEX(dp), PARAMETER :: two_pi_i = two_pi * imag

END MODULE constant_module

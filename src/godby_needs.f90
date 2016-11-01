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
!> Implements the Godby-Needs Plasmon pole model.
!!
!! The plasmon pole is defined as
!! \f{equation}{
!!   W_{GG'}(\omega) = \frac{A_{GG'}}{\omega + \tilde \omega_{GG'}} 
!!                   - \frac{A_{GG'}}{\omega - \tilde \omega_{GG'}}
!! \f}
!! Godby and Needs suggested to fit this to the values evaluated at two
!! frequencies \f$\omega = 0\f$ and \f$\omega = i \omega_{\text{p}}\f$. Then
!! we can determine the constants
!! \f{equation}{
!!   A_{GG'} = \frac12 W_{GG'}(0) \tilde \omega
!! \f}
!! and
!! \f{equation}{
!!   \tilde\omega_{GG'} = \sqrt{
!!     \frac{W_{GG'}(\omega_{\text{p}})}{W_{GG'}(0) - W_{GG'}(\omega_{\text{p}})]}
!!   } \omega_{\text{p}}~.
!! \f}
MODULE godby_needs_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE godby_needs_coeffs(N, z, u, a)

  USE kinds,       ONLY: dp
  USE constants,   ONLY: eps8

  IMPLICIT NONE

  INTEGER     :: N
  COMPLEX(dp) :: z(N), u(N)
  COMPLEX(dp) :: a(N)

  !> complex constant of 0
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

  !> large complex number
  COMPLEX(dp), PARAMETER :: large = CMPLX(10.0_dp, 0.0_dp, KIND=dp)

  !! Instead of a higher order pad\'e approximant to the self energy we just calculate two 
  !! frequency points one at w = 0 and one at w = i\Omega_{p}. and fit it to the plasmon pole
  !! model: this is the godby_needs plasmon pole.
  !! The assumed form of the function is:
  !! \epsilon^{-1}_{\G,\G'}(\q,\omega) = \frac{A_{GG'(q)}}{\omega - \omegatilde + idelta} -
  !!                                     \frac{A_{GG'(q)}}{\omega + \omegatilde - idelta}
  a(:) = zero 
  !! Currently using the same criterion as in SaX
  !! this essentially checks if the real part of the pole
  !! is smaller than the imaginary part of the pole and if so
  !! kills the pole...
  !! for inclusion of plasmon pole parameters...
  !! might be more appropriate to think about the conditions some more.
  !! this essentially says:
  !! \sqrt(a + ib) = c + id
  !! a + ib = c**2 - d**2 +2icd
  !! if a < 0 :: d > c
  !! so they kill that pole...
  !! essentially we cast aside any heavily damped oscillations
  !! (which would not effect the real part of the selfenergy anyway...
  !! a(1) = \tilde(\omega) a(2) = R

  IF (ABS(REAL(u(1) - u(2))) < eps8) THEN 
    ! We zero the weight of the pole and place the pole way out
    ! on the real axis to avoid numerical instability.
    ! although this isn't really that far out....
    a(1) = large
    a(2) = zero

  ELSE IF (REAL(u(2) / (u(1) - u(2))) < 0.0_dp) THEN
    ! case for wings having been zerod
    a(1) = large
    a(2) = zero 

  ELSE
    ! \tilde{\omega}:
    a(1) = z(2) * SQRT(REAL(u(2) / (u(1) - u(2))))
    !(A_{GG'qq}):
    a(2) = -u(1) * a(1)**2

  END IF

END SUBROUTINE godby_needs_coeffs

  !> construct the screened Coulomb potential
  FUNCTION godby_needs_model(freq, coeff) RESULT (res)

    USE kinds,     ONLY: dp
    USE constants, ONLY: eps8

    !> the frequency at which we want to determine the Coulomb potential
    COMPLEX(dp), INTENT(IN) :: freq

    !> the strength and position of the pole
    COMPLEX(dp), INTENT(IN) :: coeff(2)

    !> the Coulomb potential at the given frequency
    COMPLEX(dp) res

    res = coeff(2) / (freq**2 - coeff(1)**2)

  END FUNCTION godby_needs_model

END MODULE godby_needs_module

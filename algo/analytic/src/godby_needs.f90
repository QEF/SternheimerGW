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
!!   A_{GG'} = \frac12 W_{GG'}(0) \tilde \omega_{GG'}
!! \f}
!! and
!! \f{equation}{
!!   \tilde\omega_{GG'} = \sqrt{
!!     \frac{W_{GG'}(\omega_{\text{p}})}{W_{GG'}(0) - W_{GG'}(\omega_{\text{p}})}
!!   } \omega_{\text{p}}~.
!! \f}
MODULE godby_needs_module

  IMPLICIT NONE

CONTAINS

  !> evaluate the Godby-Needs coefficients
  SUBROUTINE godby_needs_coeffs(omega_p, coulomb)

    USE kinds,     ONLY: dp
    USE constants, ONLY: eps8

    !> the Plasmon frequency at which the Coulomb potential is evaluated
    REAL(dp),    INTENT(IN)    :: omega_p

    !> *on input* the screened Coulomb potential at \f$\omega = 0\f$ and
    !! \f$\omega = i \omega_{\text{p}}\f$ <br>
    !! *on output* the Godby-Needs coefficients used to evaluate the
    !! screened Coulomb potential at different frequencies
    COMPLEX(dp), INTENT(INOUT) :: coulomb(:,:,:)

    !> constant of 1/2
    REAL(dp),    PARAMETER :: half = 0.5_dp

    !> complex constant of 0
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    !> counter on the G vectors
    INTEGER ig, igp

    !> set the current element to 0
    LOGICAL set_zero

    !> complex work variable
    COMPLEX(dp) work

    !
    ! sanity check of the input
    !
    ! plasmon frequency positive
    IF (omega_p < 0) CALL errore(__FILE__, "plasmon frequency must be positive", 1)
    !
    IF (SIZE(coulomb, 3) /= 2) &
      CALL errore(__FILE__, "must provide exactly 2 frequencies", 1)

    DO igp = 1, SIZE(coulomb, 2)
      DO ig = 1, SIZE(coulomb, 1)
        !
        ! work = W(0) - W(wp)
        work = coulomb(ig, igp, 1) - coulomb(ig, igp, 2)
        set_zero = ABS(REAL(work)) < eps8
        !
        IF (.NOT.set_zero) THEN
          !
          ! work = W(wp) / (W(0) - W(wp))
          work = coulomb(ig, igp, 2) / work
          set_zero = REAL(work) < eps8
          !
        END IF
        !
        IF (.NOT.set_zero) THEN
          !         _____________
          !        /    W(wp)
          ! w~ =  / ------------   wp
          !      V  W(0) - W(wp)
          !
          coulomb(ig, igp, 2) = SQRT(work) * omega_p
          !     1
          ! A = - W(0) w~
          !     2
          coulomb(ig, igp, 1) = half * coulomb(ig, igp, 1) * coulomb(ig, igp, 2)
          !
        ELSE
          !
          coulomb(ig, igp, 1) = zero
          coulomb(ig, igp, 2) = zero
          !
        END IF
        !
      END DO ! ig
    END DO ! igp

  END SUBROUTINE godby_needs_coeffs

  !> construct the screened Coulomb potential
  FUNCTION godby_needs_model(freq, coeff) RESULT (res)

    USE kinds,      ONLY: dp
    USE constants,  ONLY: eps8

    !> the frequency at which we want to determine the Coulomb potential
    COMPLEX(dp), INTENT(IN) :: freq

    !> the strength and position of the pole
    COMPLEX(dp), INTENT(IN) :: coeff(2)

    !> the Coulomb potential at the given frequency
    COMPLEX(dp) res

    !> complex constant of 1
    COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

    !> complex constant of 0
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    !          A        A
    ! W(w) = ------ + ------
    !        w~ + w   w~ - w
    IF (ABS(coeff(1)) > eps8) THEN
      !
      res = coeff(1) * (one / (coeff(2) + freq) + one / (coeff(2) - freq))
      !
    ELSE
      res = zero
    END IF

  END FUNCTION godby_needs_model

END MODULE godby_needs_module

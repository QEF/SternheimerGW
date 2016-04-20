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
SUBROUTINE godby_needs_coeffs (N, z, u, a)
USE kinds,                     ONLY : DP
USE mp_world,   ONLY : mpime
USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps4, pi, eps8

IMPLICIT NONE

integer  :: p, N, i 
complex(DP) :: z(N), u(N)
complex(DP) :: a(N)
real(DP) :: ar, ai

!Instead of a higher order pad\'e approximant to the self energy we just calculate two 
!frequency points one at w = 0 and one at w = i\Omega_{p}. and fit it to the plasmon pole
!model: this is the godby_needs plasmon pole.
!The assumed form of the function is:
!\epsilon^{-1}_{\G,\G'}(\q,\omega) = \frac{A_{GG'(q)}}{\omega - \omegatilde + idelta} -
!                                    \frac{A_{GG'(q)}}{\omega + \omegatilde - idelta}
     a(:) = DCMPLX(0.0d0, 0.0d0)
!Currently using the same criterion as in SaX
!this essentially checks if the real part of the pole
!is smaller than the imaginary part of the pole and if so
!kills the pole...
!for inclusion of plasmon pole parameters...
!might be more appropriate to think about the conditions some more.
!this essentially says:
!\sqrt(a + ib) = c + id
!a + ib = c**2 - d**2 +2icd
!if a < 0 :: d > c
!so they kill that pole...
!essentially we cast aside any heavily damped oscillations
!(which would not effect the real part of the selfenergy anyway...
!a(1) = \tilde(\omega) a(2) = R

   if(real(u(1)-u(2)).eq.(0.0d0)) then 
       !We zero the weight of the pole and place the pole way out
       !on the real axis to avoid numerical instability.
       !although this isn't really that far out....
           a(1) = 10.0
          !a(2) = -u(1)*a(1)**2 
           a(2) = 0.0d0
          !write(1000+mpime,'("WINGS ZEROD")') 
   else if(real(u(2)/(u(1)-u(2))).lt.(0.0d0)) then
    !case for wings having been zerod
           a(1) = 10.0
           a(2) = 0.0
          !a(2) = -u(1)*a(1)**2 
          !write(1000+mpime,'("less than 0")') 
   else
!\tilde{\omega}:
!     a(1) = z(2)*SQRT(real(u(2)/(u(1)-u(2))))
!(A_{GG'qq}):
!     a(2) = -((u(1)*a(1))/DCMPLX(2.0d0,0.0d0))
!new math
     a(1) = z(2)*SQRT(real(u(2)/(u(1)-u(2))))
     a(2) = -u(1)*a(1)**2
   endif
!Condition for catching nan's!
     do p = 1, N
        ar = real(a(p))
        ai = aimag(a(p))
        if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
           write(1000+mpime,'(2f12.7)') (z(i),i=1,N)
           write(1000+mpime,'(2f12.7)') (u(i),i=1,N)
           write(1000+mpime,'(2f12.7)') (a(i),i=1,N)
           !a(1) = DCMPLX(20.0d0,0.0d0)
           !a(2) = DCMPLX(0.0d0,0.0d0)
           a(1) = 10.0
           a(2) = -u(1)*a(1)**2 
        endif
     enddo
END SUBROUTINE godby_needs_coeffs

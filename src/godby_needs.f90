SUBROUTINE godby_needs_coeffs (N, z, u, a)
USE kinds,                     ONLY : DP
USE mp_global,   ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime
USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps4, pi, eps8

IMPLICIT NONE

complex(DP) :: z(N), u(N)
complex(DP) :: a(N)
real(DP) :: ar, ai
integer  :: p, N,i 

!Instead of a higher order pad\'e approximant to the self energy we just calculate two 
!frequency points one at w = 0 and one at w = i\Omega_{p}. and fit it to the plasmon pole
!model: this is the godby_needs plasmon pole.
!The sign here is wrong on the imaginary part...
!The assumed form of the function is:
!\epsilon^{-1}_{\G,\G'}(\q,\omega) = \frac{A_{GG'(q)}}{\omega - \omegatilde + idelta} -
!                                    \frac{A_{GG'(q)}}{\omega + \omegatilde - idelta}
!Except we solve only for the correlation energy:
!(\epsilon^{-1} - 1)
     a(:) = DCMPLX(0.0d0, 0.0d0)
!pole position(omegatilde):
!real frequency
    if(N.ne.2) then
      write(6,'("Only two freqs. allowed with godby needs")')
      STOP
    endif

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
!The exact meaning of all this eludes me...
!essentially we cast aside any heavily damped oscillations
!(which would not effect the real part of the selfenergy anyway...
!a(1) = \tilde(\omega) a(2) = R
   if(real(u(2)/(u(1)-u(2))).lt.0) then
       !We zero the weight of the pole and place the pole way out
       !on the real axis to avoid numerical instability.
       !although this isn't really that far out....
           a(1) = 20.0
           a(2) = -((u(1)*a(1))/DCMPLX(2.0d0,0.0d0))
!   else if (real(abs(u(1)-u(2))).lt.eps8) then
!           write(1000+mpime,'("zeroing eps")')
!          a(1) = 10.0d0
!          a(2) = 0.0
   else
!\tilde{\omega}:
     a(1) = z(2)*SQRT(real(u(2))/(u(1)-u(2)))
!(A_{GG'qq}):
     a(2) = -((u(1)*a(1))/DCMPLX(2.0d0,0.0d0))
   endif

!Condition for catching nan's!
     do p = 1, N
        ar = real(a(p))
        ai = aimag(a(p))
        if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
           write(1000+mpime,*) (z(i),i=1,N)
           write(1000+mpime,*) (u(i),i=1,N)
           write(1000+mpime,*) (a(i),i=1,N)
           a(1) = DCMPLX(10.0d0,0.0d0)
           a(2) = DCMPLX(0.0d0,0.0d0)
        endif
     enddo
END SUBROUTINE godby_needs_coeffs

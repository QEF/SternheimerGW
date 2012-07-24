SUBROUTINE godby_needs_coeffs (N, diag, z, u, a)
USE kinds,                     ONLY : DP
USE mp_global,   ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime
implicit none
complex(DP) :: z(N), u(N)
complex(DP) :: a(N)
real(DP) :: ar, ai
integer  :: p, N,i 
logical  :: diag

!Instead of a higher order pad\'e approximant to the self energy we just calculate two 
!frequency points one at w = 0 and one at w = i\Omega_{p}. and fit it to the plasmon pole
!model: this is the godby_needs plasmon pole.
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

    if(diag) then
      a(1) = z(2)*SQRT((u(2))/((u(1))-(u(2))))
    else
      a(1) = z(2)*SQRT((u(2))/(u(1)-u(2)))
    endif

!weight(A_{GG'qq})
     a(2) = -((u(1)*a(1))/DCMPLX(2.0d0,0.0d0))
     do p = 1, N
        ar = real(a(p))
        ai = aimag(a(p))
        if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
           write(1000+mpime,*) (z(i),i=1,N)
           write(1000+mpime,*) (u(i),i=1,N)
           write(1000+mpime,*) (a(i),i=1,N)
           a(2) = DCMPLX(0.0d0,0.0d0)
        endif
     enddo
     !write(1000+mpime,*)z, u
END SUBROUTINE godby_needs_coeffs

!SUBROUTINE godby_needs_eval(N, a, w, padapp)
!USE kinds,    ONLY : DP
!implicit none
 !----------------------------------------------------
 ! a(1:N) - coefficients for the PPM
 ! w      - point at which we need the approximant
 ! output
 ! 
 ! padapp - value of the approximant at the point w
 !----------------------------------------------------
! complex(DP) :: a(N)
! cone   = (1.0d0, 0.0d0)
! ci     = (0.0d0, 0.0d0)
! padapp = (a(2)/( w - a(1) + ci*eta )) - ((a(2))/(w + a(1) -ci*eta))
!END SUBROUTINE godby_needs_eval

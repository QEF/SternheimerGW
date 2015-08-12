  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
  ! General note:
  ! Lebegue, Arnaud, Alouani, and Blochel [PRB 67, 155208 (2003)]
  ! state that when they use Pade of order N = 12 (resulting in
  ! numerator or order (N-2)/2 = 5 and denominator N/2 = 6),
  ! they obtain extremely stable fits and the quasiparticle energies
  ! are essentially identical to those obtained using the contour
  ! deformation method.
  !
  ! using this sub:
  !
  ! integer :: N
  ! complex(DP) :: z(N), u(N), a(N), w, padapp
  !
  ! call pade_coeff ( N, z, u, a)
  ! call pade_eval ( N, z, a, w, padapp)
  !
  !-----------------------------------------------------------
  subroutine pade_coeff ( N, z, u, a)
  !-----------------------------------------------------------
  ! N-point Pade' approximant - find the Pade' coefficients
  !
  ! This subroutine uses the recursive algorithm described in
  ! HJ Vidberg and JW Serene, "Solving the Eliashberg equations
  ! by means of N-point Pade' approximants", J Low Temp Phys
  ! 29, 179 (1977). The notation adopted here is the same as
  ! in the above manuscript.
  !
  ! input
  !
  ! N      - order of the Pade' approximant
  ! z(1:N) - points at which the original function is known
  ! u(1:N) - values of the function at the z points
  !
  ! output
  !
  ! a(1:N) - coefficients of the continued fraction
  !-----------------------------------------------------------
  USE kinds,         ONLY : DP
  implicit none
  integer :: N
  complex(DP) :: z(N), u(N)
  complex(DP) :: g(N,N), a(N)
  real(DP) :: ar, ai
  complex(DP) :: tmp1, tmp2

  !complex(selevted_real_kind(18,50))?

  integer :: i, j, p
  !
  do p = 1, N
    if (p.eq.1) then
      do i = 1, N
         g (p,i) = u(i)
      enddo
    else
      do i = p, N
         tmp1 = g(p-1,p-1)/g(p-1,i)
         tmp2 = g(p-1,i)/g(p-1,i)
         g (p,i) = ( tmp1 - tmp2 ) / ( z(i) - z(p-1) )
      enddo
    endif
    a(p) = g (p,p)
  ! check whether a(p) is not NaN
    ar = real(a(p))
    ai = aimag(a(p))
    if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
          a(:) = (0.0d0, 0.0d0)
    endif
    !
  enddo
  !
  end subroutine pade_coeff

  !-----------------------------------------------------------
  subroutine pade_eval ( N, z, a, w, padapp)
  !-----------------------------------------------------------
  ! N-point Pade' approximant - evaluate the Pade' approximant
  !
  ! This subroutine uses the recursive algorithm described in
  ! HJ Vidberg and JW Serene, "Solving the Eliashberg equations
  ! by means of N-point Pade' approximants", J Low Temp Phys
  ! 29, 179 (1977). The notation adopted here is the same as
  ! in the above manuscript.
  !
  ! input
  !
  ! N      - order of the Pade' approximant
  ! z(1:N) - points at which the original function is known
  ! a(1:N) - coefficients of the continued fraction
  ! w      - point at which we need the approximant
  !
  ! output
  !
  ! padapp - value of the approximant at the point w
  !-----------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  implicit none
  integer :: N
  complex(DP) :: a(N), z(N), acap(0:N), bcap(0:N)
  complex(DP) :: w, padapp
  integer :: i
  real(DP) :: ar, ai

  !
  acap(0) = 0.d0
  acap(1) = a(1)
  bcap(0) = 1.d0
  bcap(1) = 1.d0
  !
  do i = 2, N
    acap(i) = acap(i-1) + (w-z(i-1)) * a(i) * acap(i-2)
    bcap(i) = bcap(i-1) + (w-z(i-1)) * a(i) * bcap(i-2)
  enddo
  padapp = acap(N)/bcap(N)

!Turning on pade catch.
  ar = real(padapp)
  ai = aimag(padapp)
  if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
    padapp = (0.0d0,0.0d0)
  endif
  !
  end subroutine pade_eval
  !-----------------------------------------------------------
  !

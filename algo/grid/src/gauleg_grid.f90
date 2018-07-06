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
SUBROUTINE gauleg_grid(x1, x2, x, w, n)

  USE constants,  ONLY : pi
  USE kinds,      ONLY : DP
     
  IMPLICIT NONE 
   
INTEGER, intent(in)  :: n            !number of points
INTEGER     :: i,j,m
REAL(DP)    :: x1, x2       !ranges
REAL(DP)    :: x(n), w(n)   !abcissa and weights
REAL(DP)    :: eps
Parameter (eps=3D-14) ! eps is the relative precision
REAL(DP)    :: p1, p2, p3, pp, xl, xm, z, z1

   m=(n+1)/2
   xm=0.5D0*(x2+x1)
   xl=0.5D0*(x2-x1)

! loop over the desired roots
   do i=1,m
     z=cos(pi*(i-.25D0)/(n+.5D0))
!   starting with the above approximation to the ith root, we enter the
!   main loop of refinement by Newton's method.
   1       continue
       p1=1D0
       p2=0D0
!     Loop up the recurrence relation to get the Legendre polynomial
!     evaluated at z.
       do j=1,n
         p3=p2
         p2=p1
         p1=((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
       enddo
!     p1 is now the desired Legendre polynomial. We next compute pp, its
!     derivative, by a standard rlation involving also p2, the polynomial
!     of one lower order
       pp=n*(z*p1-p2)/(z*z-1D0)
       z1=z
!     Newton's method
       z=z1-p1/pp
     if(abs(z-z1).gt.eps) goto 1
!   Scale the root to the desired interval, and put in its symmetric counter
!   part.
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
!   compute the weight and its symmetric counterpart
     w(i)=2D0*xl/((1D0-z*z)*pp*pp)
     w(n+1-i)=w(i)
   enddo
   return
END SUBROUTINE gauleg_grid

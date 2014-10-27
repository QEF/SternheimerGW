SUBROUTINE freq_gl()
  USE kinds,      ONLY : DP
  USE constants,  ONLY : RYTOEV, pi
  USE control_gw, ONLY : eta
     
  IMPLICIT NONE 
   
INTEGER     :: n            !number of points
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
END SUBROUTINE freq_gl

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
  !
  implicit none
  integer :: N
  complex(16) :: z(N), u(N)
  complex(16) :: g(N,N), a(N)
  ! g(p,i) = g_p (z_i) in the notation of Vidberg and Serene
  integer :: i, j, p
  real(16) :: ar, ai
  complex(16) :: tmp1, tmp2
  !
  do p = 1, N
    if (p.eq.1) then
      do i = 1, N
         g (p,i) = u(i)
      enddo
    else
      do i = p, N
   !     g (p,i) = ( g(p-1,p-1) - g(p-1,i) ) / &
   !               ( ( z(i) - z(p-1) ) * g (p-1,i) )
         !
         ! this seems necessary to avoid nasty NaN when
         ! still don't quite understand why the procedure
         ! becomes unstable - certainly it happens only
         ! when u(:) is very small
         !
         tmp1 = g(p-1,p-1)/g(p-1,i)
         tmp2 = g(p-1,i)/g(p-1,i)
         g (p,i) = ( tmp1 - tmp2 ) / ( z(i) - z(p-1) )
         !
      enddo
    endif
    a(p) = g (p,p)
    !
    ! check whether a(p) is not NaN
    !
    ar = real(a(p))
    ai = aimag(a(p))
!    if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
!       write(6,*) (z(i),i=1,N)
!       write(6,*) (u(i),i=1,N)
!       write(6,*) (a(i),i=1,N)
!       call error ('pade_coeff','one or more coefficients are NaN',1)
!    endif
    !
  enddo
  !
 end subroutine pade_coeff
  !
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
  implicit none
  integer :: N
  complex(16) :: a(N), z(N), acap(0:N), bcap(0:N)
  complex(16) :: w, padapp
  integer :: i
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
  !
  padapp = acap(N)/bcap(N)
  !
  end subroutine pade_eval
  !
  !-----------------------------------------------------------
  subroutine pade_poly ( N, z, a,m, num,den)
  !-----------------------------------------------------------
  ! N-point Pade' approximant -
  ! gives the coefficiens of the numerator and denominator polynomial
  !
  ! input
  !
  ! N      - order of the Pade' approximant
  ! z(1:N) - points at which the original function is known
  ! a(1:N) - coefficients of the continued fraction
  ! m      - dimension of coefficient vectors
  !
  ! output
  !
  ! num(1:m) - coefficients of the numerator polynomial 
  ! den(1:m) - coefficients of the denominator polynomial
  !-----------------------------------------------------------
  !
    implicit none
  integer :: N , m
  complex(16) :: a(N), z(N)
  complex(16) :: num(m), den(m)
  integer :: i,ii,j
  !
  complex(16) :: mat_num(N,N),mat_den(N,N),b(N)
  !
  ! prepare the coefficents of the continued fraction such that
  ! (z-z_i)*a'(i)=z*a'(i)-z(i)*a'(i) -> z*a(i)+b(i)
  do i=1,N-1
     b(i) = -z(i)*a(i+1)
  enddo
  !
  ! aux matrix in which the coeffiecients are recursively
  ! accumulated, such that c_i=mat(i,m+2)
  ! Could be done more efficiently.
  mat_num(:,:) = 0._16
  mat_den(:,:) = 0._16
!  ! the constant coefficient (c_0) is easiest treated as a a special case
  mat_num(1,1) = 0._16
  mat_num(1,2) = a(1)
!  !
  mat_den(1,1) = 1._16
  mat_den(1,2) = 1._16
  do j=3,N-1
     mat_num(1,j) = b(j-2)*mat_num(1,j-2)+mat_num(1,j-1)
     mat_den(1,j) = b(j-2)*mat_den(1,j-2)+mat_den(1,j-1)
  enddo
  num(1) = mat_num(1,N-1)
  den(1) = mat_den(1,N-1)
  !
  ! Now loop over all remaining coefficients
  do i=2,m
     do j=3,N-1
        mat_num(i,j) = b(j-2)*mat_num(i,j-2)+a(j-2+1)*mat_num(i-1,j-2)+mat_num(i,j-1)
        mat_den(i,j) = b(j-2)*mat_den(i,j-2)+a(j-2+1)*mat_den(i-1,j-2)+mat_den(i,j-1)
     enddo
     num(i) = mat_num(i,N-1)
     den(i) = mat_den(i,N-1)
  enddo
  !
end subroutine pade_poly
  !-----------------------------------------------------------

  !
  !-----------------------------------------------------------
  subroutine pade_poly_old ( N, z, a,m, poly)
  !-----------------------------------------------------------
  ! N-point Pade' approximant -
  ! gives the coefficiens of the denominator polynomial
  !
  ! This subroutine uses the recursive algorithm described in
  ! HJ Vidberg and JW Serene, "Solving the Eliashberg equation
  ! by means of N-point Pade' approximants", J Low Temp Phys
  ! 29, 179 (1977). The notation adopted here is the same as
  ! in the above manuscript.
  !
  ! input
  !
  ! N      - order of the Pade' approximant
  ! z(1:N) - points at which the original function is known
  ! a(1:N) - coefficients of the continued fraction
  ! m      - dimension of coefficient vector
  !
  ! output
  !
  !poly(1:m) - coefficients of the denominator polynomial
  !-----------------------------------------------------------
  !
    implicit none
  integer :: N , m
  complex(16) :: a(N), z(N)
  complex(8) :: poly(m), acap(0:N), bcap(0:N)
  complex(8) :: w
  integer :: i,ii,j
  !
  complex(8) :: mat(m,m),b(m)
  !

  ! solve linear system where coefficients are the unknowns
  mat(:,:) = 0.
  do i=1,m
     w = cmplx(real(i),0._8)
     bcap(0) = 1.d0
     bcap(1) = 1.d0
     do ii = 2, N
        bcap(ii) = bcap(ii-1) + (w-z(ii-1)) * a(ii) * bcap(ii-2)
     enddo
     !
     b(i) = bcap(N)
     !
     do j=1,m
        mat(i,j) = w**(j-1)
     enddo
  enddo
  !
  call invertc(m,mat)
  !
  ! matrix product explicitly, matmul didnt work (?)
  poly(:) = 0.
  do i=1,m
     do j=1,m
        poly(i) = poly(i) + mat(i,j)*b(j)
     enddo
  enddo
  !
  !
end subroutine pade_poly_old
  !-----------------------------------------------------------

subroutine test_pade(nw,pade_w)
  !
  !
  implicit none
  !
  ! Arguments
  integer       :: nw    ! number of imag frequencies
  complex(16)   :: pade_w(nw) ! values of imag frequencies
  !
  ! internal
  complex(16)   :: a_pade(nw)
  complex(16)   :: padapp, u(nw),w, p1, p2, p3, p4, w1, w2, w3, w4
  integer       :: i1,i2,iw
  complex(8), external :: rational_function
  !
  ! evaluate function at imaginary frequencies
  do iw=1,nw
     w = pade_w(iw)
print *, w
     u(iw) = rational_function(w)
  enddo
  !
  call pade_coeff( nw, pade_w, u, a_pade)
  !
  ! evaluate functiona real freuquencies
  do iw=1,3000
     w = (iw-1)*0.01
     call pade_eval( nw, pade_w, a_pade, w, padapp)
     !
     write(123,'(e12.6,2x,e12.6,2x,e12.6)') real(w) , real(padapp) , aimag(padapp)
     write(124,'(e12.6,2x,e12.6,2x,e12.6)') real(w) , real(rational_function(w)), aimag(rational_function(w))
     !
  enddo
stop
end subroutine test_pade


complex(8) function rational_function(w)
  complex(16) :: w, val
  real(16)    :: pol1,pol2,pol3,pol4
  !
  pol1 = 10+1.2354
  pol2 = 10+2.4764
  pol3 = 10+4.05421
  pol4 = 10+7.05421
  !
  w = w +cmplx(0.,0.1)
  rational_function = 1/(w**2-pol1**2)+1/(w**2-pol2**2) +1/(w**2-pol3**2)+1/(w**2-pol4**2)
  !
end function rational_function


program pade_tester
  !
  implicit none
  integer :: nw, i1, iw, nw1,m, i2, nskip, nw0, nw_save, task, count0
  complex(16), allocatable   :: pade_w(:), u(:), a_pade(:),num(:),den(:)
!  complex(8), allocatable :: num(:),den(:)
  logical :: lexist, sym_4, find_poles
  character(1028) :: name, name1,name2,name4
  complex(16)   :: padapp,w, Im
  real(8) :: dummy1, dummy2
  complex(8) :: dummyc1,dummyc2, dummyc3,dummyc4, win,pin
  complex(8), allocatable ::  poly(:), poles(:)
  !
  Im =cmplx(0._8,1._8)
  !
  sym_4 =.false.
  find_poles=.true.
  nskip = 0
  !
  ! here pade from file
  ! loop over three file
  do i1=1,3
     write(name4,'(I10)') i1
     name1 = 'pol_pade_'//trim(adjustl(name4))
     inquire(file=name1,exist=lexist)
     if(.not.lexist) exit
     open(unit=526,file=name1,form='unformatted')
     !
     nw = 0
     count0=0
     do while(.true.) 
        read(526,end=1) dummyc1, dummyc2
        ! skip some input points
        do i2=1,nskip
           read(526,end=1) dummyc1, dummyc2
        enddo
        if(real(dummyc1).eq.0._8) count0=count0+1
        nw = nw +1
     enddo
1    print *, 'found frequencies in file:', nw
     rewind(526)
     !
     !
     nw_save = nw
     ! perform pade with 5 different conditions:
     ! 1:  plain pade
     ! 2:  using time reversal symmetry, i.e. eps(-w) = eps^*(w)
     ! 3:  using conjugation symmetry, i.e. eps(w^*) = eps^*(w) and eps((-w^*))=eps(w)
     ! 4:  use pade polynomial with only even powers (this is equivalent to 2)
     ! 5:  use pade polynomial with only even powers and conjugation symmetry
     do task=1,5
        nw1=nw_save
        !
        write(name2,'(I10)') task
        !
        if(task.eq.1) then
           nw = nw1
        elseif(task.eq.2) then
           if(count0.gt.1) then
              print *,  'too many real zeroes cycle task 2'
              cycle
           endif
           nw = nw1*2-1
        elseif(task.eq.3) then
           if(count0.gt.1) then
              print *,  'too many real zeroes cycle task 3'
              cycle
           endif
           nw = nw1*4-3
        elseif(task.eq.4) then
           nw = nw1
        elseif(task.eq.5) then
           nw = nw1*2-1
        endif
        !
        allocate(pade_w(nw))
        allocate(u(nw))
        allocate(a_pade(nw))
        !
        ! read in data
        do iw=1, nw1
           read(526) win,pin
           do i2=1,nskip
              read(526,end=1) dummyc1, dummyc2
           enddo
           pade_w(iw) = cmplx(real(win),aimag(win),16)
           u(iw)      = cmplx(real(pin),aimag(pin),16)
           if(task.eq.4.or.task.eq.5)  pade_w(iw) =  pade_w(iw)**2
           !
!           if(iw.eq.1) cycle
           if(real(win).eq.0._8) cycle
           !
           if(task.eq.2.or.task.eq.3) then
              pade_w(nw1+iw-1)  = -cmplx(real(win), aimag(win),16)  ! w -> -w
              u(nw1+iw-1)  =  cmplx(real(pin),-aimag(pin),16)       ! u -> u*
           endif
           !
           if(task.eq.3) then
              pade_w(2*nw1+iw-2)  =  cmplx(real(win),-aimag(win),16)   ! w -> w*
              u(2*nw1+iw-2)  =  cmplx(real(pin), -aimag(pin),16)       ! u -> u*
              !
              pade_w(3*nw1+iw-3)  = -cmplx(real(win),-aimag(win),16)       ! w -> -w*
              u(3*nw1+iw-3)  =  cmplx(real(pin),aimag(pin),16)        ! u -> u
           endif
           !
           if(task.eq.5) then
              pade_w(nw1+iw-1)  =  cmplx(real(win),-aimag(win),16)   ! w -> w*
              u(nw1+iw-1)  =  cmplx(real(pin), -aimag(pin),16)       ! u -> u*
              pade_w(nw1+iw-1)  = pade_w(nw1+iw-1)**2
           endif
           !
        enddo
        rewind(526)
        !
        call pade_coeff( nw, pade_w, u, a_pade)
        !
        name1 = 'pol_'//trim(adjustl(name4))//'_sym_'//trim(adjustl(name2))
        open(unit=123,file=name1)
        ! evaluate functiona real freuquencies
        do iw=1,3001
           w = cmplx((iw-1)*0.01,0.3)!/13.605698066
           !
           if(task.le.3) then
              call pade_eval( nw, pade_w, a_pade,w, padapp )
           else
              call pade_eval( nw, pade_w, a_pade,w**2, padapp )
           endif
           !
           write(123,'(e12.6,2x,e12.6,2x,e12.6)') real(w) , real(padapp) , aimag(padapp)
           !
        enddo
        !
        ! epilogue 
        write(123,*) '#  '
        write(123,*) '# used Pade frequencies in [eV]'
        do iw=1,nw
           write(123,'(A,e12.6,2x,e12.6)') '#', real(pade_w(iw)), aimag(pade_w(iw))
        enddo
        write(123,*) '# used Pade points'
        do iw=1,nw
           write(123,'(A,e12.6,2x,e12.6)') '#', real(u(iw)), aimag(u(iw))
        enddo
        close(123)
        !
        !    if(.not.find_poles) then
        !       deallocate(pade_w)
        !       deallocate(u)
        !       deallocate(a_pade)
        !       cycle
        !    endif
        !
        ! find poles
        !
        if(mod(nw,2).eq.0) then
           m = nw/2
        else
           m = (nw-1)/2
        endif
        !
        !     allocate(poles(m))
        !     allocate(poly(m+1)) ! number of coefficients is n+1 
        allocate(num(m+1))
        allocate(den(m+1))
        !
        call pade_poly ( nw, pade_w, a_pade, m+1, num, den)
        !
        ! write denominator coefficients to file (to be used my Mathematica)
        name1 = 'den_coeff_'//trim(adjustl(name4))//'_sym_'//trim(adjustl(name2))
        open(unit=123,file=name1)
        do i2=1,m+1
           write(123,*) real(den(i2)), aimag(den(i2))
           ! in case padd with zeroes for absent odd coefficients
           if(task.gt.3) write(123,*)  0. , 0. 
        enddo
        close(123)
        !write numerator coefficients to file (to be used my Mathematica)                                                    
        name1 = 'num_coeff_'//trim(adjustl(name4))//'_sym_'//trim(adjustl(name2))
        open(unit=123,file=name1)
        do i2=1,m+1
           write(123,*) real(num(i2)), aimag(num(i2))
           ! in case padd with zeroes for absent odd coefficients                                                               
           if(task.gt.3) write(123,*)  0. , 0.
        enddo
        close(123)
        !
        name1 = 'alt_pol_'//trim(adjustl(name4))//'_sym_'//trim(adjustl(name2))
        open(unit=123,file=name1)
        do iw=1,3001
           w = cmplx((iw-1)*0.01,0.3)
           if(task.le.3) then
              call poly_eval (m+1,num,den,w,padapp)
           else
              call poly_eval (m+1,num,den,w**2,padapp)
           endif
           write(123,'(e12.6,2x,e12.6,2x,e12.6)') real(w) , real(padapp) , aimag(padapp)
        enddo
        close(123)
        !
        !     call find_root(m+1,den,poles)
        !     !
        !     write(name4,'(I10)') i1
        !     name1 = 'pade_poles_'//trim(adjustl(name4))
        !     open(unit=123,file=name1,status='replace')
        !     !
        !     ! write poles to file
        !     do i2=1,m
        !        !
        !        write(123,*) real(poles(i2)), aimag(poles(i2))
        !        !
        !     enddo
        !     close(123)
        !
        deallocate(pade_w)
        deallocate(u)
        deallocate(a_pade)
        !deallocate(poles)
        deallocate(num)
        deallocate(den)
     enddo! tasks
     !
  enddo! 3 files
  !
  call system('/bin/bash mathwrap.sh')
  !
end program pade_tester

subroutine find_root(n,p,r)
  !
  ! find roots r of polynomial of order n-1 and coefficients p
  ! P(x) = p(1) + p(2)*x + p(3)*x^2 + ... + p(n)*x^(n-1)
  !
  !
  implicit none
  !
  ! Arguments
  integer :: n
  complex(16) :: p(n)
  complex(8) :: r(n-1)
  !
  !internal
  integer :: i , info, dum
  complex(8) :: m(n-1,n-1), dummy, work(10*(n-1)), p1(n), t(4*n)
  real(8)    :: rwork(n-1)
  logical     :: dummyl
  !
!  ! construct matrix of which P is the characteristic polynomial
!  m(:,:) = cmplx(0._8,0._8)
!  !
!  ! first row
!  m(1,n-1) = -p(1)/p(n)
!  ! all other rows
!  do i=2,n-1
!     m(i,n-1) = -p(i)/p(n)
!     m(i,i-1) = 1.
!  enddo
!  !
!  ! diagonalize m
!  !
!  call ZGEES( 'N', 'N', dummy, n-1, m, n-1, dum, r, dummy, &
!              1, work, 10*(n-1), rwork, dummyl, info )
  !
!  ! USE slatec library this depends on files
!  !  cpzero.f, cpevl.f, i1mach.f
!  !
!  ! reverse order of coefficients
!  do i=1,n
!     p1(i) = p(n-i+1)
!  enddo
!  !
!  call CPZERO (n, p1, r, t, 0, rwork)
!  !
!print *, rwork
end subroutine find_root

subroutine invertc(n,A)
!********************************************************************
! Invert complex matrix A of dimension n*n using lapack stuff
! Written by H. Huebener February 2011
!*********************** INPUT **************************************
! integer n        : dimension
!******************** INPUT AND OUTPUT *******************************
! complex    A(n,n) : to be inverted->inverted
!*********************************************************************
!
! Modules
! Argument types and dimensions
      integer       ::  n
      complex(8)   ::  A(n,n)
! Internal variables and arrays
      integer             :: ierror,ipiv(n), lwork
      real(8),pointer    :: work(:)
      real(8)            :: temp


      call zgetrf( n, n, A, n, ipiv, ierror )
      if( ierror.eq.0 ) then
         !workspace query
         call zgetri( n, A, n, ipiv, temp, -1, ierror )
         lwork = temp ! dimension of workspace
         allocate(work(lwork*2))
         call zgetri( n, A, n, ipiv, work, lwork, ierror )
      else
         stop 'Terminating due to failed LU decomp'
      endif
      if (ierror.ne.0) then
         stop 'Terminating due to failed inversion'
      endif
      deallocate(work)
  end subroutine invertc
  !
  !
  subroutine poly_eval (m,num,den,w,padapp)
  !-----------------------------------------------------------
  ! evaluate pade approximation at w based on the polynomial indices
  ! num and den
  !
  ! This subroutine uses the recursive algorithm described in
  ! HJ Vidberg and JW Serene, "Solving the Eliashberg equation
  ! by means of N-point Pade' approximants", J Low Temp Phys
  ! 29, 179 (1977). The notation adopted here is the same as
  ! in the above manuscript.
  !
  ! input
  !
  ! N      - order of the Pade' approximant
  ! z(1:N) - points at which the original function is known
  ! a(1:N) - coefficients of the continued fraction
  ! m      - dimension of coefficient vector
  !
  ! output
  !
  !poly(1:m) - coefficients of the denominator polynomial
  !-----------------------------------------------------------
  !
   implicit none
   integer ::  m
   complex(16) ::  w, padapp
   complex(16) :: num(m),den(m),den1,num1
   integer :: i,ii,j
  !
  !
   num1=0._16
   den1=0._16
  ! its ok to let the sum to m, eventhough the order of the numerator is
  ! m-1 in case of odd pade order, because in that case num(m) = 0 automatically
   do i=1,m
      num1  = num1 + w**(i-1)*num(i)
      den1  = den1 + w**(i-1)*den(i)
   enddo
  !
   padapp = num1/den1
  !
end subroutine poly_eval

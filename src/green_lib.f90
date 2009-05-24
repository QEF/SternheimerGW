  !
  !-----------------------------------------------------------------------
  subroutine haydock ( g2kin, vr, un, a, b, a_term, b_term)
  !-----------------------------------------------------------------------
  !
  use gspace
  use constants
  use parameters
  implicit none
  !
  complex(kind=DP) :: vr (nr)
  real(kind=DP) :: g2kin (ngm)
  integer :: i, j, k, ir
  !
  ! notation as in Haydock - review (pag 227)
  ! e.g.: bnpunp = b_{n+1} * u_{n+1}, hun = H * u_n
  complex(DP) :: unp (ngm), un (ngm), unm (ngm), hun (ngm), bnpunp (ngm)
  !
  real(DP) :: an, bn, bnp, a(nstep+2), b(nstep+2), bnp2
  complex(DP) :: ZDOTC, DZNRM2
  integer :: ig, istep, nterm 
  complex(DP) :: Z, aux
  real(DP) :: D(nstep+1), E(nstep), WORK, tmp
  integer :: INFO, N
  !
  ! variables for the terminator
  !
  real(DP) :: a_term, b_term
  !

  !
  ! initialize the recursion 
  !
  a = zero
  b = zero
  bn = zero
  unm = czero
  !
  do istep = 0, nstep
    !
    ! basic Haydock scheme [Eqs. (5.2)-(5.13) of SSP]
    !
    call h_psi_c ( un, hun, g2kin, vr) 
    an = DREAL ( ZDOTC (ngm, un, 1, hun, 1) )
    !
    !  bnpunp = hun - an * un - bn * unm
    !  unp = bunpunp / norm(bnpunp)
    !
    !  NOTE it is crytical to enter ZAXPY with complex coefficients
    !  otehrwise unbelievable things may happen...
    ! 
    aux = dcmplx(-an,0)
    call ZAXPY (ngm, aux, un , 1, hun, 1)
    aux = dcmplx(-bn,0)
    call ZAXPY (ngm, aux, unm, 1, hun, 1)
    bnp = DREAL ( ZDOTC (ngm, hun, 1, hun, 1) )
    bnp = sqrt (bnp)
    call DSCAL (2 * ngm, one / bnp, hun, 1)
    !
    call ZCOPY (ngm, un , 1, unm, 1)  ! unm <- un
    call ZCOPY (ngm, hun, 1, un , 1)  ! un  <- unp
    !
    a (istep+1) = an
    b (istep+2) = bnp 
    bn = bnp
    !
    if (istep.gt.0) then
     !
     ! define the terminator using the extrema of the eigenvalue spectrum
     ! (this should be much more stable than trying to fit the fluctuating
     ! coefficients to constant values) - Giannozzi et al, Appl. Num. Math. 
     ! 4, 273 (1988) - Eq. (7.1)
     !
     ! eigenvalues of the tridiagonal Hamiltonian at this iteration
     !
!    write(6,'(/4x,a)') repeat('-',67)
!    write(6,'(4x,"Eigenvalues (eV) of the tridiagonal matrix at iteration",i5)') &
!      istep + 1
!    write(6,'(4x,a/)') repeat('-',67)
     !
     N = istep + 1
     D = zero
     E = zero
     D(1:N) = a(1:N)
     E(1:N-1) = b(2:N)
     call DSTEV( 'N', N, D, E, Z, 1, WORK, INFO )
!    write(6,'(4x,7f9.3)') D(1:N) * ryd2ev
     !
     a_term = 0.5d0  * ( D(1) + D(N) ) 
     b_term = 0.25d0 * ( D(N) - D(1) ) 
     ! 
!    write(44,*) a_term, b_term
     !
    endif
    !  
  enddo
! !
! ! write down the recursion coefficients
! !  
! write(6,'(/4x,a)') repeat('-',67)
! write(6,'(4x,"Recursion coefficients")') 
! write(6,'(4x,a/)') repeat('-',67)
! do istep = 0, nstep
!   write(6,'(4x,4f15.5)') a(istep+1), b(istep+1), a_term, b_term
! enddo
  !
  return
  end subroutine haydock
  !
  !-----------------------------------------------------------------------
  subroutine recfrac ( w, a, b, a_term, b_term, gr)
  !-----------------------------------------------------------------------
  !
  ! the recursive fraction
  !
  ! Note that the use of the terminator does remove some wild oscillations
  ! at high energy, but it cannot fix in any way the location of the poles.
  ! I have the impression that in order to converge the Green's function
  ! to the energy E, we need to apply the Hamiltonian (# steps in Haydock)
  ! the number of times required to reach the eigenvalue close to E.
  !
  use constants
  use parameters
  implicit none
  !
  complex(DP) :: gr
  real(DP) :: a(nstep+2), b(nstep+2), w, a_term, b_term
  integer :: i
  integer, parameter :: nterm = 300
  ! nterm is the number of fractions used to generate the terminator
  real(DP) :: a1, b2
  !
  gr = czero
  !
  ! the terminator
  !
  do i = nterm+nstep,nstep+1,-1
    gr = w + cmplx(zero, eta) - a_term - b_term ** two * gr
    gr = cone / gr
  enddo
  !
  ! the actual fraction
  !
  do i = nstep,0,-1
    gr = w + cmplx(zero, eta) - a(i+1) - b(i+2) ** two * gr
    gr = cone / gr
  enddo
  !
  ! gr in 1/Ry now
  !
  return
  end subroutine recfrac
  !
  !-----------------------------------------------------------------------
  subroutine normalize ( phi, norm )
  !-----------------------------------------------------------------------
  !
  ! simply normalize a vector in G-space
  !
  use constants
  use parameters
  use gspace
  implicit none
  !
  complex(DP) :: phi(ngm)
  real(DP) :: norm
  complex(DP) :: ZDOTC
  !
  norm = DREAL ( ZDOTC (ngm, phi, 1, phi, 1) )
  norm = sqrt(norm)
  call DSCAL (2 * ngm, one / norm, phi, 1)
  !
  return
  end subroutine normalize
  !

  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine bcgsolve_all (Ax, e, b, x, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_all, anorm, nbnd, g2kin, vr, evq, cw)
  !-----------------------------------------------------------------------
  !
  !   Iterative solution of the linear system:
  !
  !                 ( h - e + w + i * eta ) * x = b         
  !
  !   where h is a complex hermitian matrix, e, w, and eta are
  !   real scalar, x and b are complex vectors
  !
  !   I use the biorthogonal conjugate gradient method as described in
  !   D. A. H. Jacobs, A generalization of the conjugate-gradient method to
  !   solve complex systems, IMA Journal of Numerical Analysis 6, 447 (1986).
  !   Most subsequent relevant CG methods are based on this paper.
  !
  !   Felix Giustino, Berkeley Aug 2008
  !   
  !   As the paper does not treat preconditioning, I followed the general
  !   procedure described by J. R. Shewchuk (An introduction to the
  !   conjugate gradient method without the agonizing pain, Carnegie Mellon 
  !   University, 1994, pp. 39-40) to work out the preconditioned Jacobs algorithm   
  !   (cf. my hand-written notes)
  !   I use the Teter-Payne-Allan preconditioner.
  !
  !   Felix Giustino, Oxford Jan-May 2009
  !
  !   I did the profiling and the most time-consuming part (by a factor 9:1)
  !   is the application of the Hamiltonian, and in particular the FFT forth
  !   and back to do the product V*psi. 
  !
  !   rp  = r preconditioned [ inv(M) * r ]
  !   rt  = r tilda
  !   rtp = r tilda preconditioned [ inv(M) * rt ]
  !
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP, nbnd_occ, eta, maxter
  use gspace, only : nr, ngm
  use constants
  implicit none
  integer :: i
  complex(DP) :: vr(nr), cw
  real(DP) :: g2kin (ngm)
  complex(kind=DP) :: evq (ngm, nbnd_occ)
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             ik      ! input: the k point

  real(kind=DP) :: &
             e(nbnd), &        ! input: the actual eigenvalue
             anorm,   &        ! output: the norm of the error in the solution
             h_diag(ndmx,nbnd), & ! input: preconditioner
             ethr               ! input: the required precision
  complex(kind=DP) :: &
             x (ndmx, nbnd), & ! output: the solution of the linear syst
             y (ndmx, nbnd), & ! output: the solution of the linear syst
             b (ndmx, nbnd), & ! input: the known term
             bt (ndmx, nbnd)    ! input: the known term
  external Ax, &                ! input: the routine computing Ax
           cg_psi               ! input: the routine to apply the preconditioner
  !
  ! local variables
  !
  real(DP) :: eu(nbnd)
  integer :: iter, ibnd, jbnd
  complex(kind=DP), allocatable :: r  (:,:), q  (:,:), p  (:,:), pold  (:,:)
  complex(kind=DP), allocatable :: rt (:,:), qt (:,:), pt (:,:), ptold (:,:)
  complex(kind=DP), allocatable :: rp (:,:), rtp (:,:)
  complex(kind=DP), allocatable :: aux1 (:,:), aux2 (:,:)
  complex(kind=DP) :: a, c, beta, alpha, ZDOTC
  logical :: conv (nbnd), conv_all
!
!  ! timing setup
!  ! 
!  integer :: time_array(8), isteptime
!  real(DP) :: time_old, time_now, steptime(10)
!  call date_and_time(values=time_array)
!  time_old = time_array (5) * 3600 + time_array (6) * 60 + time_array (7) + 0.001 * time_array (8)
!  !
  !
  allocate ( r (ndmx,nbnd), q (ndmx,nbnd), p (ndmx,nbnd), pold (ndmx ,nbnd) )    
  allocate ( rt(ndmx,nbnd), qt(ndmx,nbnd), pt(ndmx,nbnd), ptold(ndmx ,nbnd) )    
  allocate ( rp(ndmx,nbnd), rtp(ndmx,nbnd) )    
  allocate ( aux1 (ndmx,nbnd), aux2 (ndmx,nbnd) )    
  !
  conv = .false.
  conv_all = .false.
  iter = 0
  !

  !
  do while ( iter.lt.maxter .and. .not.conv_all )
    !
    iter = iter + 1
    !
    if (iter .eq. 1) then
       !
       ! r  = b - A*x
       ! rt = conjg ( r )
       !
       call Ax ( x, r, e, ik, nbnd, g2kin, vr, evq, cw)  
       !
       call ZAXPY (ndim*nbnd, -cone, b, 1, r, 1)
       call ZSCAL (ndim*nbnd, -cone, r, 1)
       ! 
       rt = conjg ( r )
       !
       ! p  = inv(M) * r
       ! pt = conjg ( p )
       !
       do ibnd = 1, nbnd
         call ZCOPY (ndim, r  (1, ibnd), 1, p  (1, ibnd), 1)
         call cg_psi (ndmx, ndim, 1, p  (1,ibnd), h_diag(1,ibnd) )
       enddo
       pt = conjg ( p )
       !
    endif
    !
    ! get out if all bands are converged 
    !
    conv_all = .true.
    do ibnd = 1, nbnd
       anorm = sqrt ( abs ( ZDOTC (ndim, r(1,ibnd), 1, r(1,ibnd), 1)  ) )
       if (anorm.lt.ethr) conv(ibnd) = .true.
!      write(6,'(4x,i5,e15.5)') iter, anorm
       conv_all = conv_all .and. conv(ibnd)
    enddo
    if ( conv_all ) exit
    !
    ! Here we calculate q = A * p and qt = A^\dagger * pt
    ! In order to avoid building the Hamiltonian for each band, 
    ! we apply the Hamiltonian to all of the nonconverged bands.
    ! The trick is to pack the nonconverged bands in a aux array,
    ! apply H, and then unpack.
    !
    !  q  = A * p 
    !  qt = A^\dagger * pt
    !
!
! packing does not work: when I do this the bands converge
! until very close to the threshold, and then the unconverged
! ones start to diverge badly. Without the packing the algorithm
! seems very stable and the gain from preconditioning is about
! a factor 3.5 in the number of iterations
!

!    jbnd = 0
!    do ibnd = 1, nbnd
!      if ( .not.conv(ibnd) ) then
!         jbnd = jbnd + 1
!         call ZCOPY (ndim, p  (1, ibnd), 1, aux1 (1, jbnd), 1)
!         call ZCOPY (ndim, pt (1, ibnd), 1, aux2 (1, jbnd), 1)
!         eu (jbnd) = e (ibnd)
!      endif
!    enddo
!    !
!    call Ax ( aux1, q , eu, ik, jbnd, g2kin, vr, evq, + eta )
!    call Ax ( aux2, qt, eu, ik, jbnd, g2kin, vr, evq, - eta )

    !
    ! no packing
    !

    call Ax ( p , q , e, ik, nbnd, g2kin, vr, evq, cw )
    call Ax ( pt, qt, e, ik, nbnd, g2kin, vr, evq, conjg(cw) )

    !
    jbnd = 0
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         !
!        ! remember only q and qt carry the packed jbnd index
!        !
!        jbnd = jbnd + 1
         jbnd = ibnd
         !
         !  rp  = inv(M) * r
         !
         ! M must be real and symmetric in order to decompose as M = transp(E) * E
         ! the TPA preconditioner is real and diagonal, so ok
         !
         call ZCOPY (ndim, r  (1, ibnd), 1, rp  (1, ibnd), 1)
         call cg_psi (ndmx, ndim, 1, rp  (1,ibnd), h_diag(1,ibnd) )
         !
         ! alpha = <rt|rp>/<pt|q>
         ! [ the denominator is stored for subsequent use in beta ]
         !
         a = ZDOTC (ndim, rt(1,ibnd), 1, rp(1,ibnd), 1) 
         c = ZDOTC (ndim, pt(1,ibnd), 1, q (1,jbnd), 1) 
         alpha = a / c 
         !
         !  x  = x  + alpha        * p
         !  r  = r  - alpha        * q 
         !  rt = rt - conjg(alpha) * qt
         !
         call ZAXPY (ndim,  alpha,        p  (1,ibnd), 1, x  (1,ibnd), 1)
         call ZAXPY (ndim, -alpha,        q  (1,jbnd), 1, r  (1,ibnd), 1)
         call ZAXPY (ndim, -conjg(alpha), qt (1,jbnd), 1, rt (1,ibnd), 1)
         !
         !  rp  = inv(M) * r
         !  rtp = inv(M) * rt
         !
         call ZCOPY (ndim, r  (1, ibnd), 1, rp  (1, ibnd), 1)
         call ZCOPY (ndim, rt (1, ibnd), 1, rtp (1, ibnd), 1)
         call cg_psi (ndmx, ndim, 1, rp  (1,ibnd), h_diag(1,ibnd) )
         call cg_psi (ndmx, ndim, 1, rtp (1,ibnd), h_diag(1,ibnd) )
         !
         !  beta = - <qt|rp>/<pt|q>
         !
         a = ZDOTC (ndim, qt(1,jbnd), 1, rp(1,ibnd), 1)
         beta = - a / c 
         !
         ! pold  = p
         ! ptold = pt
         !
         call ZCOPY (ndim, p  (1, ibnd), 1, pold  (1, ibnd), 1)
         call ZCOPY (ndim, pt (1, ibnd), 1, ptold (1, ibnd), 1)
         !
         !  p  = rp  +       beta  * pold
         !  pt = rtp + conjg(beta) * ptold
         !
         call ZCOPY (ndim, rp  (1, ibnd), 1, p  (1, ibnd), 1)
         call ZCOPY (ndim, rtp (1, ibnd), 1, pt (1, ibnd), 1)
         call ZAXPY (ndim,       beta,  pold  (1,ibnd), 1, p (1,ibnd), 1)
         call ZAXPY (ndim, conjg(beta), ptold (1,ibnd), 1, pt(1,ibnd), 1)
         !
       endif
    enddo
    !
  enddo
  if (iter.lt.maxter) then 
!    write(6,'(4x,"bcgsolve_all: ",2i5)') ik, iter
  else
!    write(6,'(4x,"bcgsolve_all: ",2i5,e12.4)') ik, iter, anorm
  endif
  !
  kter = kter + iter
  !
  deallocate ( r , q , p , pold, rt, qt, pt, ptold, aux1, aux2)
  !
  return
  end subroutine bcgsolve_all
  !


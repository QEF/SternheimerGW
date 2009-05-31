  !
  !-----------------------------------------------------------------------
  subroutine bcgsolve_all (Ax, e, b, x, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_all, anorm, nbnd, g2kin, vr, evq)
  !-----------------------------------------------------------------------
  !
  !   iterative solution of the linear system:
  !
  !                 ( h - e + w + i * eta ) * x = b         
  !
  !   where h is a complex hermitian matrix, e, w, and eta are
  !   real scalar, x and b are complex vectors
  !
  !   we use the biorthogonal conjugate gradient method as described in
  !   D. A. H. Jacobs, A generalization of the conjugate-gradient method to
  !   solve complex systems, IMA Journal of Numerical Analysis 6, 447 (1986).
  !
  !   Note that this is the only paper wher I found the algorithm for
  !   complex and non-hermitian operators. Most available versions 
  !   consider the complex and symmetric case, which is not what we need.
  !
  !   Felix Giustino, Berkeley Aug 2008
  !
  !   I did the profiling and the most time-consuming part (by a factor 9:1)
  !   is the application of the Hamiltonian, and in particular the FFT forth
  !   and back to do the product V*psi. We cannot save much on that, better
  !   to look at preconditioning to decrease the number of iterations.
  !
  !   Preconditioner added 5 Jan 2009 - (Teter-Payne-Allan)
  !   In order to add the preconditioner I simply modified the linear
  !   system from Ax = b to (MA)x = Mb with M given by h_diag
  !   and replaced everywhere in the procedure A by MA and b by Mb.
  !   The new matrix MA is still complex, and does not need to be symmetric
  !   or Hermitian in the complex-BCG algorithm by Jacobs - so no problem here.
  !   In practice the preconditioner is applied in only three locations.
  !
  !-----------------------------------------------------------------------
  !
  use parameters, only : DP, nbnd_occ, eta, maxter
  use gspace, only : nr, ngm
  use constants
  implicit none
  integer :: i
  complex(DP) :: vr(nr)
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
       call Ax ( x, r, e, ik, nbnd, g2kin, vr, evq, + eta)  
       !
       call ZAXPY (ndim*nbnd, -cone, b, 1, r, 1)
       call ZSCAL (ndim*nbnd, -cone, r, 1)

       ! TPA preconditioner
       do ibnd = 1, nbnd
         call cg_psi (ndmx, ndim, 1, r(1,ibnd), h_diag(1,ibnd) )
       enddo

       rt = conjg ( r )
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
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         !
         if (iter.eq.1) then
           !
           ! p  = r
           ! pt = rt
           !
           call ZCOPY (ndim, r  (1, ibnd), 1, p  (1, ibnd), 1)
           call ZCOPY (ndim, rt (1, ibnd), 1, pt (1, ibnd), 1)
         else
           !
           !  beta = - <qt|r>/<pt|q>
           !
           a = ZDOTC (ndim, qt(1,ibnd), 1, r(1,ibnd), 1)
           c = ZDOTC (ndim, pt(1,ibnd), 1, q(1,ibnd), 1)
           beta = - a / c
           !
           !  p  = r  +       beta  * pold
           !  pt = rt + conjg(beta) * ptold
           !
           call ZCOPY (ndim, r  (1, ibnd), 1, p  (1, ibnd), 1)
           call ZAXPY (ndim,       beta,  pold  (1,ibnd), 1, p (1,ibnd), 1)
           call ZCOPY (ndim, rt (1, ibnd), 1, pt (1, ibnd), 1)
           call ZAXPY (ndim, conjg(beta), ptold (1,ibnd), 1, pt(1,ibnd), 1)
           !
         endif
         !
       endif
    enddo
    !

    !
    ! Here we calculate q = A * p and qt = A^\dagger * pt
    ! In order to avoid building the Hamiltonian for each band, 
    ! we apply the Hamiltonian to all of the nonconverged bands.
    ! The trick is to pack the nonconverged bands in a aux array,
    ! apply H, and then unpack.
    !
    !  q  = A * p 
    !
    jbnd = 0
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         jbnd = jbnd + 1
         call ZCOPY (ndim, p  (1, ibnd), 1, aux1  (1, jbnd), 1)
         eu (jbnd) = e (ibnd)
       endif
    enddo

    call Ax ( aux1, aux2, eu, ik, jbnd, g2kin, vr, evq, + eta )

    ! TPA preconditioner
    do ibnd = 1, jbnd
      call cg_psi (ndmx, ndim, 1, aux2(1,ibnd), h_diag(1,ibnd) )
    enddo

    jbnd = 0
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         jbnd = jbnd + 1
         call ZCOPY (ndim, aux2 (1, jbnd), 1, q  (1, ibnd), 1)
       endif
    enddo

    !
    !  qt = A^\dagger * pt
    !
    jbnd = 0
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         jbnd = jbnd + 1
         call ZCOPY (ndim, pt (1, ibnd), 1, aux1 (1, jbnd), 1)
         eu (jbnd) = e (ibnd)
       endif
    enddo

    ! TPA preconditioner
    do ibnd = 1, jbnd
      call cg_psi (ndmx, ndim, 1, aux1(1,ibnd), h_diag(1,ibnd) )
    enddo

    call Ax ( aux1, aux2, eu, ik, jbnd, g2kin, vr, evq, - eta )

    jbnd = 0
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         jbnd = jbnd + 1
         call ZCOPY (ndim, aux2 (1, jbnd), 1, qt (1, ibnd), 1)
       endif
    enddo
    !

    !
    do ibnd = 1, nbnd
       if ( .not.conv(ibnd) ) then
         !
         ! alpha = <rt|r>/<pt|q>
         !
         a = ZDOTC (ndim, rt(1,ibnd), 1, r(1,ibnd), 1)
         c = ZDOTC (ndim, pt(1,ibnd), 1, q(1,ibnd), 1)
         alpha = a / c
         !
         !  x  = x  + alpha        * p
         !  r  = r  - alpha        * q 
         !  rt = rt - conjg(alpha) * qt
         !
         call ZAXPY (ndim,  alpha,        p  (1,ibnd), 1, x  (1,ibnd), 1)
         call ZAXPY (ndim, -alpha,        q  (1,ibnd), 1, r  (1,ibnd), 1)
         call ZAXPY (ndim, -conjg(alpha), qt (1,ibnd), 1, rt (1,ibnd), 1)
         !
         ! pold  = p
         ! ptold = pt
         !
         call ZCOPY (ndim, p  (1, ibnd), 1, pold  (1, ibnd), 1)
         call ZCOPY (ndim, pt (1, ibnd), 1, ptold (1, ibnd), 1)
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
  deallocate ( r , q , p , pold, rt, qt, pt, ptold, aux1, aux2)
  !
  return
  end subroutine bcgsolve_all
  !

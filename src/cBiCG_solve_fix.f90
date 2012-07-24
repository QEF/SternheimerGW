SUBROUTINE cbcg_solve_fix(h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, cw, tprec)
!
!-----------------------------------------------------------------------
!
!   Iterative solution of the linear system:
!
!                 ( h - e + w + i * eta ) * x = b
!                 ( h - cw + i * eta ) * G(G,G') = -\delta(G,G')
!
!   where h is a complex hermitian matrix, e, w, and eta are
!   real scalar, x and b are complex vectors

USE kinds,       ONLY: DP
USE mp_global,   ONLY: intra_pool_comm, mpime
USE mp,          ONLY: mp_sum
USE control_gw,  ONLY: maxter_green

implicit none

! first I/O variables

logical :: tprec

integer ::   ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             npol, & ! input: number of components of the wavefunctions
             ik      ! input: the k point

real(DP) :: &
             anorm,   & ! output: the norm of the error in the solution
             ethr,    & ! input: the required precision
             h_diag(ndmx,nbnd) ! input: an estimate of ( H - \epsilon )

!COMPLEX(DP) :: h_diag(ndmx*npol,nbnd) ! input: an estimate of ( H - \epsilon )

  complex(DP) :: &
             dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx*npol, nbnd)   ! input: the known term

  logical :: conv_root ! output: if true the root is converged

  external h_psi       ! input: the routine computing h_psi

  external cg_psi      ! input: the routine computing cg_psi

  !
  !  here the local variables
  !

!HL upping iterations to get convergence with green_linsys?

  !integer, parameter :: maxter = 200
  !integer, parameter :: maxter = 600
  !the maximum number of iterations
  integer :: iter, ibnd, lbnd
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged
  !HL NB GWTC: g t h hold -> SGW: r q p pold
  complex(DP), allocatable :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space

  COMPLEX(DP)  :: cw

!HL need to introduce gt tt ht htold for BICON
! also gp grp for preconditioned systems

  complex(DP), allocatable :: gt (:,:), tt (:,:), ht (:,:), htold (:,:)
  complex(DP), allocatable :: gp (:,:), gtp (:,:)
  complex(DP) ::  dcgamma, dclambda, alpha, beta
  !  the ratio between rho
  !  step length
  complex(DP), external :: zdotc
  !HL (eigenvalue + iw) 

  complex(DP) :: e(nbnd)

  ! the scalar product
  real(DP), allocatable :: rho (:), rhoold (:), a(:), c(:), astar(:), cstar(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  
!HL 
      call start_clock ('cgsolve')
!     call start_clock ('cbcgsolve')

  allocate ( g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol ,nbnd) )
  allocate ( gt(ndmx*npol,nbnd), tt(ndmx*npol,nbnd), ht(ndmx*npol,nbnd), &
             htold(ndmx*npol, nbnd) )
  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))
  allocate (rho(nbnd),rhoold(nbnd))

 ! WRITE(6,*) g,t,h,hold
 ! Initialize

  kter_eff = 0.d0

  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
 
  conv_root = .false.

  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)

  gt = (0.d0,0.d0)
  tt = (0.d0,0.d0)
  ht = (0.d0,0.d0)
  htold = (0.d0,0.d0)
  gp(:,:) = (0.d0, 0.0d0)
  gtp(:,:) = (0.d0, 0.0d0)

!HL PH.x does explicit do loop: 
!I'd prefer to mimic SGW approach here, but Quesspro Guys like there goto etc
!and I want to stick with their convergence parameters/mindset etc. 

  do iter = 1, maxter_green
    ! kter = kter + 1
    ! g    = (-PcDv\Psi) - (H \Delta\Psi)
    ! gt   = conjg( g)
    ! r    = b - Ax 
    ! rt   = conjg ( r )
    ! write(6, '("cBiCG")')
     if (iter .eq. 1) then
        !r = b - A* x
        !rt = conjg (r) 
        call h_psi (ndim, dpsi, g, e, cw, ik, nbnd)
       !call h_psi (ndim, dpsi, g, e, ik, nbnd)
        do ibnd = 1, nbnd
           call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
           call zscal (ndim, (-1.0d0, 0.0d0), g(1,ibnd), 1)
           gt(:,ibnd) = conjg ( g(:,ibnd) )
        ! p   =  inv(M) * r
        ! pt  =  conjg ( p )
           call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
           if(tprec) call cg_psi(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
           ht(:,ibnd) = conjg( h(:,ibnd) )
        enddo
        IF (npol==2) THEN
           do ibnd = 1, nbnd
              call zaxpy (ndim, (-1.d0,0.d0), d0psi(ndmx+1,ibnd), 1, &
                                              g(ndmx+1,ibnd), 1)
              gt(:,ibnd) = conjg ( g(:,ibnd) )
           enddo
        END IF
     endif

! I've moved:  compute preconditioned residual vector and convergence check
! into the first iteration loop. Like in SGW.
! p  = inv(M) * r
! pt = conjg ( p )
!    lbnd = 0
!     do ibnd = 1, nbnd
!        if (conv (ibnd) .eq.0) then
!           lbnd = lbnd+1
!           call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
!           call cg_psi(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
!           ht(:, ibnd) = conjg (h(:,ibnd))
!       HL
!           not too sure what this rho is up to
!           rho(lbnd) = zdotc (ndmx*npol, h(1,ibnd), 1, g(1,ibnd), 1)
!        endif
!     enddo

     lbnd = nbnd
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)


#ifdef __PARA
!HL rho   call mp_sum(  rho(1:lbnd) , intra_pool_comm )
#endif

!     do ibnd = nbnd, 1, -1
!        write(6,*) ibnd, anorm
!        if (conv(ibnd).eq.0) then
!HL           rho(ibnd)=rho(lbnd)
!HL           lbnd = lbnd -1
!HL           anorm = sqrt (rho (ibnd) )
!HL switching convergence test so that its the same as in SGW
!           anorm = sqrt ( abs ( ZDOTC (ndim, g(1,ibnd), 1, g(1,ibnd), 1)  ) )
!           if (anorm.lt.ethr) conv (ibnd) = 1
!        endif
!     enddo
!
! "get out if all bands are converged."

     conv_root = .true.
     do ibnd = 1, nbnd
        anorm = sqrt ( abs ( ZDOTC (ndim, g(1,ibnd), 1, g(1,ibnd), 1)))
!        write(mpime+600,*) iter, anorm
        if (anorm.lt.ethr) conv (ibnd) = 1
        conv_root = conv_root.and.(conv (ibnd).eq.1)
     enddo

     if (conv_root) goto 100

!  " compute the step direction h. Conjugate it to previous step"
!   HL actually ... We do the conjugation steps with h a leetle bit differently.
!     lbnd = 0
!     do ibnd = 1, nbnd
!        if (conv (ibnd) .eq.0) then
!          change sign to h
!           call dscal (2 * ndmx * npol, - 1.d0, h (1, ibnd), 1)
!           if (iter.ne.1) then
!              dcgamma = rho (ibnd) / rhoold (ibnd)
!              call zaxpy (ndmx*npol, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
!           endif
! here hold is used as auxiliary vector in order to efficiently compute t = A*h
! it is later set to the current (becoming old) value of h
!           lbnd = lbnd+1
!           call zcopy (ndmx*npol, h (1, ibnd), 1, hold (1, lbnd), 1)
!           call zcopy (ndmx*npol, ht (1, ibnd), 1, htold (1, lbnd), 1)
!           eu (lbnd) = e (ibnd)
!        endif
!     enddo
! compute t = A*h
! SGW: Here we calculate q = A * p and qt = A^\dagger * pt
! GWTC: Here we calculate t = A * h and tt = A^\dagger * ht
!****************** THIS IS THE MOST EXPENSIVE PART**********************!


!    write(6,*) h(:,:)

     call h_psi (ndim, h, t, e(1), cw, ik, nbnd)
     call h_psi (ndim, ht, tt, e(1), conjg(cw), ik, nbnd)
      
!     write(6,'("Search Directions")')
!     write(6,*) t(:,:)
!     stop
!     call h_psi (ndim, h, t, e(1), ik, nbnd)
!     call h_psi (ndim, ht, tt, e(1), ik, nbnd)
!     compute the coefficients a and c for the line minimization
!     compute step length lambda
!     lbnd=0
!     do ibnd = 1, nbnd
!        if (conv (ibnd) .eq.0) then
!           lbnd=lbnd+1
!           lbnd = ibnd
!           old regular CG stuff
!           HL a(lbnd) = zdotc (ndmx*npol, h(1,ibnd), 1, g(1,ibnd), 1)
!           HL c(lbnd) = zdotc (ndmx*npol, h(1,ibnd), 1, t(1,lbnd), 1)
!           call zcopy(ndmx*npol, g(1,ibnd), 1, gp(1,ibnd), 1)
!           call cg_psi(ndmx, ndim,1, gp(1,ibnd), h_diag(1,ibnd)) 
!        end if
!     end do

#ifdef __PARA
!HL   call mp_sum(  a(1:lbnd), intra_pool_comm )
!     call mp_sum(  c(1:lbnd), intra_pool_comm )
#endif

     lbnd=0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then

        !HL packing? lbnd=lbnd+1

          lbnd = ibnd

        !HL  dclambda is the PH step size I want alpha...
        !HL  dclambda = CMPLX( - a(lbnd) / c(lbnd), 0.d0,kind=DP)
        !
        ! alpha = <rt|rp>/<pt|q>
        ! [ the denominator is stored for subsequent use in beta ]
        !

         call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
         if (tprec) call cg_psi (ndmx, ndim, 1, gp(1,ibnd), h_diag(1,ibnd) )
         a(ibnd) = ZDOTC (ndim, gt(1,ibnd), 1, gp(1,ibnd), 1)
         c(ibnd) = ZDOTC (ndim, ht(1,ibnd), 1, t (1,lbnd), 1)
         alpha = a(ibnd) / c(ibnd)
        ! write(mpime+600,*) a(ibnd), c(ibnd), alpha

        !
        !  x  = x  + alpha        * p
        !  r  = r  - alpha        * q
        !  rt = rt - conjg(alpha) * qt
        !

         call ZAXPY (ndmx*npol,  alpha,        h  (1,ibnd), 1, dpsi  (1,ibnd), 1)
         call ZAXPY (ndmx*npol, -alpha,        t  (1,lbnd), 1, g  (1,ibnd), 1)
         call ZAXPY (ndmx*npol, -conjg(alpha), tt (1,ibnd), 1, gt (1,ibnd), 1)

        !  rp  = inv(M) * r
        !  rtp = inv(M) * rt

         call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, gt (1, ibnd), 1, gtp (1, ibnd), 1)
        if (tprec) call cg_psi (ndmx, ndmx*npol, 1, gp  (1,ibnd), h_diag(1,ibnd) )
        if (tprec) call cg_psi (ndmx, ndmx*npol, 1, gtp (1,ibnd), h_diag(1,ibnd) )

        !
        !  beta = - <qt|rp>/<pt|q>
        !

         a(ibnd) = ZDOTC (ndmx*npol, tt(1,ibnd), 1, gp(1,ibnd), 1)
         beta = - a(ibnd) / c(ibnd)
        ! write(mpime+600,*) a(ibnd), beta
        !
        ! pold  = p
        ! ptold = pt
        !

         call ZCOPY (ndmx*npol, h  (1, ibnd), 1, hold  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, ht (1, ibnd), 1, htold (1, ibnd), 1)

        !
        !  p  = rp  +       beta  * pold
        !  pt = rtp + conjg(beta) * ptold
        !

         call ZCOPY (ndmx*npol, gp  (1, ibnd), 1, h  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, gtp (1, ibnd), 1, ht (1, ibnd), 1)
         call ZAXPY (ndmx*npol,       beta,  hold  (1,ibnd), 1, h (1,ibnd), 1)
         call ZAXPY (ndmx*npol, conjg(beta), htold (1,ibnd), 1, ht(1,ibnd), 1)

        !  write(6,*)h, ht
        !  move to new position
        !HL call zaxpy (ndmx*npol, dclambda, h(1,ibnd), 1, dpsi(1,ibnd), 1)
        !   call zaxpy (ndmx*npol, alpha, h(1,ibnd), 1, dpsi(1,ibnd), 1)
        !   update to get the gradient
        !   g=g+lam
        !   call zaxpy (ndmx*npol, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
        !   save current (now old) h and rho for later use
        !   call zcopy (ndmx*npol, h(1,ibnd), 1, hold(1,ibnd), 1)
        !   rhoold (ibnd) = rho (ibnd)
        endif
     enddo
  enddo

!write(6,*)d0psi(:,2)
!write(6,*) g(:,2)
!write(6,*) h(:,1)

100 continue
  kter = kter_eff
  deallocate (rho, rhoold)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)
  deallocate (gt, tt, ht, htold)
  deallocate (gtp, gp)
  call stop_clock ('cgsolve')
  return
END SUBROUTINE cbcg_solve_fix
 

SUBROUTINE cbcg_solve_coul(h_psi, cg_psi, e, d0psi, dpsi, dpsic, h_diag, &
           ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, niters, alphabeta,tprec)
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
USE control_gw,  ONLY: maxter_green
USE units_gw,    ONLY: iunresid, lrresid, iunalphabeta, lralphabeta

IMPLICIT NONE
!
! first I/O variables
!
  integer ::   ndmx,  & ! input: the maximum dimension of the vectors
               ndim,  & ! input: the actual dimension of the vectors
               kter,  & ! output: counter on iterations
               nbnd,  & ! input: the number of bands
               npol,  & ! input: number of components of the wavefunctions
               ik,    & ! input: the k point
               nrec   ! for composite rec numbers

  integer :: niters(nbnd)
  real(DP) :: anorm,   &        ! output: the norm of the error in the solution
              ethr,    &        ! input: the required precision
              h_diag(ndmx,nbnd) ! input: an estimate of ( H - \epsilon )
!Frommer paper defines beta as \frac{\rho_{k}}{\rho_{k-1}}
  COMPLEX(DP)    :: beta_old 
  complex(DP)    :: dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
                    d0psi (ndmx*npol, nbnd)   ! input: the known term
  COMPLEX(DP)    :: alphabeta(2, nbnd, maxter_green+1)
  COMPLEX(DP)    :: dpsic(ndmx, nbnd, maxter_green+1)
  logical :: conv_root ! output: if true the root is converged
  external h_psi       ! input: the routine computing h_psi
  external cg_psi      ! input: the routine computing cg_psi
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
  complex(DP), external :: zdotc
  complex(DP) :: e(nbnd), eu(nbnd)
  !the scalar product
  real(DP), allocatable :: rho (:), a(:), c(:)
  real(DP) :: kter_eff
  LOGICAL :: tprec

  call start_clock ('cgsolve')
  allocate ( g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol ,nbnd) )
  allocate ( gt(ndmx*npol,nbnd), tt(ndmx*npol,nbnd), ht(ndmx*npol,nbnd), &
             htold(ndmx*npol, nbnd) )
  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))
  allocate (rho(nbnd))
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

!seed frequency should be cw=0
  cw = (0.0d0, 0.0d0)
  do iter = 1, maxter_green
        ! kter = kter + 1
        ! g    = (-PcDv\Psi) - (H \Delta\Psi)
        ! gt   = conjg( g)
        ! r    = b - Ax 
        ! rt   = conjg ( r )
     if (iter .eq. 1) then
        !r = b - A* x
        !rt = conjg (r) 
        call h_psi (ndim, dpsi, g, e, cw, ik, nbnd)
        do ibnd = 1, nbnd
!initial residual should be r = b
           !call davcio (d0psi(:,1), lrresid, iunresid, iter, +1)
           !dpsic(:, ibnd, iter) = d0psi(:,ibnd)
           call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
           call zscal (ndim, (-1.0d0, 0.0d0), g(1,ibnd), 1)
           if(tprec) call cg2_psi(ndmx, ndim, 1, g(1,ibnd), h_diag(1,ibnd) )
           !dpsic(:, ibnd, iter) = d0psi(:,ibnd)
           dpsic(:, ibnd, iter) = g(:,ibnd)
        ! p   =  inv(M) * r
        ! pt  =  conjg ( p )
           call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
        !  gt(:,ibnd) = conjg (g(:,ibnd) )
        !  ht(:,ibnd) = conjg( h(:,ibnd) )
        !not necessary to choose tilde
           gt(:,ibnd) =  (g(:,ibnd) )
           ht(:,ibnd) =  (h(:,ibnd) )
        enddo
     endif!iter.eq.1

     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd).eq.0) then
            lbnd = lbnd+1
            rho(lbnd) = abs(ZDOTC (ndim, g(1,ibnd), 1, g(1,ibnd), 1))
        endif
     enddo

     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)

     do ibnd = nbnd, 1, -1
        if (conv(ibnd).eq.0) then
            rho(ibnd) = rho(lbnd)
            lbnd = lbnd-1
            anorm = sqrt(rho(ibnd))
            if (anorm.lt.ethr) conv (ibnd) = 1
        !band dependent number of iterations
            niters(ibnd) = iter
        endif
     enddo

     conv_root = .true.
     do ibnd = 1, nbnd
        conv_root = conv_root.and.(conv (ibnd).eq.1)
     enddo

     if (conv_root) goto 100
! compute t = A*h
! we only apply hamiltonian to unconverged bands.
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv(ibnd).eq.0) then
            lbnd = lbnd + 1
            call zcopy(ndmx*npol, h(1,ibnd),  1, hold(1,  lbnd), 1)
            call zcopy(ndmx*npol, ht(1,ibnd), 1, htold(1, lbnd), 1)
            eu(lbnd) = e(ibnd)
        endif
     enddo
!****************** THIS IS THE MOST EXPENSIVE PART**********************!
     if(tprec) then
        do ibnd =1, lbnd
             call ZCOPY (ndmx*npol, hold  (1, ibnd), 1, gp  (1, ibnd), 1)
             call ZCOPY (ndmx*npol, htold (1, ibnd), 1, gtp (1, ibnd), 1)
             call cg2_psi (ndmx, ndim, 1, gp(1,ibnd), h_diag(1,ibnd))
             call cg2_psi (ndmx, ndim, 1, gtp(1,ibnd), h_diag(1,ibnd))
        enddo 
        call h_psi (ndim, gp,   t, eu(1), cw, ik, lbnd)
        call h_psi (ndim, gtp, tt, eu(1), conjg(cw), ik, lbnd)
     else
        call h_psi (ndim, hold,   t, eu(1), cw, ik, lbnd)
        call h_psi (ndim, htold, tt, eu(1), conjg(cw), ik, lbnd)
     endif
     if(tprec) then
       do ibnd =1, lbnd
             call cg2_psi (ndmx, ndim, 1, t(1,ibnd),  h_diag(1,ibnd))
             call cg2_psi (ndmx, ndim, 1, tt(1,ibnd), h_diag(1,ibnd))
       enddo 
     endif

     lbnd=0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
            lbnd=lbnd+1
        !alpha = <\tilde{r}|M^{-1}r>/<\tilde{u}|A{u}>
        ![ the denominator is stored for subsequent use in beta ]
            call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
            a(lbnd) = ZDOTC (ndim, gt(1,ibnd), 1, gp(1,ibnd), 1)
            c(lbnd) = ZDOTC (ndim, ht(1,ibnd), 1, t (1,lbnd), 1)
        endif
    enddo

     lbnd=0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1 
           alpha = a(lbnd)/c(lbnd)
           alphabeta(1,ibnd,iter) = alpha
        !  x  = x  + alpha        * u
         call ZAXPY (ndmx*npol,  alpha,        h  (1,ibnd), 1, dpsi  (1,ibnd), 1)
        !  r  = r  - alpha       * Au
        !  \tilde{r} = \tilde{r} - conjg(alpha) * A^{H}\tilde{u}
         call ZAXPY (ndmx*npol, -alpha,        t  (1,lbnd), 1, g  (1,ibnd), 1)
         call ZAXPY (ndmx*npol, -conjg(alpha), tt (1,lbnd), 1, gt (1,ibnd), 1)
        !rp  = inv(M) * r
        !rtp = inv(M) * rt
         call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, gt (1, ibnd), 1, gtp (1, ibnd), 1)
!        nrec = iter+1
!        call davcio (g(:,1), lrresid, iunresid, nrec, +1)
         dpsic(:, ibnd, iter+1) = g(:,ibnd)
        ! if (mpime.eq.0) then
        !     write(503,*)g(:,1)
        ! endif
         a(lbnd) = ZDOTC (ndmx*npol, tt(1,lbnd), 1, gp(1,ibnd), 1)
         beta = - a(lbnd) / c(lbnd)
!        alphabeta(2) = beta
!        call davcio (alphabeta, lralphabeta, iunalphabeta, iter, +1)
         alphabeta(2,ibnd,iter) = beta
        ! u_{old}  = u
        ! \tilde{u}_{old} = \tilde{u}
         call ZCOPY (ndmx*npol, h  (1, ibnd), 1, hold  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, ht (1, ibnd), 1, htold (1, ibnd), 1)
        !new search directions
        !  u  = M^{-1}r  +  beta  * u_old
        !  \tilde{u} = M^{-1}\tilde{r} + conjg(beta) * \tilde{u}_old
         call ZCOPY (ndmx*npol, gp  (1, ibnd), 1, h  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, gtp (1, ibnd), 1, ht (1, ibnd), 1)
         call ZAXPY (ndmx*npol,       beta,  hold  (1,ibnd), 1, h (1,ibnd), 1)
         call ZAXPY (ndmx*npol, conjg(beta), htold (1,ibnd), 1, ht(1,ibnd), 1)
        endif
     enddo!do ibnd
  enddo!iter

100 continue
  kter   =  kter_eff
  deallocate (rho)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)
  deallocate (gt, tt, ht, htold)
  deallocate (gtp, gp)
  call stop_clock ('cgsolve')
  return
END SUBROUTINE cbcg_solve_coul

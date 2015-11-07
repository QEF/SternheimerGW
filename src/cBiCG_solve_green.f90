!-----------------------------------------------------------------------
! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt.
!-----------------------------------------------------------------------
SUBROUTINE cbcg_solve_green(h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
    ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, cw, niters, tprec)
!-----------------------------------------------------------------------
!
!   Iterative solution of the linear system:
!
!                 ( h - e + w + i * eta ) * x = b
!                 ( h - cw + i * eta ) * G(G,G') = -\delta(G,G')
!
!   where h is a complex hermitian matrix, e, w, and eta are
!   real scalar, x and b are complex vectors
!
  USE kinds,       ONLY : DP
  USE mp_global,   ONLY : intra_pool_comm
  USE mp,          ONLY : mp_sum
  USE control_gw,  ONLY : maxter_green
  USE units_gw,    ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
  USE gwsigma,     ONLY : sigma_x_st, sigma_c_st

  implicit none
!first I/O variables
!Frommer paper defines beta as \frac{\rho_{k}}{\rho_{k-1}}
  complex(DP), allocatable :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  complex(DP)  :: cw
  complex(DP) :: alphabeta(2), beta_old
  complex(DP) :: alpha1, beta1
  complex(DP) :: dpsi  (ndmx*npol, nbnd), & ! output: the solution of the linear syst
                 d0psi (ndmx*npol, nbnd)    ! input: the known term
  complex(DP), allocatable :: gt (:,:), tt (:,:), ht (:,:), htold (:,:)
  complex(DP), allocatable :: gp (:,:), gtp (:,:)
  complex(DP) ::  dcgamma, dclambda, alpha, beta
  complex(DP), external :: zdotc
  complex(DP) :: e(nbnd)

  real(DP),    allocatable :: rho (:), rhoold (:), a(:), c(:), astar(:), cstar(:)
  real(DP) :: kter_eff
  real(DP)  :: anorm,   &        ! output: the norm of the error in the solution
               ethr,    &        ! input: the required precision
               h_diag(ndmx,nbnd) ! input: an estimate of ( H - \epsilon )

  integer :: iter, ibnd, lbnd
  integer , allocatable :: conv (:)
  integer   :: ndmx,  & ! input: the maximum dimension of the vectors
               ndim,  & ! input: the actual dimension of the vectors
               kter,  & ! output: counter on iterations
               nbnd,  & ! input: the number of bands
               npol,  & ! input: number of components of the wavefunctions
               ik,    & ! input: the k point
               niters,& ! number of iterations for this BiCG min
               nrec     ! for composite rec numbers

  logical  :: tprec
  logical   :: conv_root ! output: if true the root is converged

  external h_psi       ! input: the routine computing h_psi
  external cg_psi      ! input: the routine computing cg_psi
  
  call start_clock ('cgsolve')
  allocate ( g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol ,nbnd) )
  allocate ( gt(ndmx*npol,nbnd), tt(ndmx*npol,nbnd), ht(ndmx*npol,nbnd), &
             htold(ndmx*npol, nbnd) )
  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))
  allocate (rho(nbnd),rhoold(nbnd))
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
  do iter = 1, maxter_green
     if (iter .eq. 1) then
!r = b - A* x
!rt = conjg (r) 
        call h_psi (ndim, dpsi, g, e, cw, ik, nbnd)
        do ibnd = 1, nbnd
!initial residual should be r = b
!          call davcio (d0psi(:,1), lrresid, iunresid, iter, +1)
           call davcio (d0psi(1:(2*sigma_c_st%ngmt),1), lrresid, iunresid, iter, +1)
           call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
           call zscal (ndim, (-1.0d0, 0.0d0), g(1,ibnd), 1)
!p   =  inv(M) * r
!pt  =  conjg ( p )
           call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
!           gt(:,ibnd) = g(:,ibnd) 
!           ht(:,ibnd) = h(:,ibnd)
           call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gt  (1, ibnd), 1)
           call ZCOPY (ndmx*npol, h  (1, ibnd), 1, ht  (1, ibnd), 1)

        enddo
        IF (npol==2) THEN
           do ibnd = 1, nbnd
              call zaxpy (ndim, (-1.d0,0.d0), d0psi(ndmx+1,ibnd), 1, &
                                              g(ndmx+1,ibnd), 1)
              gt(:,ibnd) = dconjg ( g(:,ibnd) )
           enddo
        END IF
     endif!iter.eq.1
     lbnd = nbnd
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
! "get out if all bands are converged."
     conv_root = .true.
     do ibnd = 1, nbnd
        anorm = sqrt ( abs ( ZDOTC (ndim, g(1,ibnd), 1, g(1,ibnd), 1)  ) )
        if (anorm.lt.ethr) conv (ibnd) = 1
        conv_root = conv_root.and.(conv (ibnd).eq.1)
     enddo
     if (conv_root) goto 100
!****************** THIS IS THE MOST EXPENSIVE PART**********************!
     call h_psi (ndim, h, t, e(1), cw, ik, nbnd)
     call h_psi (ndim, ht, tt, e(1), conjg(cw), ik, nbnd)
     lbnd=0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
          lbnd = ibnd
!alpha = <\tilde{r}|M^{-1}r>/<\tilde{u}|A{u}>
         a(ibnd) = ZDOTC (ndim, gt(1,ibnd), 1, g(1,ibnd), 1)
         c(ibnd) = ZDOTC (ndim, ht(1,ibnd), 1, t (1,lbnd), 1)
         alpha = a(ibnd) / c(ibnd)
         alphabeta(1) = alpha
!x  = x  + alpha        * u
!r  = r  - alpha       * Au
!\tilde{r} = \tilde{r} - conjg(alpha) * A^{H}\tilde{u}
         call ZAXPY (ndmx*npol,  alpha,        h  (1,ibnd), 1, dpsi  (1,ibnd), 1)
         call ZAXPY (ndmx*npol, -alpha,        t  (1,lbnd), 1, g  (1,ibnd), 1)
         call ZAXPY (ndmx*npol, -dconjg(alpha), tt (1,ibnd), 1, gt (1,ibnd), 1)
!rp  = inv(M) * r
!rtp = inv(M) * rt
        call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
        call ZCOPY (ndmx*npol, gt (1, ibnd), 1, gtp (1, ibnd), 1)
!Transformed:
         nrec = iter+1
         call davcio (g(1:(2*sigma_c_st%ngmt), 1), lrresid, iunresid, nrec, +1)
         a(ibnd) = ZDOTC (ndmx*npol, tt(1,ibnd), 1, gp(1,ibnd), 1)
         beta = - a(ibnd) / c(ibnd)
         alphabeta(2) = beta
         call davcio (alphabeta, lralphabeta, iunalphabeta, iter, +1)
!u_{old}  = u
!\tilde{u}_{old} = \tilde{u}
         call ZCOPY (ndmx*npol, h  (1, ibnd), 1, hold  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, ht (1, ibnd), 1, htold (1, ibnd), 1)
!new search directions
!u  = M^{-1}r  +  beta  * u_old
!\tilde{u} = M^{-1}\tilde{r} + conjg(beta) * \tilde{u}_old
         call ZCOPY (ndmx*npol, gp  (1, ibnd), 1, h  (1, ibnd), 1)
         call ZCOPY (ndmx*npol, gtp (1, ibnd), 1, ht (1, ibnd), 1)
         call ZAXPY (ndmx*npol,       beta,  hold  (1,ibnd), 1, h (1,ibnd), 1)
         call ZAXPY (ndmx*npol, conjg(beta), htold (1,ibnd), 1, ht(1,ibnd), 1)
        endif
     enddo!do ibnd
  enddo!iter
100 continue
  niters =  iter
  kter   =  kter_eff
  deallocate (rho, rhoold)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)
  deallocate (gt, tt, ht, htold)
  deallocate (gtp, gp)
  call stop_clock ('cgsolve')
  return
END SUBROUTINE cbcg_solve_green
 

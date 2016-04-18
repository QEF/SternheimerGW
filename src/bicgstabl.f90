!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
SUBROUTINE cbicgstabl(h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, cw, lmresloc, tprec)
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
USE mp_global,   ONLY: intra_pool_comm
USE mp,          ONLY: mp_sum
USE control_gw,  ONLY: maxter_green
USE units_gw,    ONLY: iunresid, lrresid, iunalphabeta, lralphabeta
!USE mp_global,   ONLY: inter_pool_comm, intra_pool_comm, mp_global_end, mpime, &
!                        nproc_pool, nproc, me_pool, my_pool_id, npool
!USE mp,          ONLY : mp_barrier, mp_bcast, mp_sum

implicit none

!tprec true or false conditioning?
logical :: tprec

integer ::   ndmx, &   ! input: the maximum dimension of the vectors
             ndim, &   ! input: the actual dimension of the vectors
             kter, &   ! output: counter on iterations
             nbnd, &   ! input: the number of bands
             npol, &   ! input: number of components of the wavefunctions
             ik,   &   ! input: the k point
             niters, &  ! number of iterations for this BiCG min
             nrec, &    ! for composite rec numbers
             lmresloc     ! set number of loops over minimal residuals...

real(DP) :: &
             anorm,   &        ! output: the norm of the error in the solution
             ethr,    &        ! input: the required precision
             h_diag(ndmx,nbnd) ! input: an estimate of ( H - \epsilon )


COMPLEX(DP) :: alphabeta(2), beta_old

complex(DP) :: &
             dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx*npol, nbnd)   ! input: the known term

!GMRES part
  complex(DP) :: tau(lmresloc + 1,lmresloc + 1)
  complex(DP) :: sigma(lmresloc + 1)
  complex(DP) :: gmgamma(lmresloc + 1)
  complex(DP) :: gmgammapp(lmresloc + 1)
  complex(DP) :: gmgammap(lmresloc + 1)
  complex(DP) :: gschmidt

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
  integer :: iter, ibnd, lbnd, iterj, iteri
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

  complex(DP), allocatable :: gt (:,:), hn (:,:) 
  complex(DP), allocatable :: gp (:,:), gtp (:,:)
  complex(DP) ::  dcgamma, dclambda, alpha, beta, omega
  !  the ratio between rho
  !  step length
  complex(DP), external :: zdotc
  !HL (eigenvalue + iw) 

  complex(DP) :: e(nbnd)
  complex(DP) :: rho, rhoold


!min residual vector
!  complex(DP), allocatable :: s (:,:)

  ! the scalar product
  ! why is this real??
  !real(DP), allocatable :: a(:), c(:)
    complex(DP), allocatable :: a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  call start_clock ('cgsolve')
  allocate ( g(ndmx*npol,lmresloc+1), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol, nbnd) )

  allocate ( gt(ndmx*npol,nbnd), hn(ndmx*npol,nbnd))
             
  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))

  kter_eff = 0.d0

  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
  conv_root = .false.

  g        = dcmplx(0.d0,0.d0)
  t        = dcmplx(0.d0,0.d0)
  h        = dcmplx(0.d0,0.d0)
  hold     = dcmplx(0.d0,0.d0)
  gt       = dcmplx(0.d0,0.d0)
  hn       = dcmplx(0.d0,0.d0)
  gp(:,:)  = dcmplx(0.d0,0.0d0)
  gtp(:,:) = dcmplx(0.d0,0.0d0)
  rho      = dcmplx(1.0d0,0.0d0)
  rhoold   = dcmplx(1.0d0,0.0d0)
  omega    = dcmplx(1.0d0,0.0d0)
  alpha    = dcmplx(0.0d0,0.0d0)
  sigma(:) = dcmplx(1.0d0, 0.0d0)

  do iter = 1, maxter_green
  ! r    = b - Ax 
  ! rt   = conjg ( r )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!cBICG PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! iterj allows for the bicgstab(2,3,4,...) generalizations.
  ! a choice of lmresloc = 1 corresponds to the standard bicgstab.
  !H. A. Van Der Corst Siam J. Sci. Stat. Comput.
  !Vol. 13, No. 2, pp. 631-644
            if (iter .eq. 1) then
               !initialize iter 0 stuff
               !r = b - A* x
               !rt = conjg (r) 
               !call h_psi (ndim, dpsi, g, e, cw, ik, nbnd)
               do ibnd = 1, nbnd
               !initial residual should be r = b
               ! I/O for multishift.
               ! call davcio (d0psi(:,1), lrresid, iunresid, iter, +1)
                  call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
                  call zscal (ndim, (-1.0d0, 0.0d0), g(1,ibnd), 1)
    ! copy r -> u, i.e. u = r_{0}
    ! call zcopy (ndim, g (1, ibnd), 1, h (1, ibnd), 1)
    ! set \tilde{r} = r^{*}
    ! gt(:,ibnd)    = conjg ( g(:,ibnd) )
                  gt(:,ibnd) = g(:,ibnd)
    !for preconditioning we solve the system MAx = Mb
    !x shouldn't change...
    !            gt = Mb
    !need to precondition both sides....
                  if(tprec) call cg_psi(ndmx, ndim, 1, g(1,ibnd), h_diag(1,ibnd) )
                  if(tprec) call cg_psi(ndmx, ndim, 1, gt(1,ibnd), h_diag(1,ibnd) )

    !              write(600+mpime, *) g(:,ibnd)
    !              write(600+mpime, '("finished intialization.")')
               enddo
               IF (npol==2) THEN
                  do ibnd = 1, nbnd
                     call zaxpy (ndim, (-1.d0,0.d0), d0psi(ndmx+1,ibnd), 1, &
                                                     g(ndmx+1,ibnd), 1)
                     gt(:,ibnd) = g(:,ibnd)
                  enddo
               END IF
            endif
  
     rhoold  = -omega*rhoold
     do iterj = 1, lmresloc
        do ibnd = 1, nbnd
                rho     = ZDOTC (ndim, gt(1,ibnd), 1, g(1,ibnd), 1)
                beta    = alpha*(rho/rhoold)
                rhoold  = rho
                do iteri = 1, iterj 
          !\hat{u}_{i} =  \hat{r}_{i} - \beta \hat{u}_{i}
                   hold(:,:) = dcmplx(0.0d0, 0.d0)
                   call ZCOPY (ndmx*npol,  g(1, iteri),      1, hold (1, iteri), 1)
                   call ZAXPY (ndmx*npol,  -beta, h(1,iteri), 1, hold(1,iteri), 1)
                   call ZCOPY (ndmx*npol,  hold(1, iteri),   1, h (1, iteri), 1)
                enddo
            lbnd = nbnd
            kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
        enddo 

      !"get out if all bands are converged."
            conv_root = .true.
            do ibnd = 1, nbnd
               anorm = sqrt ( abs ( ZDOTC (ndim, g(1,ibnd), 1, g(1,ibnd), 1)  ) )
               !write(600+mpime, '("norm")')
               !write(600+mpime,*) anorm, ibnd, ethr
               if (anorm.lt.ethr) conv (ibnd) = 1
               conv_root = conv_root.and.(conv (ibnd).eq.1)
            enddo

            if (conv_root) goto 100
     !****************** THIS IS THE MOST EXPENSIVE PART**********************!
     ! HLstab
     ! t = u_{j+1} = A u_{j}
     ! write(600+mpime, '("applying H to search direction.")')
        do ibnd = 1, nbnd
     ! MA x
           call h_psi (ndim, h(:,ibnd), hn(:,ibnd), e(1), cw, ik, nbnd)
           if(tprec) call cg_psi(ndmx, ndim, 1, hn(1,ibnd), h_diag(1,ibnd))
        enddo

            lbnd=0
            do ibnd = 1, nbnd
               if (conv (ibnd) .eq.0) then
                   lbnd = ibnd
        ! alpha = <\tilde{r}|M^{-1}r>/<\tilde{u}|A{u}>
        ! the denominator is stored for subsequent use in beta
        ! HLSTAB in sleijpen notation
        ! gamma = (u_{j+1}, r_{0}) = (t,gp)
        ! call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
        ! a(ibnd) = ZDOTC (ndim, g(1,ibnd), 1, gt(1,ibnd), 1)
                a(ibnd) = rho
        ! c(ibnd) = ZDOTC (ndim, ht(1,ibnd), 1, t (1,lbnd), 1)
        ! HL-STAB slight difference here
        ! c(ibnd) = ZDOTC (ndim, hn(1,ibnd), 1, gt(1,ibnd), 1)
                c(ibnd) = ZDOTC (ndim, gt(1,ibnd), 1, hn(1,ibnd), 1)
        !do i need to take straight dot products here?
        !c(ibnd) = ZDOTC (ndim, hn(1,ibnd), 1, gt(1,ibnd), 1)
                alpha = a(ibnd) / c(ibnd)
        !HLSTAB do i = 0; iterj:
        !          ri  = ri  - alpha*hn
        !       end
               do iteri = 1, iterj
                  call ZAXPY (ndmx*npol, -alpha, hn(1,ibnd), 1, g(1,ibnd), 1)
               enddo
!!!FOR MINIMAL RESIDUAL SCHEMA WE NEED TO UPDATE THE RESIDUALS AND THE SOLUTION VECTORS here
!!!but we need to hold off on the solution vectors u0
!!!Also we introduce the min res part:
           !t  =  r_{j+1} = A r_{j}
           !t = Hg = Ar i.e. the min resid. part
            call h_psi (ndim, g(:,iterj), t(:,iterj), e(1), cw, ik, nbnd)
           !do we only want to precondtion the real part?
            if(tprec) call cg_psi(ndmx, ndim, 1, t(1,ibnd), h_diag(1,ibnd))

           !MA r
           !x_{0}  = x_{0}  + alpha * u_{0}
               call ZAXPY (ndmx*npol, alpha, h(1,ibnd), 1, dpsi (1,ibnd), 1)
               endif
            enddo!do ibnd
     enddo!iterj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Min Residual Part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do ibnd = 1, nbnd
     do iterj = 2, lmresloc + 1
        if (lmresloc.gt.1)  then
            do iteri = 2, iterj
               tau(iteri, iterj) = (1.0d0/sigma(iteri))*ZDOTC(ndim, g(1,iterj), g(1,iteri), 1)
!orthogonalizes each residual:
               call ZAXPY (ndmx*npol, -tau(iteri, iterj),  g(1,iteri), 1, g(1,iterj), 1)
            enddo
        endif

!should correspond to van der vorst's:sleipjen:SGW
!w=(t,s)(t,t):(\sigma_{1} = (\r_{1},\r_{1})), \gammap=(1/\sigma_{j})(\r_{0}, \r_{j}): 
!sigma=(t,t) gmgammap=(1/sigma(j))*(ro,t)
!sigma(iterj)   = ZDOTC(ndmx*npol, g(1,iterj), 1, g(1,iterj), 1)
        sigma(iterj)   = ZDOTC(ndmx*npol, t(1,ibnd), 1, t(1,ibnd), 1) 
!should be r_{1}! i.e. t 
!gmgammap(iterj) = (1.0d0/sigma(iterj)) * ZDOTC(ndim, g(1,ibnd), 1, g(1,iterj), 1)
!       gmgammap(iterj) = (1.0d0/sigma(iterj)) * ZDOTC(ndim, g(1,ibnd), 1, t(1,ibnd), 1)
        gmgammap(iterj) = (1.0d0/sigma(iterj)) * ZDOTC(ndim, t(1,ibnd), 1, g(1,ibnd), 1)
     enddo
   enddo

!w = (t,s)(t,t) = (\sigma_{j} = (\r_{j}, \r_{j})), \gammap = (1/\sigma_{j})(\r_{0}, \r_{0})
     gmgamma(lmresloc+1) = gmgammap(lmresloc+1)
     omega      = gmgamma(lmresloc+1)


!write(600,'(6f12.7)') gmgamma(:), omega
!ok what's goin on here??
     if (lmresloc.gt.1)  then
         do iterj = lmresloc+1, 2
            gschmidt = dcmplx(0.0d0, 0.0d0)
         do iteri = iterj + 1, lmresloc+1
            gschmidt = gschmidt + tau(iterj, iteri)*gmgamma(iteri)
         enddo
            gmgamma(iterj) = gmgammap(iterj) - gschmidt
         enddo
     endif


     if(lmresloc.gt.1) then
        do iterj = 1, lmresloc-1
           do iteri = iterj + 1, lmresloc-1
              gschmidt = gschmidt + tau(iterj, iteri)*gmgamma(iteri + 1)
           enddo
              gmgammapp(iterj) = gmgamma(iterj+1) + gschmidt
        enddo
     endif

!The second update phase:
   do ibnd = 1, nbnd
!      x  = x  + gmgamma(1)  * \hat{r}_{0}
       call ZAXPY (ndim,  (gmgamma(lmresloc+1)), g(1,ibnd),  1, dpsi(1,ibnd), 1)
!     \hat{u}_{0}  = \hat{u}_{0}  +  beta  * u_old
       call ZAXPY (ndim, -(gmgamma(lmresloc+1)), hn(1,ibnd), 1, h(1,ibnd), 1)
!     \hat{r}_{0}  = \hat{r}_{0}  - gammap_{l}*r_{lmresloc}
       call ZAXPY (ndim, -(gmgammap(lmresloc+1)), t(1,ibnd), 1, g(1,ibnd), 1)
      !and again here we get to a point and then stop progressing...
      !call ZAXPY (ndmx*npol, (gmgamma(lmresloc+1)), hn(1,ibnd), 1, h(1,ibnd), 1)
     if(lmresloc.gt.1) then
        do iterj = 2, lmresloc
!      \hat{u_0}  = \hat{u_0} + gmgamma(lmresloc)*u_{j}
       call ZAXPY (ndmx*npol, -(gmgamma(iterj)),   hn  (1,iterj), 1, h (1,iterj), 1)
!      x_{0}  = x_{0} + gmgamma''_{j}  * r_{j}
       call ZAXPY (ndmx*npol,  gmgammapp(iterj), g(1,iterj),  1, dpsi  (1,iterj), 1)
!      r_{0}  = r_{0}  - gamma'_{j}*r_{j}
       call ZAXPY (ndmx*npol, -(gmgammap(iterj)), t(1,iterj), 1, g  (1,ibnd), 1)
        enddo
     endif
   enddo
  enddo!maxter

100 continue

  niters =  iter
  kter   =  iter !kter_eff
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)
  deallocate (gt, hn)
  deallocate (gtp, gp)
  call stop_clock ('cgsolve')
  return
END SUBROUTINE cbicgstabl
 

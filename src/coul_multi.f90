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
SUBROUTINE coul_multi(h_psi, cg_psi, e, d0psi, x_sig, h_diag, w_ryd, nfreq, &
           ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, niters, tprec)
!-----------------------------------------------------------------------
!
!   Iterative solution of the linear system:
!
!                 ( h - e + w + i * eta ) * x = b
!                 ( h - cw + i * eta ) * G(G,G') = -\delta(G,G')
!
!   where h is a complex hermitian matrix, e, w, and eta are
!   real scalar, x and b are complex vectors
USE control_gw,       ONLY : maxter_coul, alpha_pv
USE freq_gw,          ONLY : nfs
USE kinds,            ONLY : DP
USE linear_op_module, ONLY : linear_op
USE qpoint,           ONLY : ikqs
USE units_gw,         ONLY : iunresid, lrresid, iunalphabeta, lralphabeta

IMPLICIT NONE
!
! first I/O variables
!
  integer :: nfreq, iw, iwp
  integer ::   ndmx,  & ! input: the maximum dimension of the vectors
               ndim,  & ! input: the actual dimension of the vectors
               kter,  & ! output: counter on iterations
               nbnd,  & ! input: the number of bands
               npol,  & ! input: number of components of the wavefunctions
               ik,    & ! input: the k point
               nrec     ! for composite rec numbers

  complex(DP), allocatable  :: g (:,:), t (:,:), h (:,:), hold (:,:)
  complex(DP), allocatable  :: gt (:,:), tt (:,:), ht (:,:), htold (:,:)
!  complex(DP), allocatable  :: gp (:,:), gtp (:,:)
  complex(DP), allocatable  :: u_sig (:,:,:), u_sig_old(:,:,:)
  complex(DP)               :: alphabeta(2, nbnd)

  complex(DP)               :: cw, cwpsi
  complex(DP)               :: dcgamma, dclambda, alpha, beta
  complex(DP)               :: e(nbnd), eu(nbnd)
  complex(DP)               :: alpha_old(nbnd), beta_old(nbnd) , beta_sig(nbnd, nfreq), alpha_sig(nbnd, nfreq)
  complex(DP)               :: pi_coeff (nbnd, nfreq), pi_coeff_old (nbnd, nfreq), pi_coeff_new(nbnd, nfreq)
  complex(DP)               :: w_ryd(nfs)
  complex(DP)               :: dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
                               d0psi (ndmx*npol, nbnd)   ! input: the known term
  complex(DP)               :: x_sig (ndmx*npol, nbnd, nfreq)
  complex(DP), parameter    :: cone = (1.0d0,0.0d0), czero=(0.0d0, 0.0d0)
  complex(DP), external     :: zdotc
  complex(DP)               :: freq(nfreq)
  complex(DP)               :: omega(nbnd)

  real(DP)                  :: anorm,   &               ! output: the norm of the error in the solution
                               ethr,    &               ! input: the required precision
                               h_diag(ndmx,nbnd)        ! input: an estimate of ( H - \epsilon )

  real(DP), allocatable     :: b(:)
  real(DP), allocatable     :: rho (:), a(:), c(:)
  real(DP)                  :: kter_eff
  integer                   :: iter, ibnd, lbnd, ngvecs
  integer , allocatable     :: conv (:)
  integer                   :: niters(nbnd)
  integer                   :: ios, itol

  logical                   :: tprec
  logical                   :: conv_root ! output: if true the root is converged

  external h_psi       ! input: the routine computing h_psi
  external cg_psi      ! input: the routine computing cg_psi




  call start_clock ('cgsolve')
  allocate ( g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol ,nbnd) )
  allocate ( gt(ndmx*npol,nbnd), tt(ndmx*npol,nbnd), ht(ndmx*npol,nbnd), &
             htold(ndmx*npol, nbnd) )
!  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))
  allocate (rho(nbnd))
  allocate(u_sig(ndmx, nbnd, nfreq), u_sig_old(ndmx, nbnd, nfreq))
  allocate(b(nbnd))

  itol = 0
  kter_eff = 0.d0
  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo

  conv_root = .false.
  g    = (0.d0,0.d0)
  gt   = (0.d0,0.d0)
  t    = (0.d0,0.d0)
  tt   = (0.d0,0.d0)
  h    = (0.d0,0.d0)
  ht   = (0.d0,0.d0)
  hold = (0.d0,0.d0)
  htold = (0.d0,0.d0)
 ! gp(:,:)  = (0.d0, 0.0d0)
!  gtp(:,:) = (0.d0, 0.0d0)
  pi_coeff_new(:,:) = (0.0d0,0.0d0)

  x_sig(:,:,:) = czero
  u_sig(:,:,:) = czero
  u_sig_old(:,:,:) = czero
  alpha = czero
  beta  = czero

  do iter  = 1, maxter_coul
! kter = kter + 1
! g    = (-PcDv\Psi) - (H \Delta\Psi)
! gt   = conjg( g)
! r    = b - Ax 
! rt   = conjg ( r )
     cwpsi = (0.0d0, 0.0d0)
     if (iter .eq. 1) then
! r = b - A* x
! rt = conjg (r) 
        do ibnd = 1, nbnd
! initial residual should be r = b
! call davcio (d0psi(:,1), lrresid, iunresid, iter, +1)
! dpsic(:, ibnd, iter) = d0psi(:,ibnd)
         do iw = 1, nfreq
            u_sig(:,ibnd,iw) = d0psi(:,ibnd)
         enddo
         call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1) !calculate |b|^2 
         IF(itol==1) b(ibnd) = abs(ZDOTC (ndmx*npol, d0psi(1,ibnd), 1, d0psi(1,ibnd), 1)) 
         call zscal (ndim, (-1.0d0, 0.0d0), g(1,ibnd), 1)
! MERGE
! call zcopy (ndmx*npol, g (1, ibnd), 1, dpsic(1, ibnd, iter), 1)
! p   =  inv(M) * r
! pt  =  conjg ( p )
         call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
         gt(:,ibnd) = conjg (g(:,ibnd) )
         ht(:,ibnd) = conjg( h(:,ibnd) )
! not necessary to choose tilde
!         call zcopy (ndmx*npol, g(1, ibnd), 1, gt(1, ibnd), 1)
!         call zcopy (ndmx*npol, h(1, ibnd), 1, ht(1, ibnd), 1)
        enddo
        alpha_old(:)        = (1.d0,0.0d0)
        beta_old(:)         = (1.d0,0.0d0)
        pi_coeff_old (:,:)  = (1.d0,0.0d0)
        pi_coeff (:,:)      = (1.d0,0.0d0)
     endif !iter.eq.1

     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd).eq.0) then
            lbnd = lbnd+1
            rho(lbnd) = abs(zdotc (ndim, g(1,ibnd), 1, g(1,ibnd), 1))
            if(itol==1) rho(lbnd)=rho(lbnd)/b(ibnd)
        endif
     enddo

     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
     do ibnd = nbnd, 1, -1
        if (conv(ibnd).eq.0) then
            rho(ibnd) = rho(lbnd)
            lbnd = lbnd-1
            anorm = sqrt(rho(ibnd))
            if (anorm.lt.ethr) conv (ibnd) = 1
            niters(ibnd) = iter
        endif
     enddo

     conv_root = .true.
     do ibnd = 1, nbnd
        conv_root = conv_root.and.(conv (ibnd).eq.1)
     enddo
     if (conv_root) goto 100
     !if (iter.eq.maxter_coul .and. .not.conv_root) then
     !   do ibnd=1, nbnd
     !      if(conv(ibnd)/=1)then
     !        WRITE( 6, '(5x,"kpoint",i4," ibnd",i4, &
     !        &               "solve_linter: root not converged ",e10.3)') &
     !        &                ik , ibnd, anorm
     !      end if
     !   end do
     !end if 
! compute t = A*h
! we only apply hamiltonian to unconverged bands.
     lbnd    = 0
     do ibnd = 1, nbnd
        if (conv(ibnd).eq.0) then
            lbnd = lbnd + 1
            call zcopy(ndmx*npol, h(1,ibnd),  1, hold(1,  lbnd), 1)
            call zcopy(ndmx*npol, ht(1,ibnd), 1, htold(1, lbnd), 1)
            eu(lbnd) = e(ibnd)
        endif
     enddo
!****************** THIS IS THE MOST EXPENSIVE PART ****************** !
    omega(:lbnd) = eu(:lbnd) + cwpsi
    CALL linear_op(ikqs(ik), ndim, -omega(:lbnd), alpha_pv, hold(:,1:lbnd), t(:,1:lbnd))
    omega(:lbnd) = eu(:lbnd) + conjg(cwpsi)
    CALL linear_op(ikqs(ik), ndim, -omega(:lbnd), alpha_pv, htold(:,1:lbnd), tt(:,1:lbnd))

    lbnd=0
    do ibnd = 1, nbnd
       if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
        !alpha = <\tilde{r}|M^{-1}r>/<\tilde{u}|A{u}>
        ![ the denominator is stored for subsequent use in beta ]
           a(lbnd) = zdotc (ndmx*npol, gt(1,ibnd), 1, g(1,ibnd), 1)
           c(lbnd) = zdotc (ndmx*npol, ht(1,ibnd), 1, t (1,lbnd), 1)
       endif
    enddo
    lbnd=0
    do ibnd = 1, nbnd
       if (conv (ibnd).eq.0) then
           lbnd=lbnd+1 
           alpha = a(lbnd)/c(lbnd)
           do iw = 1, nfreq
                   if (iw.le.int((nfreq+1)/2)) then
                       cw = w_ryd(iw)
                   else
                       cw = (-1.0d0,0.0d0)*w_ryd(iw-nfs+1)
                   endif
                   pi_coeff_new(ibnd, iw) = ( cone - alpha*cw)*pi_coeff(ibnd, iw) -      &
                                           (( alpha*beta_old(ibnd))/(alpha_old(ibnd))) *  &
                                            ( pi_coeff_old(ibnd, iw)-pi_coeff(ibnd, iw))
! beta  = (pi_old/pi)**2 *beta, alpha = (pi/pi_new)*alpha
                   alpha_sig(ibnd, iw)    = (pi_coeff(ibnd, iw)/pi_coeff_new(ibnd, iw))*alpha
! x_sig =  x_sig + alpha_sig*u_sig
                   x_sig(:, ibnd, iw)     = x_sig(:,ibnd,iw) + alpha_sig(ibnd,iw)*u_sig(:,ibnd,iw)
           enddo
! x  = x  + alpha        * u
! call zaxpy (ndmx*npol,  alpha,        h  (1,ibnd), 1, dpsi  (1,ibnd), 1)
! if(tprec) call cg3_psi (ndmx, ndim, 1, dpsi(1,ibnd), h_diag(1,ibnd))
! r  = r  - alpha       * Au
! \tilde{r} = \tilde{r} - conjg(alpha) * A^{H}\tilde{u}
           call zaxpy (ndmx*npol, -alpha,        t  (1,lbnd), 1, g  (1,ibnd), 1)
           call zaxpy (ndmx*npol, -conjg(alpha), tt (1,lbnd), 1, gt (1,ibnd), 1)
! rp  = inv(M) * r
! rtp = inv(M) * rt
!           call zcopy (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
!           call zcopy (ndmx*npol, gt (1, ibnd), 1, gtp (1, ibnd), 1)
! MERGE
! call zcopy (ndmx*npol, g (1, ibnd), 1, dpsic(1, ibnd, iter+1), 1)
           a(lbnd) = zdotc (ndmx*npol, tt(1,lbnd), 1, g(1,ibnd), 1)
           beta = - a(lbnd) / c(lbnd)
! MERGE
! alphabeta(2,ibnd,iter) = beta
           do iw = 1, nfreq
              beta_sig(ibnd,iw)  = (pi_coeff(ibnd,iw) / pi_coeff_new(ibnd,iw))**2.d0 * beta
!update the residual and solution vector at each shifted frequency.
              call zcopy (ndim,  u_sig(1,ibnd,iw), 1, u_sig_old(1,ibnd,iw), 1)
              call zcopy (ndim,  g(1,ibnd), 1, u_sig(1,ibnd,iw), 1)
              call zscal (ndim, (cone/pi_coeff_new(ibnd, iw)), u_sig(1,ibnd,iw), 1)
              call zaxpy (ndim,  beta_sig(ibnd,iw), u_sig_old (1,ibnd,iw), 1, u_sig (1,ibnd,iw), 1)
!store values for next iteration
              alpha_old    (ibnd)      = alpha
              beta_old     (ibnd)      = beta
              pi_coeff_old (ibnd, iw)  = pi_coeff(ibnd, iw)
              pi_coeff     (ibnd, iw)  = pi_coeff_new(ibnd, iw)
           enddo !iw
! u_{old}  = u
! \tilde{u}_{old} = \tilde{u}
           call zcopy (ndmx*npol, h  (1, ibnd), 1, hold  (1, ibnd), 1)
           call zcopy (ndmx*npol, ht (1, ibnd), 1, htold (1, ibnd), 1)
! new search directions
! u  = M^{-1}r  +  beta  * u_old
! \tilde{u} = M^{-1}\tilde{r} + conjg(beta) * \tilde{u}_old
           call zcopy (ndmx*npol, g  (1, ibnd), 1, h  (1, ibnd), 1)
           call zcopy (ndmx*npol, gt (1, ibnd), 1, ht (1, ibnd), 1)
           call zaxpy (ndmx*npol,       beta,  hold  (1,ibnd), 1, h (1,ibnd), 1)
           call zaxpy (ndmx*npol, conjg(beta), htold (1,ibnd), 1, ht(1,ibnd), 1)
        endif
     enddo !do ibnd
  enddo !iter

  ! abort if solver doesn't converge
  CALL errore(__FILE__, "Coulomb solver did not converge within given number&
                       & of iterations", maxter_coul)

100 continue
  kter   =  kter_eff
  deallocate (rho)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)
  deallocate (gt, tt, ht, htold)
!  deallocate (gtp, gp)
  deallocate(u_sig, u_sig_old)

  call stop_clock ('cgsolve')
  return
END SUBROUTINE coul_multi


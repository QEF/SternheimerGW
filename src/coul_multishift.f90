SUBROUTINE coul_multishift(ndmx, ndim, nfreq, niters, x_sig, dpsic, alphabeta, freq)
   USE kinds,       ONLY : DP
   USE wvfct,       ONLY : nbnd
   USE units_gw,    ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
   USE freq_gw,     ONLY : fpol, fiu, nwgreen, wgreen
   USE constants,   ONLY : degspin, pi, tpi, RYTOEV, eps8
   USE mp_global,   ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime, &
                           nproc_pool, nproc, me_pool, my_pool_id, npool
   USE control_gw,  ONLY : eta, tr2_green, maxter_green

IMPLICIT NONE
!coefficient of quadratic form
  INTEGER :: nfreq, iw, iwp
  COMPLEX(DP)   :: alpha(nbnd), beta(nbnd)
  complex(kind=DP), allocatable :: u_sig (:,:,:), u_sig_old(:,:,:), r(:,:)
  complex(kind=DP) :: alpha_old(nbnd), beta_old(nbnd) , beta_sig(nbnd, nfreq), alpha_sig(nbnd, nfreq)
!pi coefficients for each frequency up to nfreqgreen.
  complex(kind=DP) :: pi_coeff (nbnd, nfreq), pi_coeff_old (nbnd, nfreq), pi_coeff_new(nbnd, nfreq)
  real(DP) :: w_ryd(nfreq)

REAL(DP) :: anorm(nwgreen)

!variable for reading in the stored alpha beta coefficients.
  COMPLEX(DP)                 :: alphabeta(2, nbnd, maxter_green+1)
  COMPLEX(DP)                 :: dpsic(ndmx, nbnd, maxter_green+1)
  COMPLEX(DP)                 :: x_sig (ndmx, nbnd, nfreq)
  COMPLEX(DP), PARAMETER      :: cone = (1.0d0,0.0d0), czero=(0.0d0, 0.0d0)
  complex(DP), external       :: zdotc
  complex(DP)                 :: freq(nfreq)
  integer ::   ndmx, & ! input: the maximum dimension of the vectors
               ndim, & ! input: the actual dimension of the vectors
               ngvecs,&
               iter,&
               nrec,&
               ibnd
  integer :: ios
  integer :: niters(nbnd)

  ALLOCATE(r(ndmx, nbnd))
  ALLOCATE(u_sig(ndmx, nbnd, nfreq), u_sig_old(ndmx, nbnd, nfreq))
! w_ryd(:) = wgreen(:)/RYTOEV
  w_ryd(:) = freq(:)
!Green shifted system.
     u_sig(:,:,:) = czero
     x_sig(:,:,:) = czero
     r(:,:)       = czero
     u_sig_old(:,:,:) = czero
     nrec = 0
     alpha(:) = czero
     beta(:) = czero

     do iter = 1, maxval(niters)-1
!now need to read residual vectors which have been written to disk in green_linsys:
!          call davcio (alphabeta, lralphabeta, iunalphabeta, iter, -1)
!          alpha = alphabeta(1)
           if(iter.eq.1) then
!initialize variables for shifted system.
!             call davcio (r, lrresid, iunresid, iter, -1)
              do iw = 1, nfreq 
               do ibnd = 1, nbnd
                  alpha(ibnd)      = alphabeta(1,ibnd,iter)
                  u_sig(:,ibnd,iw) = dpsic(:,ibnd,iter)
               enddo
              enddo

              alpha_old(:)      = (1.d0,0.0d0)
              beta_old(:)       = (1.d0,0.0d0)
              pi_coeff_old (:,:)  = (1.d0,0.0d0)
              pi_coeff (:,:)      = (1.d0,0.0d0)
           endif!iter.eq.1

!new pi coefficients.
          do iw = 1, nfreq
            do ibnd = 1 , nbnd
!-alpha because we are solve (H-w^{+}):
               pi_coeff_new(ibnd, iw) = (cone - alpha(ibnd)*DCMPLX(w_ryd(iw), eta))*pi_coeff(ibnd, iw) - &
                                       ((alpha(ibnd)*beta_old(ibnd))/(alpha_old(ibnd)))*(pi_coeff_old(ibnd, iw) - pi_coeff(ibnd, iw))

!beta = (pi_old/pi)**2 *beta, alpha = (pi/pi_new)*alpha
               alpha_sig(ibnd, iw)    = ( pi_coeff(ibnd, iw)/pi_coeff_new(ibnd, iw))*alpha(ibnd)

! x_sig = x_sig + alpha_sig*u_sig
               x_sig(:, ibnd, iw) = x_sig(:,ibnd, iw) + alpha_sig(ibnd, iw) * u_sig(:,ibnd,iw)
            enddo!ibnd
          enddo!iw

!       call davcio (r, lrresid, iunresid, nrec, -1)
          r(:,:) = dpsic(:,:,iter)
!alpha=k, beta = k+1
          beta(ibnd)  = alphabeta(2,nbnd,iter)

!update the u's for the shifted systems
          do iw = 1, nfreq
             do ibnd = 1, nbnd
!u_sig(iw) = (1/\pi_sig)*rp - beta_{sig} * u_{sig} 
                 beta_sig(ibnd,iw)     = ( pi_coeff (ibnd,iw) / pi_coeff_new(ibnd,iw) )**2.d0 * beta(ibnd)
!update the residual and solution vector at each shifted frequency.
                 call ZCOPY (ndim, u_sig (1,nbnd,iw), 1, u_sig_old(1,nbnd,iw), 1)
                 call ZCOPY (ndim, r (1,nbnd), 1, u_sig (1,nbnd,iw), 1)
                 call ZSCAL (ndim, (cone / pi_coeff_new(ibnd, iw)), u_sig(1,ibnd,iw), 1)
                 call ZAXPY (ndim, beta_sig(ibnd,iw), u_sig_old (1,ibnd,iw), 1, u_sig (1,ibnd,iw), 1)
!store values for next iteration
                 alpha_old(ibnd)         = alpha(ibnd)
                 beta_old(ibnd)          = beta(ibnd)
                 pi_coeff_old(ibnd, iw)  = pi_coeff(ibnd, iw)
                 pi_coeff(ibnd, iw)      = pi_coeff_new(ibnd, iw)
             enddo
          enddo!iw
     enddo!iter
  DEALLOCATE(r)
  DEALLOCATE(u_sig)
  DEALLOCATE(u_sig_old) 
END SUBROUTINE coul_multishift

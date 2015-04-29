SUBROUTINE green_multishift_im(ndmx, ndim, nfreq, niters, ngvecs, w_ryd, x_sig)
   USE kinds,       ONLY : DP
   USE units_gw,    ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
   USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax, nwgreen, wgreen
   USE constants,   ONLY : degspin, pi, tpi, RYTOEV, eps8
   USE control_gw,  ONLY : eta, tr2_green
   USE ener,        ONLY : ef

IMPLICIT NONE
!coefficient of quadratic form
  COMPLEX(DP)   :: alpha, beta
!complex(kind=DP), allocatable, INTENT(OUT) :: x_sig (:,:)
  COMPLEX(DP), INTENT(OUT) :: x_sig (ndmx,nfreq)
  COMPLEX(kind=DP), allocatable :: u_sig (:,:), r(:), u_sig_old(:,:), r_sig(:,:)
  COMPLEX(kind=DP) :: alpha_old, beta_old , beta_sig(nfreq), alpha_sig(nfreq)
!pi coefficients for each frequency up to nfreqgreen.
  COMPLEX(kind=DP) :: pi_coeff (nfreq), pi_coeff_old (nfreq), pi_coeff_new(nfreq)
!variable for reading in the stored alpha beta coefficients.
  COMPLEX(DP)                 :: alphabeta(2)
  COMPLEX(DP), PARAMETER      :: cone = (1.0d0,0.0d0), czero=(0.0d0, 0.0d0)
  COMPLEX(DP), external       :: zdotc
  REAL(DP) :: w_ryd(nwgreen)
  REAL(DP) :: mu
  REAL(DP) :: anorm(nwgreen)

  INTEGER  :: nfreq, iw, iwp
  INTEGER  :: ndmx, & ! input: the maximum dimension of the vectors
              ndim, & ! input: the actual dimension of the vectors
              ngvecs,&
              niters,&
              iter,&
              nrec 
  INTEGER  :: ios

!ALLOCATE(x_sig(ndmx,nfreq), r(ndmx))
  ALLOCATE(r(ndmx))
  ALLOCATE(u_sig(ndmx,nfreq), u_sig_old(ndmx,nfreq), r_sig(ndmx,nfreq))
!Green shifted system.
     u_sig(:,:) = czero
     x_sig(:,:) = czero
     r(:)       = czero
     u_sig_old(:,:) = czero
     nrec = 0
     alpha = czero
     beta = czero

     do iter = 1, niters-1
!now need to read residual vectors which have been written to disk in green_linsys:
           call davcio (alphabeta, lralphabeta, iunalphabeta, iter, -1)
           alpha = alphabeta(1)
           if(iter.eq.1) then
!initialize variables for shifted system.
              x_sig(:,:)        = czero
              call davcio (r, lrresid, iunresid, iter, -1)
              do iw = 1, nfreq 
                 u_sig(:,iw) = r(:)
                 r_sig(:,iw) = r(:)
              enddo
              alpha_old         = (1.d0,0.0d0)
              beta_old          = (1.d0,0.0d0)
              pi_coeff_old (:)  = (1.d0,0.0d0)
              pi_coeff (:)      = (1.d0,0.0d0)
           endif!iter.eq.1
!new pi coefficients...
!with rhs as delta <rt,r> is always one on the first iteration.
         do iw = 1, nfreq
!-alpha because we are solve (H-w^{+}):
! conjg means something...
            pi_coeff_new(iw) = (cone - alpha*DCMPLX( 0.0d0 , w_ryd(iw)))*pi_coeff(iw) - &
                              ((alpha*beta_old)/(alpha_old))*(pi_coeff_old(iw) - pi_coeff(iw))
!beta = (pi_old/pi)**2 *beta, alpha = (pi/pi_new)*alpha
            alpha_sig(iw)    = (pi_coeff(iw)/pi_coeff_new(iw))*alpha
! x_sig = x_sig + alpha_sig*u_sig
            x_sig(:,iw)      =  x_sig(:,iw) + alpha_sig(iw) * u_sig(:,iw)
         enddo!iw
!update residuals:
         nrec = iter + 1
         call davcio (r, lrresid, iunresid, nrec, -1)
!alpha=k, beta = k+1
          beta  = alphabeta(2)
!update the u's for the shifted systems
        do iw = 1, nfreq
!u_sig(iw) = (1/\pi_sig)*rp - beta_{sig} * u_{sig} 
          beta_sig(iw)     = ( pi_coeff (iw) / pi_coeff_new(iw) )**2.d0 * beta
!update the residual and solution vector at each shifted frequency.
          call ZCOPY (ndim, u_sig (1,iw), 1, u_sig_old(1,iw), 1)
          call ZCOPY (ndim, r (1), 1, u_sig (1, iw), 1)
          call ZSCAL (ndim, (cone / pi_coeff_new(iw)), u_sig(1,iw), 1)
          call ZAXPY (ndim, beta_sig(iw), u_sig_old (1,iw), 1, u_sig (1,iw), 1)
!store values for next iteration
          alpha_old  = alpha
          beta_old   = beta
          pi_coeff_old (iw) = pi_coeff (iw)
          pi_coeff     (iw) = pi_coeff_new(iw)
        enddo!iw
     enddo!iter

  DEALLOCATE(r)
  DEALLOCATE(u_sig)
  DEALLOCATE(u_sig_old) 
  DEALLOCATE(r_sig)
END SUBROUTINE green_multishift_im

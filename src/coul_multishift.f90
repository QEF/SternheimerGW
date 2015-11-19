  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt.
  !-----------------------------------------------------------------------
SUBROUTINE coul_multishift(ndmx, ndim, nfreq, niters, x_sig, dpsic, alphabeta, freq)
   USE kinds,       ONLY : DP
   USE wvfct,       ONLY : nbnd
   USE units_gw,    ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
   USE freq_gw,     ONLY : fpol, fiu, nwgreen, wgreen
   USE constants,   ONLY : degspin, pi, tpi, RYTOEV, eps8
   USE control_gw,  ONLY : eta, tr2_green,  maxter_coul

IMPLICIT NONE
!coefficient of quadratic form
  INTEGER :: nfreq, iw, iwp
  complex(DP)   :: alpha(nbnd), beta(nbnd)
  complex(kind=DP), allocatable :: u_sig (:,:,:), u_sig_old(:,:,:), r(:,:)
  complex(kind=DP) :: alpha_old(nbnd), beta_old(nbnd) , beta_sig(nbnd, nfreq), alpha_sig(nbnd, nfreq)
!pi coefficients for each frequency up to nfreqgreen.
  complex(kind=DP) :: pi_coeff (nbnd, nfreq), pi_coeff_old (nbnd, nfreq), pi_coeff_new(nbnd, nfreq)
  complex(DP)      :: w_ryd(nfreq)
!  real(DP)     :: h_diag(ndmx, nbnd)

  real(DP) :: anorm(nwgreen)

!variable for reading in the stored alpha beta coefficients.
  complex(DP)                 :: alphabeta(2, nbnd, maxter_coul+1)
  complex(DP)                 :: dpsic(ndmx, nbnd, maxter_coul+1)
  complex(DP)                 :: x_sig (ndmx, nbnd, nfreq)
  complex(DP), PARAMETER      :: cone = (1.0d0,0.0d0), czero=(0.0d0, 0.0d0)
  complex(DP), external       :: zdotc
  complex(DP)                 :: freq(nfreq)
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             ngvecs,&
             iter,&
             nrec,&
             ibnd
  integer :: ios
  integer :: niters(nbnd)

  ALLOCATE(u_sig(ndmx, nbnd, nfreq), u_sig_old(ndmx, nbnd, nfreq))
  w_ryd(:) = freq(:)
!coul shifted system.
  u_sig(:,:,:) = czero
  x_sig(:,:,:) = czero
  u_sig_old(:,:,:) = czero
  nrec = 0
  alpha(:) = czero
  beta(:) = czero

  do iter = 1, maxval(niters)-1
! now need to read residual vectors which have 
! been written to disk in green_linsys:
! call davcio (alphabeta, lralphabeta, iunalphabeta, iter, -1)
! alpha = alphabeta(1)
           do ibnd = 1, nbnd
              if(iter.lt.niters(ibnd)) then
                 alpha(ibnd) = alphabeta(1,ibnd,iter)
              endif
           enddo
           if(iter.eq.1) then
! initialize variables for shifted system.
! call davcio (r, lrresid, iunresid, iter, -1)
              do iw = 1, nfreq 
               do ibnd = 1, nbnd
                  u_sig(:,ibnd,iw) = dpsic(:,ibnd,iter)
               enddo
              enddo
              alpha_old(:)        = (1.d0,0.0d0)
              beta_old(:)         = (1.d0,0.0d0)
              pi_coeff_old (:,:)  = (1.d0,0.0d0)
              pi_coeff (:,:)      = (1.d0,0.0d0)
           endif!iter.eq.1
!new pi coefficients.
           do iw = 1, nfreq
             do ibnd = 1, nbnd
!-alpha because we are solve (H-w^{+}):
                if(iter.lt.niters(ibnd)) then
                   pi_coeff_new(ibnd, iw) = (cone - alpha(ibnd)*w_ryd(iw))*pi_coeff(ibnd, iw) - &
                                            ((alpha(ibnd)*beta_old(ibnd))/(alpha_old(ibnd)))   * &
                                             (pi_coeff_old(ibnd, iw)-pi_coeff(ibnd, iw))
!beta  = (pi_old/pi)**2 *beta, alpha = (pi/pi_new)*alpha
                   alpha_sig(ibnd, iw)    = (pi_coeff(ibnd, iw)/pi_coeff_new(ibnd, iw))*alpha(ibnd)
!x_sig =  x_sig + alpha_sig*u_sig
                   x_sig(:, ibnd, iw)     = x_sig(:,ibnd,iw) + alpha_sig(ibnd,iw) * u_sig(:,ibnd,iw)
                endif
             enddo!ibnd
           enddo!iw
!call davcio (r, lrresid, iunresid, nrec, -1)
!alpha=k, beta = k+1
           do ibnd = 1, nbnd
              if(iter.lt.niters(ibnd)) beta(ibnd) = alphabeta(2,ibnd,iter)
           enddo
!update the u's for the shifted systems
          do iw = 1, nfreq
            do ibnd = 1, nbnd
               if(iter.lt.niters(ibnd)) then
!u_sig(iw) = (1/\pi_sig)*rp - beta_{sig} * u_{sig}
                 beta_sig(ibnd,iw)     = (pi_coeff(ibnd,iw) / pi_coeff_new(ibnd,iw))**2.d0 * beta(ibnd)
!update the residual and solution vector at each shifted frequency.
                 call zcopy (ndim, u_sig(1,ibnd,iw), 1, u_sig_old(1,ibnd,iw), 1)
                 call zcopy (ndim, dpsic (1,ibnd, iter+1), 1, u_sig(1,ibnd,iw), 1)
                 call zscal (ndim, (cone / pi_coeff_new(ibnd, iw)), u_sig(1,ibnd,iw), 1)
                 call zaxpy (ndim, beta_sig(ibnd,iw), u_sig_old (1,ibnd,iw), 1, u_sig (1,ibnd,iw), 1)
!store values for next iteration
                 alpha_old(ibnd)         = alpha(ibnd)
                 beta_old(ibnd)          = beta(ibnd)
                 pi_coeff_old(ibnd, iw)  = pi_coeff(ibnd, iw)
                 pi_coeff(ibnd, iw)      = pi_coeff_new(ibnd, iw)
                endif
              enddo
           enddo!iw
  enddo!iter
  deallocate(u_sig, u_sig_old)
END SUBROUTINE coul_multishift

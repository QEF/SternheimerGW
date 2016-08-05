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
SUBROUTINE green_multishift(ndmx, ndim, nfreq, niters, ngvecs, x_sig)

   USE constants,     ONLY : degspin, pi, tpi, RYTOEV, eps8
   USE control_gw,    ONLY : eta, tr2_green
   USE freq_gw,       ONLY : fiu, nfs, nfsmax, nwgreen, wgreen
   USE kinds,         ONLY : DP
   USE timing_module, ONLY : time_green_multishift
   USE units_gw,      ONLY : iunresid, lrresid, iunalphabeta, lralphabeta

IMPLICIT NONE
!coefficient of quadratic form
  INTEGER :: ndmx ! input: the maximum dimension of the vectors
  INTEGER :: nfreq, iw, iwp
  COMPLEX(DP)   :: alpha, beta
!complex(kind=DP), allocatable, INTENT(OUT) :: x_sig (:,:)
  COMPLEX(DP), INTENT(OUT) :: x_sig (ndmx,nfreq)
  complex(kind=DP), allocatable :: u_sig (:,:), r(:), u_sig_old(:,:), r_sig(:,:)
  complex(kind=DP) :: alpha_old, beta_old , beta_sig(nfreq), alpha_sig(nfreq)
!pi coefficients for each frequency up to nfreqgreen.
  complex(kind=DP) :: pi_coeff (nfreq), pi_coeff_old (nfreq), pi_coeff_new(nfreq)
  real(DP) :: w_ryd(nwgreen)
!HLA should use anorm keep track of divergent elements...
REAL(DP) :: anorm(nwgreen)
!variable for reading in the stored alpha beta coefficients.
  COMPLEX(DP)                 :: alphabeta(2)
  COMPLEX(DP), PARAMETER      :: cone = (1.0d0,0.0d0), czero=(0.0d0, 0.0d0)
  complex(DP), external       :: zdotc
  integer ::   ndim, & ! input: the actual dimension of the vectors
               ngvecs,&
               niters,&
               iter,&
               nrec 
  integer :: ios

  CALL start_clock(time_green_multishift)

!ALLOCATE(x_sig(ndmx,nfreq), r(ndmx))
  ALLOCATE(r(ndmx))
  ALLOCATE(u_sig(ndmx,nfreq), u_sig_old(ndmx,nfreq), r_sig(ndmx,nfreq))
  w_ryd(:) = wgreen(:)/RYTOEV
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
                !u_sig(:,iw) = r(:)
                !r_sig(:,iw) = r(:)
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
            pi_coeff_new(iw) = (cone - alpha*DCMPLX(w_ryd(iw), eta))*pi_coeff(iw) - &
                              ((alpha*beta_old)/(alpha_old))*(pi_coeff_old(iw) - pi_coeff(iw))
! Conjugation so we can take sensible matrix elements
! @
!            pi_coeff_new(iw) = (cone - alpha*DCMPLX(w_ryd(iw), -1.0d0*eta))*pi_coeff(iw) - &
!                              ((alpha*beta_old)/(alpha_old))*(pi_coeff_old(iw) - pi_coeff(iw))
!beta = (pi_old/pi)**2 *beta, alpha = (pi/pi_new)*alpha
            alpha_sig(iw)    = (pi_coeff(iw)/pi_coeff_new(iw))*alpha
! x_sig = x_sig + alpha_sig*u_sig
            x_sig(1:ndim,iw)      =  x_sig(1:ndim,iw) + alpha_sig(iw) * u_sig(1:ndim,iw)
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

  CALL stop_clock(time_green_multishift)

END SUBROUTINE

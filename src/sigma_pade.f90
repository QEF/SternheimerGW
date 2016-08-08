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
subroutine sigma_pade(sigma_band_c, sigma_band_con, wsigwin, nwsigwin)
  use kinds,                only : dp
  use io_global,            only : stdout, ionode_id, ionode
  use io_files,             only : prefix
  use constants,            only : e2, fpi, RYTOEV, tpi, pi
  use gwsigma,              only : ngmsig, nbnd_sig
  use ener,                 only : ef
  use freq_gw,              only : fpol, fiu, nfs, nwsigma, wsigma
  use wvfct,                only : nbnd, npw, npwx, g2kin
  use gvecw,                only : ecutwfc
  use klist,                only : lgauss
  use control_gw,           only : lgamma, eta, godbyneeds, padecont, cohsex, modielec, &
                                   do_diag_g, do_diag_w, trunc_2d, nbnd_occ, double_grid
  implicit none

  integer                  :: nwsigwin
  complex (dp)             :: sigma_band_c(nbnd_sig, nbnd_sig, nwsigma)
  complex (dp)             :: sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin)
  complex(dp), allocatable :: z(:), u(:), a(:)
  complex(dp), allocatable :: z2(:), u2(:), a2(:)
  real(dp)                 :: w_ryd(nwsigma), w_ryd2(nwsigwin), wsigwin(nwsigwin)
  real(dp)                 :: ehomo, elumo, mu
  integer                  :: ig, igp, nw, iw, ibnd, jbnd, ios, &
                              ipol, ik0, ir,irp, counter

!nwsigma is the number of points we have calculated sigma at.
!nwsigwin is the number of points in the sigma window we want
    allocate ( z(nwsigma), a(nwsigma), u(nwsigma))
!Sigma should be calculated on a uniform grid from [0:wsigma]
!We then exploit the symmetry of the selfenergy on the imaginary axis.
!Re(Sig(w)) = Re(Sig(-w)) and Im(Sig(w)) = -Im(Sig(-w))
    allocate ( z2(2*nwsigma-1), a2(2*nwsigma-1), u2(2*nwsigma-1))
    CALL get_homo_lumo (ehomo, elumo)
    if(.not.lgauss) then
      mu = ehomo + 0.5d0*(elumo-ehomo)
    else
      mu = ef
    endif
   !mu = ehomo
   !write(6, '(f14.7)'), mu 
    w_ryd(:)  = wsigma(:)/RYTOEV
    w_ryd2(:) = wsigwin(:)/RYTOEV
if (double_grid) then
    do ibnd =1 , nbnd_sig
        do jbnd = 1, nbnd_sig
            do iw = 1, nwsigma-1
               z2(iw) = dcmplx(mu, -w_ryd(iw+1))
               u2(iw) = conjg(sigma_band_c (ibnd, jbnd, iw+1))
              !z2(iw) = dcmplx(0.0d0, -w_ryd(iw+1))
              !u2(iw) = conjg(sigma_band_c(ibnd, jbnd, iw+1))
              !u2(iw) = sigma_band_c (ibnd, jbnd, iw+1)
            enddo
            do iw = 1, nwsigma 
               z2(iw+nwsigma-1) = dcmplx(mu, w_ryd(iw))
               u2(iw+nwsigma-1) = sigma_band_c (ibnd, jbnd, iw)
               !z2(iw+nwsigma-1) = dcmplx(0.0d0, w_ryd(iw))
               !u2(iw+nwsigma-1) = sigma_band_c (ibnd, jbnd, iw)
            enddo
            call pade_coeff(2*nwsigma-1, z2, u2, a2)
            do iw = 1, nwsigwin
               IF(w_ryd2(iw).lt.mu) THEN
                  call pade_eval(2*nwsigma-1, z2, a2, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               else
                  call pade_eval(2*nwsigma-1, z2, a2, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               endif
            enddo
        enddo
    enddo
else
    do ibnd =1 , nbnd_sig
        do jbnd = 1, nbnd_sig
            do iw = 1, nwsigma
               z(iw) = dcmplx(mu, w_ryd(iw))
               u(iw) = sigma_band_c (ibnd, jbnd, iw)
            enddo
            call pade_coeff(nwsigma, z, u, a)
            do iw = 1, nwsigwin
               if (w_ryd2(iw).lt.mu) then
                  call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               else
                  call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               endif
            enddo
        enddo
    enddo
endif
    deallocate ( z,a,u )
end subroutine sigma_pade

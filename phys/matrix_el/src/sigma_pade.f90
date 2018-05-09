!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
subroutine sigma_pade(sigma_band_c, sigma_band_con, mu, w_ryd, w_ryd2, nwsigwin)

  use control_gw,           only : eta, double_grid
  use freq_gw,              only : nwsigma
  use gwsigma,              only : nbnd_sig
  use kinds,                only : dp

  implicit none

  integer                  :: nwsigwin
  complex (dp)             :: sigma_band_c(nbnd_sig, nbnd_sig, nwsigma)
  complex (dp)             :: sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin)
  complex(dp), allocatable :: z(:), u(:), a(:)
  complex(dp), allocatable :: z2(:), u2(:), a2(:)
  real(dp)                 :: w_ryd(nwsigma), w_ryd2(nwsigwin)
  real(dp)                 :: mu
  integer                  :: iw, ibnd, jbnd

!nwsigma is the number of points we have calculated sigma at.
!nwsigwin is the number of points in the sigma window we want
    allocate ( z(nwsigma), a(nwsigma), u(nwsigma))
!Sigma should be calculated on a uniform grid from [0:wsigma]
!We then exploit the symmetry of the selfenergy on the imaginary axis.
!Re(Sig(w)) = Re(Sig(-w)) and Im(Sig(w)) = -Im(Sig(-w))
    allocate ( z2(2*nwsigma-1), a2(2*nwsigma-1), u2(2*nwsigma-1))
   !mu = ehomo
   !write(6, '(f14.7)'), mu 
if (double_grid) then
    do ibnd =1 , nbnd_sig
        do jbnd = 1, nbnd_sig
            do iw = 1, nwsigma-1
               z2(iw) = CMPLX(mu, -w_ryd(iw+1), KIND=dp)
               u2(iw) = CONJG(sigma_band_c(ibnd, jbnd, iw+1))
            enddo
            do iw = 1, nwsigma 
               z2(iw+nwsigma-1) = CMPLX(mu, w_ryd(iw), KIND=dp)
               u2(iw+nwsigma-1) = sigma_band_c(ibnd, jbnd, iw)
            enddo
            call pade_coeff(2*nwsigma-1, z2, u2, a2)
            do iw = 1, nwsigwin
               IF(w_ryd2(iw).lt.mu) THEN
                  CALL pade_eval(2*nwsigma-1, z2, a2, CMPLX(w_ryd2(iw), eta, KIND=dp), sigma_band_con(ibnd, jbnd, iw))
               else
                  CALL pade_eval(2*nwsigma-1, z2, a2, CMPLX(w_ryd2(iw), eta, KIND=dp), sigma_band_con(ibnd, jbnd, iw))
               endif
            enddo
        enddo
    enddo
else
    do ibnd =1 , nbnd_sig
        do jbnd = 1, nbnd_sig
            do iw = 1, nwsigma
               z(iw) = CMPLX(mu, w_ryd(iw), KIND=dp)
               u(iw) = sigma_band_c (ibnd, jbnd, iw)
            enddo
            call pade_coeff(nwsigma, z, u, a)
            do iw = 1, nwsigwin
               if (w_ryd2(iw).lt.mu) then
                  CALL pade_eval(nwsigma, z, a, CMPLX(w_ryd2(iw), eta, KIND=dp), sigma_band_con(ibnd, jbnd, iw))
               else
                  CALL pade_eval(nwsigma, z, a, CMPLX(w_ryd2(iw), eta, KIND=dp), sigma_band_con(ibnd, jbnd, iw))
               endif
            enddo
        enddo
    enddo
endif
    deallocate ( z,a,u )
end subroutine sigma_pade

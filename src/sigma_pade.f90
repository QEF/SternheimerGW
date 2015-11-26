  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
subroutine sigma_pade(sigma_band_c, sigma_band_con, wsigwin, nwsigwin)
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout, ionode_id, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE gwsigma,              ONLY : ngmsig, nbnd_sig
  USE ener,                 ONLY : ef
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE klist,                ONLY : lgauss
  USE control_gw,           ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec, &
                                   do_diag_g, do_diag_w, trunc_2d, nbnd_occ, double_grid
IMPLICIT NONE

INTEGER                  :: ig, igp, nw, iw, ibnd, jbnd, ios, &
                            ipol, ik0, ir,irp, counter
INTEGER                  :: nwsigwin
COMPLEX (DP)             :: sigma_band_c(nbnd_sig, nbnd_sig, nwsigma)
COMPLEX (DP)             :: sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin)
COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
COMPLEX(DP), ALLOCATABLE :: z2(:), u2(:), a2(:)
REAL(DP)                 :: w_ryd(nwsigwin), w_ryd2(nwsigwin), wsigwin(nwsigwin)
REAL(DP) :: ehomo, elumo, mu

!nwsigma is the number of points we have calculated sigma at.
!nwsigwin is the number of points in the sigma window we want
    ALLOCATE ( z(nwsigma), a(nwsigma), u(nwsigma))
!Sigma should be calculated on a uniform grid from [0:wsigma]
!We then exploit the symmetry of the selfenergy on the imaginary axis.
!Re(Sig(w)) = Re(Sig(-w)) and Im(Sig(w)) = -Im(Sig(-w))
    ALLOCATE ( z2(2*nwsigma-1), a2(2*nwsigma-1), u2(2*nwsigma-1))
    CALL get_homo_lumo (ehomo, elumo)
    if(.not.lgauss) then
      mu = ehomo + 0.5d0*(elumo-ehomo)
    else
      mu = ef
    endif
   !mu = ehomo
   !write(6, '(f14.7)'), mu 
   ! mu = et(nbnd_occ(1), 1) + 0.5d0*(et(nbnd_occ(1)+1, 1) - et(nbnd_occ(1), 1))
    w_ryd(:)  = wsigma(:)/RYTOEV
    w_ryd2(:) = wsigwin(:)/RYTOEV
IF(double_grid) THEN
    DO ibnd =1 , nbnd_sig
        DO jbnd = 1, nbnd_sig
            DO iw = 1, nwsigma-1
               z2(iw) = dcmplx(mu, -w_ryd(iw+1))
               u2(iw) = conjg(sigma_band_c (ibnd, jbnd, iw+1))
              !z2(iw) = dcmplx(0.0d0, -w_ryd(iw+1))
              !u2(iw) = conjg(sigma_band_c(ibnd, jbnd, iw+1))
              !u2(iw) = sigma_band_c (ibnd, jbnd, iw+1)
            ENDDO
            DO iw = 1, nwsigma 
               z2(iw+nwsigma-1) = dcmplx(mu, w_ryd(iw))
               u2(iw+nwsigma-1) = sigma_band_c (ibnd, jbnd, iw)
               !z2(iw+nwsigma-1) = dcmplx(0.0d0, w_ryd(iw))
               !u2(iw+nwsigma-1) = sigma_band_c (ibnd, jbnd, iw)
            ENDDO
            call pade_coeff(2*nwsigma-1, z2, u2, a2)
            DO iw = 1, nwsigwin
               IF(w_ryd2(iw).lt.mu) THEN
                  call pade_eval(2*nwsigma-1, z2, a2, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               ELSE
                  call pade_eval(2*nwsigma-1, z2, a2, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               ENDIF
            ENDDO
        ENDDO
    ENDDO
ELSE
    DO ibnd =1 , nbnd_sig
        DO jbnd = 1, nbnd_sig
            DO iw = 1, nwsigma
               z(iw) = dcmplx(mu, w_ryd(iw))
               u(iw) = sigma_band_c (ibnd, jbnd, iw)
            ENDDO
            call pade_coeff(nwsigma, z, u, a)
            DO iw = 1, nwsigwin
               IF(w_ryd2(iw).lt.mu) THEN
                  call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               ELSE
                  call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               ENDIF
            ENDDO
        ENDDO
    ENDDO
ENDIF
    DEALLOCATE ( z,a,u )
end subroutine sigma_pade

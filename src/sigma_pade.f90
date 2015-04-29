subroutine sigma_pade(sigma_band_c, sigma_band_con, wsigwin, nwsigwin)
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout, ionode_id, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE gwsigma,              ONLY : ngmsig, nbnd_sig
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma
  USE control_gw,           ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec, &
                                   do_diag_g, do_diag_w, trunc_2d
IMPLICIT NONE

INTEGER                  :: ig, igp, nw, iw, ibnd, jbnd, ios, &
                            ipol, ik0, ir,irp, counter
INTEGER                  :: nwsigwin
COMPLEX (DP)             :: sigma_band_c(nbnd_sig, nbnd_sig, nwsigma)
COMPLEX (DP)             :: sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin)
COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
REAL(DP)                 :: w_ryd(nwsigwin), w_ryd2(nwsigwin), wsigwin(nwsigwin)
REAL(DP) :: ehomo, elumo, mu

!nwsigma is the number of points we have calculated sigma at.
!nwsigwin is the number of points in the sigma window we want
    ALLOCATE ( z(nwsigma), a(nwsigma), u(nwsigma))

    CALL get_homo_lumo (ehomo, elumo)
    mu = ehomo + 0.5d0*(elumo-ehomo)

    w_ryd(:)  = wsigma(:)/RYTOEV
    w_ryd2(:) = wsigwin(:)/RYTOEV

    DO ibnd =1 , nbnd_sig
        DO jbnd = 1, nbnd_sig
            DO iw = 1, nwsigma
               z(iw) = dcmplx(mu, w_ryd(iw))
               u(iw) = sigma_band_c (ibnd, jbnd, iw)
            ENDDO
            call pade_coeff(nwsigma, z, u, a)
            DO iw = 1, nwsigwin
               IF(w_ryd2(iw).lt.mu) THEN
                  call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), -eta), sigma_band_con(ibnd, jbnd, iw))
               ELSE
                  call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
               ENDIF
            ENDDO
        ENDDO
    ENDDO
    DEALLOCATE ( z,a,u )
end subroutine sigma_pade


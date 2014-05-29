subroutine sigma_pade(sigma_band_c, sigma_band_con, wsigwin, nwsigwin)
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout, ionode_id, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE gwsigma,              ONLY : ngmsig, nbnd_sig, sigma_g_ex, ngmsco, ngmsex
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

!nwsigma is the number of points we have calculated sigma at.
!nwsigwin is the number of points in the sigma window we wan
    ALLOCATE ( z(nwsigma), a(nwsigma), u(nwsigma))

    w_ryd(:)  = wsigma(:)/RYTOEV
    w_ryd2(:) = wsigwin(:)/RYTOEV

    do ibnd =1 , nbnd_sig
        do jbnd = 1, nbnd_sig
            do iw = 1, nwsigma
               z(iw) = (0.0d0, 1.0d0)*wsigma(iw)
               u(iw) = sigma_band_c (ibnd, jbnd, iw)
            enddo
            call pade_coeff(nwsigma, z, u, a)
            do iw = 1, nwsigwin
           !continuation to the real axis.
               call pade_eval(nwsigma, z, a, dcmplx(w_ryd2(iw), eta), sigma_band_con(ibnd, jbnd, iw))
            enddo
        enddo
    enddo
    DEALLOCATE ( z,a,u )
end subroutine sigma_pade


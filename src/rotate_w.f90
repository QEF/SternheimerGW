  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE rotate_w(scrcoul_g, scrcoul_g_R, gmapsym, eigv, isym, xq_ibk)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, modielec, trunc_2d
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, &
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                            wgreen, wsigma, wsigmamin, wsigmamax, &
                            deltaw, wcoulmax
  USE gwsigma,       ONLY : sigma_c_st
  USE gvect,         ONLY : g, ngm, nl
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE lr_symm_base,  ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q

  IMPLICIT NONE

  COMPLEX(DP)         ::  scrcoul_g   (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)
  COMPLEX(DP)         ::  scrcoul_g_R (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)
  COMPLEX(DP)         :: z(nfs), u(nfs), a(nfs)
  COMPLEX(DP) :: phase

  INTEGER     :: gmapsym  (ngm, nrot) 
  COMPLEX(DP) :: eigv     (ngm, nrot)  

  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: rcut, spal, zcut
  REAL(DP) :: xq_ibk(3), xq_ibz(3)

  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER :: iwim, iw, ikq
  INTEGER :: iqstart, iqstop, iqs, nkr
  INTEGER :: iq, ipol, iqrec, isym

  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limq, inv_q, found

!Rotate G_vectors for FFT.
   DO igp = 1, sigma_c_st%ngmt
      DO ig = 1, sigma_c_st%ngmt
         IF((gmapsym(ig,isym).le.sigma_c_st%ngmt).and.(gmapsym(igp,isym).le.sigma_c_st%ngmt) &
             .and.(gmapsym(ig,isym).gt.0).and.(gmapsym(ig,isym).gt.0)) then
             DO iwim = 1, nfs
               phase = eigv(ig,isym)*conjg(eigv(igp,isym))
               scrcoul_g_R(ig, igp, iwim) = scrcoul_g(gmapsym(ig,isym), gmapsym(igp,isym),iwim)*phase
             ENDDO
         ENDIF
      ENDDO
   ENDDO
   rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
   IF(.not.modielec) THEN
!SPHERICAL SCREENING
    IF(.not.trunc_2d) THEN
       DO iw = 1, nfs
         DO ig = 1, sigma_c_st%ngmt
         qg2 = (g(1,ig) + xq_ibk(1))**2 + (g(2,ig) + xq_ibk(2))**2 + (g(3,ig)+xq_ibk(3))**2
         limq = (qg2.lt.eps8) 
         IF(.not.limq) THEN
            DO igp = 1, sigma_c_st%ngmt
                         !eps^{-1}(\G,\G';\omegaa  )   v(q+G)!
               scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
            ENDDO
         ENDIF
         qg = sqrt(qg2)
         spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
         IF(.not.limq) THEN
            DO igp = 1, sigma_c_st%ngmt
               scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(spal, 0.0d0)
            ENDDO
         ELSE
            IF(iw.eq.1) THEN
               scrcoul_g_R(ig, igp, iw) = real(scrcoul_g_R(ig,igp,iw))
            ENDIF
            DO igp = 1, sigma_c_st%ngmt
               scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
            ENDDO
         ENDIF
        ENDDO!ig
       ENDDO!nfs
       ELSE
            CALL truncate_2d(scrcoul_g_R(1,1,1), xq_ibk, 2)
      ENDIF
    ENDIF
    IF(.NOT.modielec) THEN
        IF(godbyneeds) THEN
        DO ig = 1, sigma_c_st%ngmt
          DO igp = 1, sigma_c_st%ngmt 
!For godby-needs plasmon pole the algebra is done assuming real frequency*i.
!that is: the calculation is done at i*wp but we pass a real number as the freq.
               DO iw = 1, nfs
                  z(iw) = dcmplx(aimag(fiu(iw)), 0.0d0)
                  u(iw) = scrcoul_g_R(ig, igp, iw)
               ENDDO
               CALL godby_needs_coeffs(nfs, z, u, a)
               DO iw = 1, nfs 
!Just overwrite scrcoul_g_R with godby-needs coefficients.
                  scrcoul_g_R (ig, igp, iw) = a(iw)
               ENDDO
          ENDDO
        ENDDO
        ELSE IF (padecont) THEN
          DO igp = 1, sigma_c_st%ngmt
           DO ig = 1, sigma_c_st%ngmt
!Pade input points on the imaginary axis
              DO iw = 1, nfs
                 z(iw) = fiu(iw)
                 u(iw) = scrcoul_g_R (ig, igp, iw)
              ENDDO
              call pade_coeff ( nfs, z, u, a)
!Overwrite scrcoul with Pade coefficients to be passed to pade_eval.
             DO iw = 1, nfs 
                scrcoul_g_R (ig, igp, iw) = a(iw)
             ENDDO
           ENDDO !enddo on ig
        ENDDO  !enddo on igp
        ELSE IF(.not.padecont.and..not.godbyneeds) THEN
                 WRITE(6,'("No screening model chosen!")')
        ENDIF
    ENDIF
END SUBROUTINE rotate_w

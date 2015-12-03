  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE construct_w(scrcoul_g, scrcoul_pade_g, w_ryd)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, modielec, trunc_2d, do_imag
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, &
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                            wgreen, wsigma, wsigmamin, wsigmamax, &
                            deltaw, wcoulmax
  USE gwsigma,       ONLY : sigma_c_st, gcutcorr
  USE gvect,         ONLY : g, ngm, nl
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE mp_global,     ONLY : mp_global_end

  IMPLICIT NONE

  COMPLEX(DP) :: scrcoul_pade_g (gcutcorr, gcutcorr)
  COMPLEX(DP) :: z(nfs), u(nfs), a(nfs)
  COMPLEX(DP)  :: scrcoul_g    (gcutcorr, gcutcorr, nfs) 

  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: w_ryd
  REAL(DP) :: rcut, spal, zcut
  REAL(DP) :: xq_ibk(3), xq_ibz(3)

  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER  :: iwim, iw, ikq

  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limq, inv_q, found

   rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
   scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)
   if(.NOT.modielec) then
     do ig = 1, gcutcorr
        do igp = 1, gcutcorr
           do iwim = 1, nfs
               z(iwim) = fiu(iwim)
               a(iwim) = scrcoul_g (ig,igp,iwim)
           enddo
           if (padecont) then
               call pade_eval ( nfs, z, a, dcmplx(0.0d0, w_ryd), scrcoul_pade_g (ig,igp))
           else if (godbyneeds .and. do_imag) then
               scrcoul_pade_g(ig,igp) = -a(2)/(w_ryd**2+(a(1))**2)
           else if (godbyneeds .and. .not. do_imag) then
               scrcoul_pade_g(ig,igp) = a(2)/(dcmplx(w_ryd**2,0.0d0)-(a(1)-(0.0d0,1.0d0)*eta)**2)
           else
                WRITE(6,'("No screening model chosen!")')
                STOP
                CALL mp_global_end()
           endif
        enddo
     enddo
   else if (modielec) then
       do ig = 1, gcutcorr
          CALL mod_diel(ig, xq_ibk, w_ryd, scrcoul_pade_g(ig,ig), 1)
!scrcoul_pade_g(ig,ig) = mod_dielec_(xq,w_ryd,ig)*v_(q+G,truncation)
           qg2 = (g(1,ig) + xq_ibk(1))**2 + (g(2,ig) + xq_ibk(2))**2 + (g(3,ig)+xq_ibk(3))**2
           limq = (qg2.lt.eps8) 
           if(.not.limq) then
              scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
           endif
           qg = sqrt(qg2)
           spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
!Normal case using truncated coulomb potential.
           if(.not.limq) then
                 scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx(spal, 0.0d0)
           else
                 scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
           endif
        enddo
   endif
end SUBROUTINE construct_w

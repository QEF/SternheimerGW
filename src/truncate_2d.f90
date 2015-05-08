subroutine truncate_2D(scrcoul, xqloc, opt)
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE gvect,         ONLY : g, ngm, nl
  USE qpoint,        ONLY : xqloc, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,       ONLY : sigma_c_st
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag, nspin_gga

IMPLICIT NONE

  LOGICAL  :: limq
  INTEGER  :: ig, igp, iw, opt
  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: xqloc(3)     
  REAL(DP) :: rcut, spal, zcut
  REAL(DP) :: dv_rho(1)
  INTEGER   :: is
  REAL(DP)  :: at1(3,3)
  COMPLEX(DP)  :: scrcoul(sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)
!2D screening.
!Choose zcut to be 1/2*L_{z}.

zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
!from PRB 73, 205119
rcut = -2.0*pi*zcut**2
IF(opt.eq.1) then
           DO ig = 1, ngm
                qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
             IF (qg2 > 1.d-8) then
                qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
                qz   = (g(3,ig)+xqloc(3))
                spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
                dv_rho(nl(ig)) = dv_rho(nl(ig))*dcmplx(spal, 0.0d0)
             ENDIF 
           ENDDO
ELSE IF (opt.eq.2) then
     DO iw = 1, nfs
         DO ig = 1, sigma_c_st%ngmt
                qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
                limq = (qg2.lt.eps8) 
          IF(.not.limq) then
                 qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
                 qz   = sqrt((g(3,ig)+xqloc(3))**2)
                 spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
                 !if(ig.eq.5) print*,scrcoul(1,1,iw), spal, qg2
                 DO igp = 1, sigma_c_st%ngmt
                    scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
                 ENDDO
                 !if(ig.eq.5) print*,scrcoul(1,1,iw), spal, qg2
          ELSE  
                scrcoul(ig, ig, iw) = scrcoul(ig,ig,iw)*dcmplx(rcut,0.0d0)
          ENDIF
         ENDDO
     ENDDO
ENDIF
END subroutine truncate_2D

subroutine truncate_2D(scrcoul, xqloc, opt)
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE gvect,         ONLY : nrxx, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                            nl, ngm, g, nlm
  USE qpoint,        ONLY : xqloc, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
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
  COMPLEX(DP)  :: scrcoul(ngmpol, ngmpol, nfs)
!2D screening.
!Choose zcut to be 1/2*L_{z}.

!rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2))**(float(1)/float(3))

at1(:,:) = at(:,:)
zcut = 0.50d0*sqrt(at1(1,3)**2 + at1(2,3)**2 + at1(3,3)**2)*alat

!from prb 73 205119
!rcut = -2.0*pi*zcut**2
!rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2))**(float(1)/float(3))
rcut = 0.5d0*zcut

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
        DO ig = 1, ngmpol
                qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
                limq = (qg2.lt.eps8) 
          IF(.not.limq) then
           do igp = 1, ngmpol
                 qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
                 qz   = (g(3,ig)+xqloc(3))
                 spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
                 scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
                 scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx(spal, 0.0d0)
           enddo 
          ELSE  
!for omega=0,q-->0, G=0 the real part of the head of the dielectric matrix should be real
!we enforce that here:
           do igp = 1, ngmpol
              scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
           enddo
          ENDIF
        ENDDO
      ENDDO
ENDIF
END subroutine truncate_2D

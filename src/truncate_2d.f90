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
subroutine truncate_2D(scrcoul, xqloc, opt)
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE gvect,         ONLY : g, ngm, nl
  USE qpoint,        ONLY : npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,       ONLY : sigma_c_st
  USE fft_base,      ONLY : dfftp
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE control_gw,  ONLY : multishift
IMPLICIT NONE
  LOGICAL  :: limq
  INTEGER  :: ig, igp, iw, opt
  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: xqloc(3)     
  REAL(DP) :: rcut, spal, zcut
  REAL(DP) :: dv_rho(dfftp%nnr)
  INTEGER  :: is
  REAL(DP) :: at1(3,3)
  COMPLEX(DP)  :: scrcoul(sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)
!2D screening.
!Choose zcut to be 1/2*L_{z}.
!zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
!from PRB 73, 205119
rcut = -2*pi*zcut**2

IF(opt.eq.1) then
     DO ig = 1, ngm
           qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
        IF (qg2 > 1.d-8) then
           qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
           qz   = sqrt((g(3,ig)+xqloc(3))**2)
           spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
           dv_rho(nl(ig)) = dv_rho(nl(ig)) + dcmplx(fpi*e2/(tpiba2*qg2)*spal, 0.0d0)
        ENDIF 
     ENDDO
!    DO ig = 1, ngm
!       qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
!       qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
!       qz   = sqrt((g(3,ig) + xqloc(3))**2)
!       IF(qxy.gt.eps8) then
!          spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
!          DO igp = 1, sigma_c_st%ngmt
!             dvaux(nl(ig)) = dvaux(nl(ig)) + dvscf(nl(ig),1)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
!          ENDDO
!       ELSE IF(qxy.lt.eps8.and.qz.gt.eps8) then
!          spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
!          DO igp = 1, sigma_c_st%ngmt
!             dvaux(nl(ig)) = dvaux(nl(ig)) + dvscf(nl(ig),1)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
!          ENDDO
!       ELSE
!          dvaux(nl(ig)) = 0.0d0
!       ENDIF
!    ENDDO
ELSE IF (opt.eq.2) then
     DO iw = 1, nfs
         DO ig = 1, sigma_c_st%ngmt
                 qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
                 qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
                 qz   = sqrt((g(3,ig)+xqloc(3))**2)
                 limq = (qg2.lt.eps8) 
                 if(qg2.gt.eps8) then
                    spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*(cos(tpiba*qz*zcut))
                    do igp = 1, sigma_c_st%ngmt
                       scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
                    enddo
                 else  
                    do igp=1,sigma_c_st%ngmt
                       scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx(rcut,0.0d0)
                    enddo
                 endif
!          IF(qxy.gt.eps8) then
!                 spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
!                 DO igp = 1, sigma_c_st%ngmt
!                    scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
!                 ENDDO
!          ELSE IF(qxy.lt.eps8.and.qz.gt.eps8) then
!                 spal = 1.0d0 - cos(tpiba*qz*zcut) - tpiba*qz*zcut*sin(tpiba*qz*zcut)
!                 DO igp = 1, sigma_c_st%ngmt
!                    scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
!                 ENDDO
!          ELSE  
!                 do igp=1,sigma_c_st%ngmt
!                    scrcoul(ig, igp, iw) = scrcoul(ig,igp,iw)*dcmplx(rcut,0.0d0)
!                 enddo
!          ENDIF
         ENDDO
     ENDDO
ENDIF
END subroutine truncate_2D

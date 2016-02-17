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
SUBROUTINE spheric_coulomb(diel_matrix, xq_coul)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE cell_base,  ONLY : alat, tpiba2, omega
  USE gvect,      ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE disp,        ONLY : nqs, nq1, nq2, nq3

  IMPLICIT NONE

!SUBROUTINE SPHERIC CUT  COULOMB POTENTIAL Appliead to v_q.
!Spencer/Alavi truncation of the bare coulomb interaction
![PRB 77,193110 (2008)]

  LOGICAL   :: limq
  REAL(DP)  :: xq_coul(3)
  REAL(DP)  :: qg, rcut, spal
  INTEGER   :: ig, igp, iw
  INTEGER     :: unf_recl, recl, ios
 !COMPLEX(DP), INTENT(INOUT) :: diel_matrix(ngmpol,ngmpol,nfs) 
  COMPLEX(DP) :: diel_matrix(ngmpol,ngmpol,nfs) 
  REAL(DP)  :: qg2, qg2coul


rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))

DO iw = 1, nfs
   DO ig = 1, ngmpol
       qg2 = (g(1,ig) + xq_coul(1))**2 + (g(2,ig) + xq_coul(2))**2 + (g(3,ig) + xq_coul(3))**2
       !if(qg2.lt.eps8) limq =.true.
       limq = (qg2.lt.eps8) 
       IF(.not.limq) then
           DO igp = 1, ngmpol
              diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
           ENDDO
       ELSE 
           if(ig.ne.1) then
              DO igp = 1, ngmpol
                 diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
              ENDDO
           endif
       ENDIF

       qg = sqrt(qg2)
       spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)

!Normal case using truncated coulomb potential.
       if(.not.limq) then
          do igp = 1, ngmpol
              diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx(spal, 0.0d0)
          enddo
       else
!should only occur case iq->0, ig = 0 use vcut (q(0) = (4pi*e2*Rcut^{2})/2
             write(6,'("Taking Limit.")')
             write(6,*) (fpi*e2*(rcut**2))/2.0d0
             write(6,*) ig, iw
             write(6,*) g(:, ig)
             !for omega=0,q-->0, G=0 the real part of the head of the dielectric matrix should be real
             !we enforce that here:
         if(iw.eq.1) then
            diel_matrix(ig, igp, iw) = real(diel_matrix(ig,igp,iw))
         endif
         do igp = 1, ngmpol
            diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
         enddo
       endif
   ENDDO !ig
ENDDO!iw
END SUBROUTINE

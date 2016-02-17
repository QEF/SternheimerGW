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
SUBROUTINE mod_diel(ig, xq_ibk, w, inveps, screening)
  USE kinds,      ONLY : DP
  USE constants,  ONLY : e2, fpi, RYTOEV, pi, eps8
  USE gvect,      ONLY : ngm, g
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at

IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: xq_ibk(3)
  COMPLEX(DP) ::  inveps
  REAL(DP)    :: w
  REAL(DP)    :: wwp, eps0, q0, wwq, fac, z, xq(3)
  REAL(DP)    :: wwpi2, wwpj2, qxy, meff, invepsqg
  REAL(DP)    :: rhoi, rhoj , wwqi, wwji
  REAL(DP)    :: wwqi2, wwqj2, wwji2, wwq2
  REAL(DP)    :: A, B, wwmi2, wwpl2 
  REAL(DP)    :: epsq, polq
  REAL(DP)    :: qg, qg2, kf2, kf, alpha
  REAL(DP)    :: expqxy, d, x 
  INTEGER     :: screening, ig
!
!Screening selects from one of the following model dielectric functions.
!in all cases the dielectric function is diagonal in G space so we include no
!local field effects. However the model is parameterized to include the possibility
!of two modes occurring from coupling between two electron gases of different density.
!
!The 'magic dielectric function' Inkson 1972
!Silicon parameters...
!check resta for parameters Phys. Rev. B 16, 2717 2722 (1977)
!wwp    = 16.0/RYTOEV  ! plasma frequency in Ry
!eps0   = 11.4         ! static diel constant of Si
!q0     = 1.1          ! characteristic momentum of Si, a.u. from Resta
!LiCL parameters
!wwp    = 17.0/RYTOEV   ! plasma frequency in Ry
!eps0   = 11.04         ! static diel constant of 
!q0     = 1.2           ! characteristic momentum of Si, a.u. from Resta
!MoS2 there are two well defined excitation for parallel and perpendicular.
!this is an anisotropic material though!
!should have (at least) 2 different plasmons!
     wwp    = 12.0_dp/RYTOEV ! plasma frequency in Ry
     eps0   =  7.4_dp         ! static diel constant of MoS2
     q0     =  1.9_dp        ! characteristic momentum of MoS2 calculated from Resta paper..
!!!!!!
     fac    = 1.0_dp/(1.0_dp - 1.0_dp/eps0)
!effective mass annd density of electrons in layer i,j.
     meff   = 1.0_dp
!or should this be the fourier co-efficient of the density at the particular G vector?
     rhoi   = 1.0_dp
     rhoj   = 1.0_dp
!z is interlayer distance.
     z = 8.5d0
 if (screening.eq.1) then
!Standard 3D-bulk ppm:
!drhoscfs (nl(ig), 1)  = 1.d0 - wwp**2.d0/((fiu(iw) + eta)**2.d0 + wwq**2.d0)
!(W-v) = (inveps(w) - delta) v
!Wave vector:
     qg2    = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
     qg     = sqrt(tpiba2*qg2)
     fac    = 1.0_dp/(1.0_dp-1.0_dp/eps0)
     wwq    = wwp * sqrt ( fac * (1.0_dp + (qg/eps0/q0)**2.0_dp ) )
     !inveps =  - wwp**2.d0/(w**2.d0 + wwq**2.d0)

 else if (screening.eq.2) then
 !Effective bi-layer plasmon model
 !Inkson and White semicond. sci. technol. 4 1989
 !only xy screened and effective coupling along z...
     qg2     = (g(1,ig)+xq(1))*2 + (g(2,ig)+xq(2))**2
     qg      = sqrt(tpiba2*qg2)

     wwpi2   = (2.0_dp*pi)/(eps0*meff)*qg
     wwpj2   = (2.0_dp*pi)/(eps0*meff)*qg

     qxy     = sqrt(tpiba2*((xq(1)+g(1, ig))**2 + (xq(2)+g(2, ig))**2))
     alpha   = wwpi2 + wwpj2*exp(-qxy*z)

!Sterne polarizability.
!    polq    =  ((eps0*meff)/(kf*x) - sqrt((epsq * meff/(kf*qg))-1))
     epsq    = qg2
     polq    = (-2*kf*meff)/(pi*qg)*((epsq*meff)/(kf*x) - real(sqrt((epsq * meff/(kf*qg))-1)))
!vq is screened coulomb, what is epsq?? need to look at dimensions.
!looks like an energy actually.
     !b       = cosh(qxy*d) - vq*polq*sinh(qxy*z) 
!eps^{-1} = sinh(qd)/\sqrt(b^{2}-1)
     inveps  = sinh(qg*z)/sqrt(b**2 -1)
!need to use inkson's relation here... wwq = alpha*(1-inveps(q,0))^{-1}
     wwq2    = alpha/(1 - inveps)
!Skip through frequencies.
     inveps  = alpha/(w**2 - wwq2)
 else if (screening.eq.3) then
!Bilayer plasmon model 
     qg2     = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2
     qg      = sqrt(tpiba2*qg2)
     wwpi2   = ((2*pi*rhoi)/(eps0*meff))*qg
     wwpj2   = ((2*pi*rhoj)/(eps0*meff))*qg
!Still need to put in wwq modes.
!     wwqi2   =
!     wwqj2   =
     expqxy  = exp(-2*qg*z)
!correct expression for two poles:
     wwpl2   = ((wwqi2+wwqj2)/2) + sqrt(((wwqi2-wwqj2)/2)**2 + wwpi2*wwpj2*expqxy)
     wwmi2   = ((wwqi2+wwqj2)/2) - sqrt(((wwqi2-wwqj2)/2)**2 + wwpi2*wwpj2*expqxy)
!some combination to ensure appropriate simplification:
!     A       = 
!     B       = 1 - A
     inveps  =  A/(w**2 - wwpl2) + B/(w**2- wwmi2)
!else if(screening.eq.4) then
!Multilayer plasmon model gives the expression for stern 2-D dielectric response.
!analytic layered electron gas inverse dielectric function from hawyrlak and co-workers.
!bqg=cosh(qg*d)-((2*pi)/(eps0*x))*((eps0*meff)/(kf*x)-sqrt((eps0*meff/(kf*x))-1))*sinh(qxy*d)
!# \inveps(q,0) = sinh(qd)/sqrt(b**2-1)
! fqg  = sinh(x*d)/(sqrt(b(x)**2-1))
!#alpha parameter
! aqg     = ((2*pi*n0)/(eps0*meff)*x)/(1-exp(-2*x*d))
! weff2qg = aqg/(1-fqg)
!#finally the static dielectric fxn is
!invepsqg = 1 - aqg/weff2qg
! else
!     STOP
 endif
END SUBROUTINE mod_diel

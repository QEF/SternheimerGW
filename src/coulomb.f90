! Copyright (C) 2001-2008 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
SUBROUTINE coulomb(iq, igstart, igstop, scrcoul) 
!-----------------------------------------------------------------------
! This subroutine is the main driver of the COULOMB self consistent cycle
! which calculates the Screened Coulomb interaction along the imaginary axis
! a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G). 

!What we actually calculate is the dielectric matrix.
!This can then be stored and we applye the coulom operator and do all the analytic continuation 
!stuff in sigma_c. This means we can actually play around with the dielectric matrix 
!as a function of frequency which I think make's more sense anyway to store.

!The dielectric matrix is given by:
!eps_{q}^{-1}(G,G',iw) - (\delta_{GG'} = drhoscfs_{G,G',iw})


  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE gvect,      ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,    ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE constants,  ONLY : e2, fpi, RYTOEV, pi, eps8
  USE cell_base,  ONLY : alat, tpiba2, omega
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout, ionode
  USE uspp,       ONLY: okvan
  USE control_gw, ONLY : zue, convt, rec_code, modielec, eta, godbyneeds, padecont
  USE partial,    ONLY : done_irr, comp_irr
  USE modes,      ONLY : nirr, npert, npertx
  USE uspp_param, ONLY : nhm
  USE eqv,        ONLY : drhoscfs, dvbare
  USE paw_variables,    ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE recover_mod, ONLY : write_rec
  USE gwsigma,     ONLY : ngmpol
  USE qpoint,      ONLY : xq
  USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
  USE units_gw,    ONLY : iuncoul, lrcoul
  USE disp,        ONLY : nqs, nq1, nq2, nq3

! PARALLEL STUFF
  USE mp_global,   ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime
  USE mp,          ONLY : mp_barrier, mp_bcast, mp_sum
 
 !Symmetry Stuff
  USE gwsymm,    ONLY : ig_unique, ngmunique

  IMPLICIT NONE
  REAL(DP) :: tcpu, get_clock

! timing variables
  REAL(DP) :: qg2, qg2coul
  INTEGER :: ig, igp, iw, npe, irr, icounter
  INTEGER :: igstart, igstop, igpert
  COMPLEX(DP), allocatable :: drhoaux (:,:) 
  COMPLEX(DP) :: padapp, w
!HL temp variable for scrcoul to write to file.  
  COMPLEX(DP) :: cw
  COMPLEX(DP), allocatable :: z(:), u(:), a(:)
  INTEGER :: unf_recl, recl, ios
  INTEGER :: iq, screening 
  LOGICAL :: exst
!again should decide if this should be allocated globally. 
  COMPLEX(DP) :: scrcoul(ngmpol, ngmpol, nfs, 1)
!modeps and spencer-alavi vars
  REAL(DP) :: wwp, eps0, q0, wwq, fac
  REAL(DP) :: qg, rcut, spal
!for Godby needs plasmon pole.
  LOGICAL :: diag, limit
! used to test the recover file
  EXTERNAL get_clock
!Extended plasmon pole model
  REAL(DP) :: wwpi2, wwpj2, qxy, meff
  CALL start_clock ('coulomb')

!DUMMY VARIABLES
!Change in charge density

ALLOCATE (drhoscfs(nrxx , nspin_mag))    
ALLOCATE ( z(nfs), u(nfs), a(nfs) )

irr=1
scrcoul(:,:,:,:) = (0.d0, 0.0d0)
!LOOP OVER ig, unique g vectors only... these then get written into the full matrix.
!g is sorted in magnitude order...
DO ig = igstart, igstop
   qg2 = (g(1,ig_unique(ig))+xq(1))**2 + (g(2,ig_unique(ig))+xq(2))**2 + (g(3,ig_unique(ig))+xq(3))**2
    do iw = 1, nfs
!HL DEBUG BROYDEN
       !do iw = nfs, nfs
         drhoscfs(:,:) = (0.0d0, 0.0d0)
         dvbare(:)     = (0.0d0, 0.0d0)
         dvbare (nl (ig_unique(ig)) ) = (1.d0, 0.d0)
       !From the inputcard we can choose whether to use a model dielectric or 
       !do the full sternheimer treatment.
       IF (modielec) then
!      !The 'magic dielectric function' Inkson 1972
!      !Silicon parameters...
!      !check resta for parameters Phys. Rev. B 16, 2717 2722 (1977)
!      !wwp    = 18.0/RYTOEV  ! plasma frequency in Ry
!      !eps0   = 11.4         ! static diel constant of Si
!      !q0     = 1.1          ! characteristic momentum of Si, a.u. from Resta
!      !LiCL parameters
!      !wwp    = 17.0/RYTOEV   ! plasma frequency in Ry
!      !eps0   = 11.04         ! static diel constant of 
!      !q0     = 1.2           ! characteristic momentum of Si, a.u. from Resta
!      !MoS2 there are two well defined excitation for parallel and perpendicular.
!      !this is an anisotropic material though!
!      !should have (at least) 2 different plasmons!
        wwp    = 24.0/RYTOEV ! plasma frequency in Ry
        eps0   = 7.4         ! static diel constant of MoS2
        q0     = 1.90        ! characteristic momentum of MoS2 calculated from Resta paper..
        qg     = sqrt(tpiba2*qg2)
        fac    = 1.d0/(1.d0-1.d0/eps0)
        wwq    = wwp * sqrt ( fac * (1.d0 + (qg/eps0/q0)**2.d0 ) )
!        meff   = 1.0d0
!        qg      = sqrt(tpiba2*qg2)
!        if(screening.eq.1) then
!     !Standard 3D-bulk ppm:
!     !diagonal term ig = igp (all the others remain 0)
  !   drhoscfs (nl(ig), 1)  = 1.d0 - wwp**2.d0/((fiu(iw) + eta)**2.d0 + wwq**2.d0)
!     !(W-v) = (inveps(w) - delta) v
!           fac    = 1.d0/(1.d0-1.d0/eps0)
!           wwq    = wwp * sqrt ( fac * (1.d0 + (qg/eps0/q0)**2.d0 ) )
           drhoscfs (nl(ig_unique(ig)), 1)  = - wwp**2.d0/((fiu(iw) + eta)**2.d0 + wwq**2.d0)
!        else if(screening.eq.2) then
!     !Effective bi-layer plasmon model
!     !Inkson and White semicond. sci. technol. 4 1989
!           wwpi2   = (2*pi)/(eps0*meff)*qg
!           wwpj2   = (2*pi)/(eps0*meff)*qg
!           qxy     = sqrt(tpiba2((q(1)+g(1, ig_unique(ig)))**2 + (q(2)+g(2, ig_unique(ig)))**2))
!           alpha   = wwpi2 + wwpj2*exp(-qxy*z)
!     !need to use inkson's relation here... wwq = alpha*(1-inveps(q,0))^{-1}
!           wwq  = alpha / (1 - inveps(q,0))
!     !Skip through frequencies.
!           scrcoul(ig_unique(ig), igp, 1, nspin_mag) = alpha
!           scrcoul(ig_unique(ig), igp, 2, nspin_mag) = wwq
!           CYCLE
!        else if(screening.eq.3) then
     !Bilayer plasmon model 
        !   wwpi2   = (2*pi)/(eps0*meff)*qg
        !   wwpj2   = (2*pi)/(eps0*meff)*qg

        !   wwpl = ()
        !   wwmi = ()
        !else if(screening.eq.4) then
!Multilayer plasmon modelM gives the expression for stern 2-D dielectric response.
!analytic layered electron gas inverse dielectric function from hawyrlak and co-workers.
!bqg=cosh(qg*d)-((2*pi)/(eps0*x))*((eps0*meff)/(kf*x)-sqrt((eps0*meff/(kf*x))-1))*sinh(qxy*d)
!# \inveps(q,0) = sinh(qd)/sqrt(b**2-1)
! fqg  = sinh(x*d)/(sqrt(b(x)**2-1))
!#alpha parameter
! aqg     = ((2*pi*n0)/(eps0*meff)*x)/(1-exp(-2*x*d))
! weff2qg = aqg/(1-fqg)
!#finally the static dielectric fxn is
!invepsqg = 1 - aqg/weff2qg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !endif
       ELSE
         call cft3 (dvbare, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)
         CALL solve_linter (dvbare, iw, drhoscfs)
         call cft3 (drhoscfs, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
         call cft3 (dvbare, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
         if(iq.eq.1) then
            WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscfs(nl(ig_unique(ig)), 1) + dvbare(nl(ig_unique(ig)))
         endif
         DO igp = 1, ngmpol
            scrcoul(ig_unique(ig), igp, iw, nspin_mag) = drhoscfs(nl(igp),1)
         ENDDO
       ENDIF
    enddo !iw
ENDDO 
tcpu = get_clock ('GW')
DEALLOCATE (drhoscfs)
DEALLOCATE ( z, u, a )
CALL stop_clock ('coulomb')
RETURN
END SUBROUTINE coulomb

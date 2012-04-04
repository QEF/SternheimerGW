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
! and then does the pade analytic continuation on to the real axis.
! a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G). 
! The screened coulomb interaction is given by 
! W_{q}(G,G',iw) = (\delta_{GG'} + drhoscfs_{G,G',iw}) * (e2*fpi)/(q2g*tpiba2)
! W = eps^{-1}(iw) v .

  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE gvect,      ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,    ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE constants,  ONLY : e2, fpi, RYTOEV, pi, eps8
  USE cell_base,  ONLY : alat, tpiba2, omega
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout, ionode
  USE uspp,       ONLY: okvan
  USE control_gw, ONLY : zue, convt, rec_code, modielec, eta
  USE partial,    ONLY : done_irr, comp_irr
  USE modes,      ONLY : nirr, npert, npertx
  USE uspp_param, ONLY : nhm
  USE eqv,        ONLY : drhoscfs, dvbare
  USE paw_variables,    ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE recover_mod, ONLY : write_rec
  USE gwsigma,     ONLY : ngmsig
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

  REAL(DP) :: qg2
  INTEGER :: ig, igp, iw, npe, irr, icounter
  INTEGER :: igstart, igstop, igpert
 
  COMPLEX(DP), allocatable :: drhoaux (:,:) 
  COMPLEX(DP) :: padapp, w
  !HL temp variable for scrcoul to write to file.  
  COMPLEX(DP) :: cw
  COMPLEX(DP), allocatable :: z(:), u(:), a(:)
  INTEGER :: unf_recl, recl, ios
  INTEGER :: iq 
  LOGICAL :: exst
 !again should decide if this should be allocated globally. 
  COMPLEX(DP) :: scrcoul(ngmsig, ngmsig, nfs, 1)

 !COMPLEX(DP) :: scrcoul(ngmsig, ngmsig, nwcoul,1)
 !modeps and spenceralavi vars
  REAL(DP) :: wwp, eps0, q0, wwq, fac
  REAL(DP) :: qg, rcut, spal
!  REAL(DP) :: w_ryd(nwcoul)

! used to test the recover file
  EXTERNAL get_clock
  CALL start_clock ('coulomb')

!DUMMY VARIABLES
!Change in charge density

ALLOCATE (drhoscfs(nrxx , nspin_mag))    
ALLOCATE ( z(nfs), u(nfs), a(nfs) )

irr=1

scrcoul(:,:,:,:) = (0.d0, 0.0d0)

!LOOP OVER ig, unique g vectors only... these then get written into the full matrix.
DO ig = igstart, igstop
    qg2 = (g(1,ig_unique(ig))+xq(1))**2 + (g(2,ig_unique(ig))+xq(2))**2 + (g(3,ig_unique(ig))+xq(3))**2
!   if (qg2.gt.eps8) then
      do iw = 1, nfs
         drhoscfs(:,:) = (0.0d0, 0.0d0)
         dvbare(:)     = (0.0d0, 0.0d0)
         dvbare (nl (ig_unique(ig)) ) = (1.d0, 0.d0)

!    WRITE(1000+mpime,'(4x,"Frequency= ",3f7.3)') fiu(iw)*RYTOEV
!    WRITE(1000+mpime,'(4x,"Screened Coulomb: q =",3f7.3,"  G =",3f7.3)') xq(:), g(:,ig)

      !From the inputcard we can choose whether to use a model dielectric or 
      !do the full sternheimer treatment.

       IF (modielec) then
        !The 'magic dielectric function' Inkson 1972
         wwp    = 18.0/RYTOEV  ! plasma frequency in Ry
         eps0   = 11.4          ! static diel constant of Si
         q0     = 1.1           ! characteristic momentum of Si, a.u. from Resta
         qg     = sqrt(tpiba2*qg2)
         fac    = 1.d0/(1.d0-1.d0/eps0)
         wwq    = wwp * sqrt ( fac * (1.d0 + (qg/eps0/q0)**2.d0 ) )

        !diagonal term ig = igp (all the others remain 0)
        !drhoscfs (nl(ig), 1)  = 1.d0 - wwp**2.d0/((fiu(iw) + eta)**2.d0 + wwq**2.d0)
        !(W-v) = (inveps(w) - delta) v
         drhoscfs (nl(ig_unique(ig)), 1)  = - wwp**2.d0/((fiu(iw) + eta)**2.d0 + wwq**2.d0)

!         WRITE(1000+mpime, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscfs(nl(ig),1) + dvbare(nl(ig))

       ELSE

         call cft3 (dvbare, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)
         CALL solve_linter (dvbare, iw, drhoscfs)

        ! IF (convt .ne. .true.) WRITE( stdout, '(/,5x,"No convergence ")')
        ! IF (convt .ne. .true.) WRITE( 1000+mpime, '(/,5x,"No convergence ")')

         call cft3 (drhoscfs, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
         call cft3 (dvbare, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
        ! WRITE(1000+mpime, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscfs(nl(ig),1) + dvbare(nl(ig))
         if(iq.eq.1) then
         WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscfs(nl(ig_unique(ig)), 1) + dvbare(nl(ig_unique(ig)))
         endif
       ENDIF

!Generate an entire row of the irreducible Screened Coulomb Matrix with 4pi*e2*|q+G|^{-2}eps_{GG'}. 
!W(Girr, G'irr) (irr= irreducible, red=reducible)
!Then W(Gred, G'red) = \sum_{Girr, Girr')\sum_{R \in Gq} W(Girr, G'irr)

       DO igp = 1, ngmunique
          scrcoul(ig_unique(ig),ig_unique(igp),iw,nspin_mag) = drhoscfs(nl(ig_unique(igp)),1)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
       ENDDO

! Spencer/Alavi truncation of the bare coulomb interaction
! [PRB 77,193110 (2008]

        rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
        qg = sqrt( (g(1,ig_unique(ig) )+xq(1))**2.d0 + (g(2,ig_unique(ig) )+xq(2))**2.d0 + (g(3,ig_unique(ig) )+ xq(3))**2.d0 )
        spal = 1.0d0 - cos ( rcut * sqrt(tpiba2) * qg )

        !do igp = 1, ngmsig
        do igp = 1, ngmunique
           scrcoul (ig_unique(ig), ig_unique(igp), iw, nspin_mag) = scrcoul (ig_unique(ig), ig_unique(igp),iw,nspin_mag) * dcmplx ( spal, 0.0d0)
        enddo
      enddo !iw

! Analytical Continuation to the real axis.
! @10TION only diagonal elements when using the model-dielectric function.
    IF (modielec) then
       do igp = ig, ig
!   Pade input points on the imaginary axis
         do iw = 1, nfs
           z(iw) = dcmplx( 0.d0, fiu(iw))
!          u(iw) = scrcoul (ig,igp,iw,nspin_mag)
           u(iw) = scrcoul (ig_unique(ig), ig_unique(igp),iw,nspin_mag)
         enddo
!        Pade coefficients
         call pade_coeff ( nfs, z, u, a)
!        Overwrite scrcoul with Pade coefficients to be read in gw_product.
!        write(6,'("Pade Co-efficients")')
         do iw = 1, nfs
           !scrcoul (ig,igp,iw,nspin_mag) = a(iw)
            scrcoul (ig_unique(ig),ig_unique(igp),iw,nspin_mag) = a(iw)
         enddo
       enddo !enddo on igp
    ELSE
       do igp = 1, ngmunique
!   Pade input points on the imaginary axis
         do iw = 1, nfs
           z(iw) = dcmplx( 0.d0, fiu(iw))
           u(iw) = scrcoul (ig_unique(ig), ig_unique(igp),iw,nspin_mag)
         enddo
!        Pade coefficients      
         call pade_coeff ( nfs, z, u, a)
!        Overwrite scrcoul with Pade coefficients to be read in gw_product.
!        write(6,'("Pade Co-efficients")')
         do iw = 1, nfs 
           scrcoul (ig_unique(ig),ig_unique(igp),iw,nspin_mag) = a(iw)
         enddo
       enddo !enddo on igp
    ENDIF
ENDDO 
!on G
!Here writing to disc the pade co-efficients of the analytically continued screened
!coulomb interaction. These co-efficients are then read back in when the GW product
!routine is called.  
!  CALL davcio(scrcoul, lrcoul, iuncoul, iq, 1 )
!  CALL mp_barrier()

tcpu = get_clock ('GW')

DEALLOCATE (drhoscfs)
DEALLOCATE ( z, u, a )

      CALL stop_clock ('coulomb')
RETURN
END SUBROUTINE coulomb

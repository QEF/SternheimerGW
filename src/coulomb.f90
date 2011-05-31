! Copyright (C) 2001-2008 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------

SUBROUTINE coulomb(iq) 

!-----------------------------------------------------------------------
! This subroutine is the main driver of the COULOMB self consistent cycle
! which calculates the Screened Coulomb interaction along the imaginary axis
! and then does the pade analytic continuation on to the real axis.
! a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G). 
! The screened coulomb interaction is given by 
! W_{q}(G,G',iw) = (\delta_{GG'} + drhoscfs_{G,G',iw}) * (e2*fpi)/(q2g*tpiba2)
! W = eps^{-1} v 

  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat
  USE gvect,    ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth, ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE constants, ONLY : e2, fpi, RYTOEV, pi, eps8
  USE cell_base, ONLY : alat, tpiba2, omega
  USE lsda_mod, ONLY : nspin
  USE io_global,  ONLY : stdout, ionode
  USE uspp,  ONLY: okvan
  USE control_gw, ONLY : zue, convt, rec_code, modielec
  USE partial,    ONLY : done_irr, comp_irr
  USE modes,      ONLY : nirr, npert, npertx
  USE uspp_param, ONLY : nhm
  USE eqv,        ONLY : drhoscfs, dvbare
  USE paw_variables, ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE recover_mod, ONLY : write_rec
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE gwsigma,     ONLY : ngmsig
  USE qpoint,      ONLY : xq
  USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
  USE units_gw,    ONLY : iuncoul, lrcoul
  USE disp,        ONLY : nqs

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
  REAL(DP), allocatable :: scrcoul_g(:) 
  COMPLEX(DP), allocatable :: z(:), u(:), a(:)
  INTEGER :: unf_recl, recl, ios
  INTEGER :: iq 
  REAL(DP), parameter :: eta = 0.04
  LOGICAL :: exst
  !again should decide if this should be allocated globally. 
  !COMPLEX(DP) :: scrcoul(ngmsig, ngmsig, nfs,1)
  COMPLEX(DP) :: scrcoul(ngmsig, ngmsig, nwcoul,1)
  !modeps and spenceralavi vars
  REAL(DP) :: wwp, eps0, q0, wwq, fac
  REAL(DP) :: qg, rcut, spal
  REAL(DP) :: w_ryd(nwcoul)

! used to test the recover file
  EXTERNAL get_clock
  CALL start_clock ('coulomb')

!DUMMY VARIABLES
!Change in charge density

ALLOCATE (drhoscfs(nrxx , nspin_mag))    
ALLOCATE (scrcoul_g(20)) 
ALLOCATE ( z(nfs), u(nfs), a(nfs) )

!Defining iuncoul variable here temporarily for testing purposes. Need to put it in the 'right'
!opening file spot as usually done in quantum espresso: 
!recl = 2 * ngms * ngms * nfs
!recl = 1 * 10 * 10 * nfs
!recl = 20 
!unf_recl = 8 * recl
!open ( iuncoul, file = "./silicon.coul", iostat = ios, form = 'unformatted', &
!status = 'unknown', access = 'direct', recl = unf_recl)

irr=1

! Shuffle the grid so that we have the indices for k and k+q and can simply read out 
!the eigenvectors from the file generated in the NSCF step. 
!everything except q-> 0 requires grid_shuffle
!if (iq.gt.1) CALL grid_shuffle()
!WRITE(stdout,'(4x,"Number of frequencies", i4)') nfs

!HL putting this in so that we can test the pade evaluation on a denser grid.
    w_ryd(:) = wcoul(:)/RYTOEV

!LOOP OVER G
!DO ig = 1, ngmsig
DO ig = 1, 20
GOTO 124
    qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
      !write(6,*)qg2 
    if (qg2.gt.eps8) then
      do iw = 1, nfs
        drhoscfs(:,:) = (0.0d0, 0.0d0)
        dvbare(:) = (0.0d0, 0.0d0)
        dvbare (nl (ig) ) = (1.d0, 0.d0)
        WRITE(6, "")
        WRITE(stdout,'(4x,"Frequency= ",3f7.3)') fiu(iw)*RYTOEV
        WRITE(stdout,'(4x,"Screened Coulomb: q =",3f7.3,"  G =",3f7.3)') xq(:), g(:,ig)
        !From the inputcard we can choose whether to use a model dielectric or 
        !do the full sternheimer treatment.
        IF (modielec) then
           !The 'magic dielectric function'
           wwp    = 18.0/RYTOEV  ! plasma frequency in Ry
           eps0  = 11.4          ! static diel constant of Si
           q0    = 1.1           ! characteristic momentum of Si, a.u. from Resta

           qg = sqrt(tpiba2*qg2)
           fac = 1.d0/(1.d0-1.d0/eps0)
           wwq = wwp * sqrt ( fac * (1.d0 + (qg/eps0/q0)**2.d0 ) )


           !diagonal term ig = igp (all the others remain 0)
           drhoscfs (nl(ig), 1)  = 1.d0 - wwp**2.d0/((fiu(iw)+eta)**2.d0+wwq**2.d0)
           WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscfs(nl(ig),1) 
        ELSE

          call cft3 (dvbare, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)

          CALL solve_linter (dvbare, iw, drhoscfs)

          !hl- debug
          !calculating non scf system. one iteration of the linear system solver should print out
          !the dielectric constant.
          !call solve_linter_nonscf (igpert, drhoscfs)
          !hl- debug
          !scrcoul_g(:) = 0.0d0
          !scrcoul (ig,igp,iw) = (0.0d0, 0.0d0)

             IF (convt .ne. .true.) WRITE( stdout, '(/,5x,"No convergence ")')

             call cft3 (drhoscfs, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
             call cft3 (dvbare, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)

!       WRITE(stdout, '(4x,4x,"dvbare = ", 2f9.5)'), dvbare(nl (ig)) 
        WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscfs(nl(ig),1) + dvbare(nl (ig) ) 
        ENDIF

! Generate an entire row of the Screened Coulomb Matrix with 4pi*e2*|q+G|^{-2}eps_{GG'}. 
        DO igp = 1, ngmsig
         !scrcoul (ig,igp,iw,nspin_mag) = drhoscfs ( nl(igp), 1 ) * dcmplx ( e2 * fpi / (tpiba2*qg2), 0.0d0 )
         scrcoul (ig,igp,iw,nspin_mag) = drhoscfs ( nl(igp), 1 ) !* dcmplx ( e2 * fpi / (tpiba2*qg2), 0.0d0 )
        ENDDO

! Spencer/Alavi truncation of the bare coulomb interaction
! [PRB 77,193110 (2008]
! rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
! Now need to think about spencer Alavi Truncation on a reduced q grid... no longer nq1*nq2*nq3
! rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))

        rcut = (float(3)/float(4)/pi*omega*float(216))**(float(1)/float(3))
        qg = sqrt( (g(1,ig )+xq(1))**2.d0 + (g(2,ig )+xq(2))**2.d0 + (g(3,ig )+ xq(3))**2.d0 )

!HL sqrt(tpiba2) should just redefine this as a variable.

        spal = 1.0d0 - cos ( rcut * sqrt(tpiba2) * qg )

!SAME DEAL AS IN SGW I've been looking at eps^{-1} and not testing W. However the W for q->0 G=0 is close to
!  Also when using the same nqs = 216 I get identical spal values and W. 

        do igp = 1, ngmsig
!        scrcoul (ig,igp,iw,nspin_mag) = scrcoul (ig,igp,iw,nspin_mag) * dcmplx ( spal, 0.0d0)
         scrcoul (ig,igp,iw,nspin_mag) = scrcoul (ig,igp,iw,nspin_mag) 
        enddo

! Currently writing elements of inveps here should do a final
! write to disk at the end of the calculation outside coulomb routine.
      if (iw.eq.1) then
        do igp = 1, 20
           if(igp.eq.ig) then
           scrcoul_g(igp) = real(scrcoul(ig, igp, iw, nspin_mag)) + 1.0d0
           else
           scrcoul_g(igp) = real(scrcoul(ig, igp, iw, nspin_mag))
           endif
        enddo     
        write (iuncoul, rec = ig, iostat = ios) scrcoul_g
      endif  
!        rec_code=20
!        CALL write_rec('done_drhod',1,0.0_DP,-1000,.false.,1,&
!                        drhoscfs)
      enddo !iw

!    Analytical Continuation to the real axis.
124 CONTINUE
!@10TION only diagonal elements when using the model-dielectric function.
!   do igp = 1, ngmsig
     do igp = 1, 20
!     do igp = ig, ig
!     Pade input points on the imaginary axis

scrcoul (ig,igp,1,nspin_mag) = (0.07852, 0.00000)
scrcoul (ig,igp,2,nspin_mag) = (0.08350, 0.00000)
scrcoul (ig,igp,3,nspin_mag) = (0.09746, 0.00000)
scrcoul (ig,igp,4,nspin_mag) = (0.14496, -0.00008)
scrcoul (ig,igp,5,nspin_mag) = (0.17541, 0.00000)
scrcoul (ig,igp,6,nspin_mag) = (0.35131, 0.00001)
scrcoul (ig,igp,7,nspin_mag) = (0.51413, 0.00000)
scrcoul (ig,igp,8,nspin_mag) = (0.63758, 0.00000)
scrcoul (ig,igp,9,nspin_mag) = (0.86489, 0.00000)

          do iw = 1, nfs
            z(iw) = dcmplx( 0.d0, fiu(iw))
            u(iw) = scrcoul (ig,igp,iw,nspin_mag)
            write(6,*) u(iw)
          enddo

!        Pade coefficients
         call pade_coeff ( nfs, z, u, a)

!        Overwrite scrcoul with Pade coefficients to be read in gw_product.
          do iw = 1, nfs
            scrcoul (ig,igp,iw,nspin_mag) = a(iw)
          enddo

! HL THIS IS FOR TESTING ONLY
!        Pade output points on the real axis (at a distance eta)
!        (I use the same grid for simplicity - this can be changed)

    if (igp.eq.1) then

    write(6,'(3f9.5)') nwcoul

    do iw = 1, nwcoul
      call pade_eval ( nfs, z, a, dcmplx( w_ryd(iw), eta), scrcoul (ig,igp,iw,nspin_mag))
    enddo

    do iw = 1, nwcoul
       write(6,'(4x,"PADE INVEPS: ",3f9.5)') w_ryd(iw)*RYTOEV, scrcoul (1,1,iw,nspin_mag)
    enddo
    endif
    STOP

     enddo !enddo on igp
    endif !on qg2 > 0 
ENDDO !on G

      !write ( iuncoul, rec = iq, iostat = ios) scrcoul
        CALL davcio(scrcoul, lrcoul, iuncoul, iq, 1 )

        tcpu = get_clock ('GW')
  DEALLOCATE (drhoscfs)
  CALL stop_clock ('coulomb')
RETURN
END SUBROUTINE coulomb

! Copyright (C) 2001-2008 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
SUBROUTINE coulomb_ms(iq, igstart, igstop, scrcoul) 
!-----------------------------------------------------------------------
! This subroutine is the main driver of the COULOMB self consistent cycle
! which calculates the dielectric matrix by generating the density response
! to a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G).
!The dielectric matrix is given by:
!eps_{q}^{-1}(G,G',iw) = (\delta_{GG'} + drhoscfs^{scf}_{G,G',iw})

  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE gvect,      ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,    ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms, doublegrid
  USE constants,  ONLY : e2, fpi, RYTOEV, pi, eps8
  USE cell_base,  ONLY : alat, tpiba2, omega
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout, ionode
  USE uspp,       ONLY: okvan
  USE control_gw, ONLY : zue, convt, rec_code, modielec, eta, godbyneeds, padecont, solve_direct
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
  COMPLEX(DP), allocatable :: drhoscf(:,:,:)
  COMPLEX(DP) :: padapp, w
!HL temp variable for scrcoul to write to file.  
  COMPLEX(DP) :: cw
  INTEGER :: unf_recl, recl, ios
  INTEGER :: iq, screening 
  LOGICAL :: exst
!again should decide if this should be allocated globally. 
  COMPLEX(DP) :: scrcoul(ngmpol, ngmpol, nfs, 1)
!modeps and spencer-alavi vars
  REAL(DP) :: wwp, eps0, q0, wwq, fac
  REAL(DP) :: qg, rcut, spal
! used to test the recover file
  EXTERNAL get_clock
  CALL start_clock ('coulomb')

ALLOCATE (drhoscf(nrxx, nspin_mag, nfs))    
irr=1
scrcoul(:,:,:,:) = (0.d0, 0.0d0)

!LOOP OVER ig, unique g vectors only. 
!g is sorted in magnitude order.
DO ig = igstart, igstop
   qg2 = (g(1,ig_unique(ig))+xq(1))**2 + (g(2,ig_unique(ig))+xq(2))**2 + (g(3,ig_unique(ig))+xq(3))**2

      drhoscf(:,:,:) = dcmplx(0.0d0, 0.0d0)
      dvbare(:)      = dcmplx(0.0d0, 0.0d0)
      dvbare (nls(ig_unique(ig)) ) = dcmplx(1.d0, 0.d0)

      call cft3s (dvbare, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +1)
      CALL solve_direct_shift (dvbare, iw, drhoscf)
      call cft3s (dvbare, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)

      DO iw = 1, nfs
         call cft3  (drhoscf(:,:,iw), nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)

         if(iq.eq.1) then
            WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) = ", 2f9.5)'), drhoscf(nl(ig_unique(ig)), 1, iw) + dvbare(nls(ig_unique(ig)))
         endif

         DO igp = 1, ngmpol
            if(igp.ne.ig_unique(ig)) then
               scrcoul(ig_unique(ig), igp, iw, nspin_mag) = drhoscf(nl(igp),1, iw)
            else
               scrcoul(ig_unique(ig), igp, iw, nspin_mag) = drhoscf(nl(igp),1, iw) + dvbare(nls(ig_unique(ig)))
            endif
         ENDDO
      ENDDO !iw
ENDDO 

tcpu = get_clock ('GW')
DEALLOCATE (drhoscfs)
CALL stop_clock ('coulomb')
RETURN

END SUBROUTINE coulomb_ms
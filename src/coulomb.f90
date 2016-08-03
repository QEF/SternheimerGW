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
!> This subroutine is the main driver of the COULOMB self consistent cycle
!! which calculates the dielectric matrix by generating the density response
!! to a charge dvbare(nl(ig)) = 1.00 + i*0.00 at a single fourier component (G).
!! The dielectric matrix is given by:
!! eps_{q}^{-1}(G,G',iw) = (\delta_{GG'} + drhoscfs^{scf}_{G,G',iw})
!-----------------------------------------------------------------------
SUBROUTINE coulomb(iq, igstart, num_task, scrcoul) 
!-----------------------------------------------------------------------
  USE cell_base,        ONLY : alat, tpiba2, omega
  USE constants,        ONLY : e2, fpi, RYTOEV, pi, eps8
  USE control_gw,       ONLY : zue, convt, rec_code, modielec, eta, godbyneeds, padecont,&
                               solve_direct, do_epsil, do_q0_only, tinvert, niter_gw
  USE disp,             ONLY : nqs, nq1, nq2, nq3
  USE eqv_gw,           ONLY : drhoscfs, dvbare
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE freq_gw,          ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
  USE gvecs,            ONLY : nls
  USE gvect,            ONLY : ngm, g, nl
  USE gwsigma,          ONLY : sigma_c_st, gcutcorr
  USE gwsymm,           ONLY : ig_unique, ngmunique
  USE io_global,        ONLY : stdout, ionode
  USE ions_base,        ONLY : nat
  USE lsda_mod,         ONLY : nspin
  USE mp,               ONLY : mp_sum, mp_barrier
  USE mp_global,        ONLY : inter_image_comm, intra_image_comm, &
                               my_image_id, nimage, root_image
  USE mp_pools,         ONLY : me_pool, root_pool, inter_pool_comm
  USE mp_world,         ONLY : mpime
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE kinds,            ONLY : DP
  USE klist,            ONLY : lgauss
  USE partial,          ONLY : done_irr, comp_irr
  USE paw_variables,    ONLY : okpaw
  USE qpoint,           ONLY : xq
  USE solve_module,     ONLY : solve_linter
  USE units_gw,         ONLY : iuncoul, lrcoul
  USE uspp,             ONLY : okvan
  USE uspp_param,       ONLY : nhm

  IMPLICIT NONE

  !> index of the active q-point
  INTEGER, INTENT(IN) :: iq

  !> first index of the G vector evaluated on this process
  INTEGER, INTENT(IN) :: igstart

  !> number of G vector evaluated by this process
  INTEGER, INTENT(IN) :: num_task

  !> the screened coulomb interaction
  COMPLEX(dp), INTENT(OUT) :: scrcoul(gcutcorr, nfs, num_task)

  !> actual index inside the array
  INTEGER :: indx

  !> complex constant of zero
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND = dp)

  !> complex constant of one
  COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND = dp)

  REAL(DP) :: tcpu, get_clock
! timing variables
  REAL(DP) :: qg2, qg2coul
  INTEGER :: ig, igp, iw, npe, icounter
  INTEGER :: igpert, isp
  COMPLEX(DP), allocatable :: drhoaux (:,:) 
  COMPLEX(DP) :: padapp, w
!HL temp variable for scrcoul to write to file.  
  COMPLEX(DP) :: cw
  INTEGER :: unf_recl, recl, ios
  INTEGER :: screening 
  LOGICAL :: exst
!modeps and spencer-alavi vars
  REAL(DP) :: wwp, eps0, q0, wwq, fac
  REAL(DP) :: qg, rcut, spal
! used to test the recover file
  EXTERNAL get_clock

  ! we use the frequency as equivalent of the perturbation in phonon
  ALLOCATE (drhoscfs(dfftp%nnr, nspin_mag, nfs))
  IF (nspin_mag /= 1) CALL errore(__FILE__, "magnetic calculation not implemented", 1)
  scrcoul = (0.d0, 0.0d0)

!LOOP OVER ig, unique g vectors only. 
!g is sorted in magnitude order.
!WRITE(1000+mpime, '(2i4)') igstart, igstop
!WRITE(1000+mpime, *) ig_unique(:)
DO indx = 1, num_task
  ig = igstart + indx - 1
!     if (do_q0_only.and.ig.gt.1) CYCLE
      qg2 = (g(1,ig_unique(ig))+xq(1))**2+(g(2,ig_unique(ig))+xq(2))**2+(g(3,ig_unique(ig))+xq(3))**2
      IF(solve_direct) THEN
         ! q + G = 0 is treated differently
         IF(tinvert .AND. qg2 < eps8) CYCLE 
         drhoscfs = zero 
         dvbare   = zero 
         dvbare(nls(ig_unique(ig))) = one
         CALL invfft('Smooth', dvbare, dffts)
         CALL solve_linter(1, dvbare, fiu(:nfs), drhoscfs)
         CALL fwfft('Smooth', dvbare, dffts)
         do iw = 1, nfs
            CALL fwfft ('Dense', drhoscfs(:,1,iw), dffts)
            WRITE(stdout, '(4x,4x,"eps_{GG}(q,w) = ", 2f10.4)') &
              drhoscfs(nls(ig_unique(ig)),1,iw)+dvbare(nls(ig_unique(ig)))
            do igp = 1, gcutcorr
               if(igp.ne.ig_unique(ig)) then
!diagonal elements drho(G,G').
                  scrcoul(igp, iw, indx) = drhoscfs(nls(igp), 1, iw)
               else
!diagonal elements eps(\G,\G') = \delta(G,G') - drho(G,G').
                  scrcoul(igp, iw, indx) = drhoscfs(nls(igp), 1, iw) + dvbare(nls(ig_unique(ig)))
               endif
            enddo
         enddo !iw
      ELSE
        IF (qg2 < eps8) CYCLE
        if(qg2.lt.0.001.AND.lgauss) then 
          write(stdout,'("Not calculating static electric field applied to metal, cycling coulomb")')
          WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) =   0.000000   0.0000000")')
          CYCLE
        endif
        drhoscfs = zero
        dvbare   = zero
        dvbare(nls(ig_unique(ig))) = one
        CALL invfft('Smooth', dvbare, dffts)
        CALL solve_linter(niter_gw, dvbare, fiu(:nfs), drhoscfs)
        CALL fwfft('Smooth', dvbare, dffts)
        DO iw = 1, nfs
           CALL fwfft('Dense', drhoscfs(:,1,iw), dffts)
           IF(ionode) THEN
             WRITE(stdout, '(4x,4x,"inveps_{GG}(q,w) = ", 2f12.5)') &
               drhoscfs(nls(ig_unique(ig)), 1, iw) + dvbare(nls(ig_unique(ig)))
             DO igp = 1, gcutcorr
               scrcoul(igp, iw, indx) = drhoscfs(nl(igp), 1, iw)
             END DO ! igp
           END IF
        END DO ! iw
      END IF !solve_direct/solve_linter
ENDDO 
545 CONTINUE
tcpu = get_clock ('GW')
DEALLOCATE (drhoscfs)
RETURN
END SUBROUTINE coulomb

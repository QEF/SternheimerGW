!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
!> Driver routine for the solution of the linear system.
!!
!! It defines the change of the wavefunction due to a perturbing potential
!! parameterized in r and iw, i.e. \f$v_{r, iw} (r')\f$
!! Currently symmetrized in terms of mode etc. Might need to strip this out
!! and check PW for how it stores/symmetrizes charge densities.
!!
!! This routine performs the following tasks
!! 1. computes the bare potential term Delta V | psi > 
!! 2. adds to it the screening term Delta V_{SCF} | psi >
!! 3. applies P_c^+ (orthogonalization to valence states)
!! 4. calls c_bi_cgsolve_all to solve the linear system
!! 5. computes Delta rho, Delta V_{SCF}.
!----------------------------------------------------------------------------
SUBROUTINE solve_linter(dvbarein, iw, drhoscf)
!-----------------------------------------------------------------------------

  USE buffers,              ONLY : save_buffer, get_buffer
  USE check_stop,           ONLY : check_stop_now
  USE constants,            ONLY : degspin, eps8
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   lgamma, convt, nbnd_occ, alpha_mix, &
                                   rec_code_read, where_rec, flmixdpot, &
                                   high_io, maxter_coul
  USE eqv,                  ONLY : dvpsi, dpsi, evq
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE freq_gw,              ONLY : fiu, nfsmax
  USE gvecs,                ONLY : doublegrid
  USE gvect,                ONLY : nl
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, wk, nkstot, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_world,             ONLY : mpime
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE qpoint,               ONLY : npwq, nksq, igkq, ikks, ikqs
  USE timing_module,        ONLY : time_coul_solver
  USE units_gw,             ONLY : lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iudwfm, iudwfp 
  USE uspp,                 ONLY : vkb
  USE uspp_param,           ONLY : nhm
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, et

  IMPLICIT NONE

  !> the initial perturbuing potential
  COMPLEX(dp), INTENT(IN)  :: dvbarein (dffts%nnr)

  !> the index of the frequency
  INTEGER,     INTENT(IN)  :: iw

  !> the change of the scf charge
  COMPLEX(dp), INTENT(OUT) :: drhoscf(dfftp%nnr, nspin_mag)

  complex(DP), allocatable, target :: dvscfin(:,:,:)
  ! change of the scf potential 
  complex(DP), pointer :: dvscfins (:,:,:)
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:), dvscfout (:,:)
  complex(DP) :: dpsip(npwx*npol, nbnd), dpsim(npwx*npol, nbnd)
  complex(DP), allocatable :: dbecsum (:,:,:), dbecsum_nc(:,:,:,:), aux1 (:,:)
  complex(DP) :: cw
  complex(DP), allocatable :: etc(:,:)
  complex(kind=DP) :: ZDOTC
  ! Misc variables for metals
  real(DP), external :: w0gauss, wgauss
  real(DP) , allocatable :: h_diag (:,:)
  real(DP) :: thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: weight, aux_avg (2)
  !HL dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  real(DP) :: tcpu, get_clock ! timing variables
  real(DP) :: meandvb

  integer :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             is,         & ! counter on spin polarizations
             nrec, nrec1,& ! the record number for dvpsi and dpsi
             lmres         ! number of gmres iterations to include when using bicgstabl.
  integer :: irr, imode0, npe

  external ch_psi_all, cg_psi, h_psi_all, ch_psi_all_nopv
  real(kind=DP)    :: DZNRM2
  external ZDOTC, DZNRM2

  logical :: conv_root    ! true if linear system is converged

  !> special treatment is possible for the first iteration
  LOGICAL first_iteration

  !> the number of frequencies
  INTEGER num_freq

  !> loop variable for frequencies
  INTEGER ifreq

  !> complex value of 0
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND = dp)

  !> complex value of 0.5
  COMPLEX(dp), PARAMETER :: half = CMPLX(0.5_dp, 0.0_dp, KIND = dp)

  CALL start_clock (time_coul_solver)

  ! note: for now we just evaluate a single frequency at a time
  num_freq = 1

  ALLOCATE (dpsi(npwx*npol, nbnd))
 
  IF (rec_code_read > 20 ) RETURN

!HL- Allocate arrays for dV_scf (need to alter these from (dfftp%nnr, nspin_mag, npe) 
!to just (dfftp%nnr, nspin_mag).
  npe    = 1
  imode0 = 1
  irr    = 1
  ipert  = 1
  lter   = 0
  lmres  = 1

  ALLOCATE(dvscfout(dfftp%nnr, nspin_mag))  
  ALLOCATE(drhoscfh( dfftp%nnr, nspin_mag))
  ALLOCATE(etc(nbnd, nkstot))
  ALLOCATE(dbecsum((nhm * (nhm + 1)) / 2, nat, nspin_mag))

  ALLOCATE(dvscfin(dfftp%nnr, nspin_mag, num_freq))
  IF (doublegrid) THEN
    ALLOCATE(dvscfins(dffts%nnr, nspin_mag, num_freq))
  ELSE
    dvscfins => dvscfin
  ENDIF

!Complex eigenvalues
  IF (noncolin) ALLOCATE (dbecsum_nc(nhm, nhm, nat, nspin))
  ALLOCATE(aux1(dffts%nnr, npol))
  ALLOCATE(h_diag(npwx*npol, nbnd))

  iter0 =  0
  convt = .FALSE.
  where_rec = 'no_recover'

! The outside loop is over the iterations.
! niter_gw := maximum number of iterations
  
  DO kter = 1, niter_gw
    !
    first_iteration = (kter == 1)
    !
    iter       = kter + iter0
    ltaver     = 0
    lintercall = 0
    drhoscf    = zero
    drhoscfh   = zero
    dbecsum    = zero
    IF (noncolin) dbecsum_nc = zero 
    !
    DO ik = 1, nksq
      !
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      ! this is necessary for now until all instances of igkq
      ! are replaced with the appropriate igk_k
      igkq = igk_k(:,ikq)
      !
      IF (lsda) current_spin = isk(ikk)
      !
      ! read unperturbed wavefunctions psi(k) and psi(k+q)
      !
      IF (nksq > 1) THEN
        IF (lgamma) THEN
          CALL get_buffer(evc, lrwfc, iuwfc, ikk)
        ELSE
          CALL get_buffer(evc, lrwfc, iuwfc, ikk)
          CALL get_buffer(evq, lrwfc, iuwfc, ikq)
        END IF
      END IF
      !
      ! compute beta functions and kinetic energy for k-point ikq
      ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
      !
      CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb)
      CALL g2_kin(ikq)
      !
      ! compute preconditioning matrix h_diag used by cgsolve_all
      !
      CALL h_prec(ik, evq, h_diag)
      !
      ! in the first iteration we initialize the linear system
      ! and may use the multishift solver
      IF (first_iteration) THEN
        !
        !  At the first iteration dvbare_q*psi_kpoint is calculated
        !  and written to file
        !
        nrec = ik
        CALL dvqpsi_us(dvbarein, ik, .FALSE.)
        CALL save_buffer(dvpsi, lrbar, iubar, nrec)
        !
        ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
        ! Apply -P_c^+.
        !-P_c^ = - (1-P_v^):
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
        !
        ! At the first iteration initialize arrays to zero
        !
        dpsi     = zero
        dpsim    = zero
        dpsip    = zero
        dvscfin  = zero
        dvscfout = zero
        !
        ! starting threshold for iterative solution of the linear system
        !
        thresh = 1.0d-2
        !
        ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
        ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
        ! TODO replace with multishift solver
        !
        conv_root = .true.
        etc = CMPLX(et, 0.0_dp, kind = dp)
        cw  = fiu(iw) 
        !
        IF (iw == 1) THEN
          !
          ! for the Fermi energy, we may use the CG solver
          !
          CALL cgsolve_all(h_psi_all, cg_psi, et(1,ikk), dvpsi, dpsip, h_diag, & 
               npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),  &
               npol)
          CALL ZCOPY(npwx * npol * nbnd_occ(ikk), dpsip, 1, dpsi, 1)
          !
        ELSE 
          !
          ! for all other frequencies we use the BiCG solver
          !
          conv_root = .TRUE.
          CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsip, h_diag, &
               npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),   &
               npol, +cw, maxter_coul, .TRUE.)
          !
          conv_root = .TRUE.
          CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsim, h_diag, &
               npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),   &
               npol, -cw, maxter_coul, .TRUE.)
          !
          dpsi = zero
          CALL ZAXPY(npwx * npol * nbnd_occ(ikk), half, dpsim, 1, dpsi, 1)
          CALL ZAXPY(npwx * npol * nbnd_occ(ikk), half, dpsip, 1, dpsi, 1)
          !
        END IF ! frequency = Fermi energy?
        !
      ELSE ! general case iter > 1
        !
        ! After the first iteration dvbare_q*psi_kpoint is read from file
        !
        nrec = ik
        CALL get_buffer(dvpsi, lrbar, iubar, nrec)
        !
        ! frequencies are the equivalent of perturbations in PHonon
        DO ifreq = 1, num_freq
          !
          ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
          ! dvscf_q from previous iteration (mix_potential)
          !
          DO ibnd = 1, nbnd_occ (ikk)
            CALL cft_wave(ik, evc (1, ibnd), aux1, +1)
            CALL apply_dpot(dffts%nnr, aux1, dvscfins(1, 1, ifreq), current_spin)
            CALL cft_wave(ik, dvpsi (1, ibnd), aux1, -1)
          ENDDO
          !
          !  In the case of US pseudopotentials there is an additional
          !  selfconsist term which comes from the dependence of D on
          !  V_{eff} on the bare change of the potential
          !
          CALL adddvscf(ifreq, ik)
          !
          ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
          ! Apply -P_c^+.
          !-P_c^ = - (1-P_v^):
          CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
          !
          ! starting value for delta_psi is read from iudwf
          !
          IF (high_io) THEN
            CALL get_buffer(dpsip, lrdwf, iudwfp, ik)
            ! TODO replace iw with ifreq
            IF (iw > 1) CALL get_buffer(dpsim, lrdwf, iudwfm, ik)
          !
          ! restart from default
          !
          ELSE
            dpsim = zero
            dpsip = zero 
          ENDIF
          !
          ! threshold for iterative solution of the linear system
          !
          thresh = MIN(1.d-1 * SQRT(dr2), 1.d-2)
          !
          ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
          ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
          !
          conv_root = .true.
          etc = CMPLX(et, 0.0_dp, kind = dp)
          ! TODO update to ifreq
          cw  = fiu(iw) 
          !
          IF (iw == 1) THEN
            !
            ! for the Fermi energy, we may use the CG solver
            !
            CALL cgsolve_all(h_psi_all, cg_psi, et(1,ikk), dvpsi, dpsip, h_diag, & 
                 npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),  &
                 npol)
            CALL ZCOPY(npwx * npol * nbnd_occ(ikk), dpsip, 1, dpsi, 1)
            !
          ELSE 
            !
            ! for all other frequencies we use the BiCG solver
            !
            conv_root = .TRUE.
            CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsip, h_diag, &
                 npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),   &
                 npol, +cw, maxter_coul, .TRUE.)
            !
            conv_root = .TRUE.
            CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsim, h_diag, &
                 npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),   &
                 npol, -cw, maxter_coul, .TRUE.)
            !
            dpsi = zero
            CALL ZAXPY(npwx * npol * nbnd_occ(ikk), half, dpsim, 1, dpsi, 1)
            CALL ZAXPY(npwx * npol * nbnd_occ(ikk), half, dpsip, 1, dpsi, 1)
            !
          END IF ! frequency = Fermi energy?

        END DO ! ifreq

      END IF ! first_iteration

        ltaver = ltaver + lter
        lintercall = lintercall + 1

!        WRITE(1000+mpime, '(5x,"kpoint ",i4,"  ibnd ",i4, &
!          &               " solve_linter:  ", e10.3 , "iter ", i4)')  &
!          &                 ik , nbnd_occ(ikk), anorm, iter

        IF (.NOT.conv_root) WRITE(1000+mpime, '(5x,"kpoint ",i4,"  ibnd ",i4, &
                  &              " solve_linter: root not converged ", e10.3 , "iter ", i4)')  &
                  &                ik , nbnd_occ(ikk), anorm, iter

           nrec1 =  ik
         !calculates dvscf, sum over k => dvscf_q_ipert
         !incdrhoscf:  This routine computes the change of the charge density due to the
         !HL low/io
          if(high_io) then
             if(iw.gt.1) call save_buffer (dpsim, lrdwf, iudwfm, ik)
             call save_buffer (dpsip, lrdwf, iudwfp, ik)
          endif
         ! perturbation. It is called at the end of the computation of the
         ! change of the wavefunction for a given k point.
           weight = wk (ikk)

         IF (noncolin) THEN
           call incdrhoscf_w (drhoscf(1, 1) , weight, ik, &
                              dbecsum(1,1,current_spin), dpsi)
         ELSE
           call incdrhoscf_w (drhoscf(1, 1) , weight, ik, &
                              dbecsum(1,1,current_spin), dpsi(1,1))
         ENDIF

     enddo !on k-points

     call mp_sum ( drhoscf, inter_pool_comm )

     if (doublegrid) then
         do is = 1, nspin_mag
            call cinterpolate (drhoscfh(1,1), drhoscf(1,1), 1)
         enddo
     else
            call zcopy (nspin_mag*dfftp%nnr, drhoscf(1,1), 1, drhoscfh(1,1), 1)
     endif
!
! In the noncolinear, spin-orbit case rotate dbecsum
! IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, npe)
! Now we compute for all perturbations the total charge and potential
!
! HL NEED FOR ULTRASOFT:
!    call addusddens (drhoscfh, dbecsum, 0)
     call zcopy (dfftp%nnr*nspin_mag, drhoscfh(1,1), 1, dvscfout(1,1),1)

     meandvb = sqrt ((sum(real(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0 )/float(dffts%nnr)
     if (meandvb.lt.1.d-10) then 
         CALL fwfft ('Dense', dvscfout(:,1), dfftp)
         dvscfout ( nl(1), current_spin ) = (0.d0, 0.0d0)
         CALL invfft ('Dense', dvscfout(:,1), dfftp)
     endif
!
! IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
! for q->0 the Fermi level can shift.
! IF (lmetq0) call ef_shift(drhoscfh,ldos,ldoss, dos_ef, irr,npe,.false.)
     call dv_of_drho (1, dvscfout(1,1), .true.)
!    nmix_gw = 4
     if (iw.eq.1) then
!Density reponse in real space should be real at zero freq no matter what!
!just using standard broyden for the zero freq. case.
        call mix_potential(2*dfftp%nnr*nspin_mag, dvscfout, dvscfin, alpha_mix(kter), &
                           dr2, tr2_gw, iter, nmix_gw, flmixdpot, convt)
     else
    !Is the hermitian mixing scheme still okay?
        call mix_potential_c(dfftp%nnr*nspin_mag, dvscfout(1,1), dvscfin, alpha_mix(kter),& 
                           dr2, tr2_gw, iter, nmix_gw, convt)
                             
     endif
  !HL need this if we aren't dealing with electric fields:
  !lmetq0 = lgauss.and.lgamma
  !if (lmetq0) then
  !      allocate ( ldos ( dfftp%nnr  , nspin_mag) )
  !      allocate ( ldoss( dffts%nnr , nspin_mag) )
  !      call localdos ( ldos , ldoss , dos_ef )
  !endif
     if (doublegrid) then
           do is = 1, nspin_mag
              call cinterpolate (dvscfin(1,is,1), dvscfins(1,is,1), -1)
           enddo
     endif

#ifdef __PARA
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)
     call mp_sum ( aux_avg, inter_pool_comm )
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif
     tcpu = get_clock ('SGW')
     dr2 = dr2 
     FLUSH(stdout)
     rec_code=10

!     WRITE(1000+mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!           "secs   av.it.: ",f5.1)') iter, tcpu, averlt
!     WRITE(1000+mpime, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
!           &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
     if (convt) goto 155

  enddo !loop on kter (iterations)

  ! abort if solver doesn't converge
  CALL errore(__FILE__, "Iterative solver did not converge within given number&
                       & of iterations", niter_gw)

155 iter0=0

   WRITE( stdout, '(/,5x," iter # ",i4," total cpu time :",f8.1, &
        & "secs av.it.:",f5.1)') iter, tcpu, averlt
!   WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
!          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
!   WRITE(1000+mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!        "secs   av.it.: ",f5.1)') iter, tcpu, averlt

  drhoscf(:,:) = dvscfin(:,:,1)
  
 !call mp_barrier(intra_image_comm)
 
  deallocate (h_diag)
  deallocate (aux1)
  deallocate (dbecsum)
  IF (noncolin) deallocate (dbecsum_nc)
  deallocate (dvscfout)
  deallocate (drhoscfh)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  deallocate (dpsi)

  CALL stop_clock (time_coul_solver)

END SUBROUTINE solve_linter

SUBROUTINE check_all_convt(convt)
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : nproc_image, me_image, intra_image_comm
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: convt
  INTEGER,ALLOCATABLE :: convt_check(:)
  !
  IF(nproc_image==1) RETURN
  !
  ALLOCATE(convt_check(nproc_image+1))
  !
  convt_check = 1
  IF(convt) convt_check(me_image+1) = 0
  !
  CALL mp_sum(convt_check, intra_image_comm)
  !CALL mp_sum(ios, inter_pool_comm)
  !CALL mp_sum(ios, intra_bgrp_comm)
  !
!  convt = ALL(convt_check==0)
  IF(ANY(convt_check==0).and..not.ALL(convt_check==0) ) THEN
    CALL errore('check_all_convt', 'Only some processors converged: '&
               &' something is wrong with solve_linter', 1)
  ENDIF
  !
  DEALLOCATE(convt_check)
  !
END SUBROUTINE

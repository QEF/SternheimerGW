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
!
! Significantly modified to use Multishift Linear System Solver
! Henry Lambert
!-------------------------------------------------------------------------------
SUBROUTINE solve_lindir(dvbarein, drhoscf)
!-----------------------------------------------------------------------------
  !  HL
  !  Driver routine for the solution of the linear system which
  !  defines the change of the wavefunction due to a perturbing potential
  !  parameterized in r and iw. i.e. v_{r, iw} (r')
  !  It performs the following tasks:
  !   a) computes the bare potential term Delta V | psi > 
  !   b) adds to it the screening term Delta V_{SCF} | psi >
  !   c) applies P_c^+ (orthogonalization to valence states)
  !   d) calls c_bi_cgsolve_all to solve the linear system
  !   e) computes Delta rho, Delta V_{SCF}.
  !   Currently symmetrized in terms of mode etc. Might need to strip this out
  !   and check PW for how it stores/symmetrizes charge densities.
!----------------------------------------------------------------------------
  USE buffers,              ONLY : get_buffer
  USE cell_base,            ONLY : tpiba2,at,bg
  USE check_stop,           ONLY : check_stop_now
  USE constants,            ONLY : degspin, eps8
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, flmixdpot, current_iq, &
                                   ext_recover, eta, maxter_coul, maxter_green, prec_direct, &
                                   prec_shift
  USE ener,                 ONLY : ef
  USE eqv,                  ONLY : dvpsi, evq
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE freq_gw,              ONLY : fpol, fiu, nfs, nfsmax
  USE gvecs,                ONLY : nls
  USE gvect,                ONLY : ngm, g, nl
  USE gwsigma,              ONLY : sigma_c_st, ecutsco, ecutprec
  USE io_files,             ONLY : prefix
  USE io_global,            ONLY : stdout, ionode
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE kinds,                ONLY : DP
  USE klist,                ONLY : lgauss, degauss, ngauss, xk, wk, nkstot, nks, ngk, igk_k
  USE lr_symm_base,         ONLY : minus_q, irgq, nsymq, rtau 
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE mp_world,             ONLY : mpime
  USE mp_pools,             ONLY : inter_pool_comm
  USE nlcc_gw,              ONLY : nlcc_any
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE output_mod,           ONLY : fildrho, fildvscf
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE scf,                  ONLY : rho
  USE spin_orb,             ONLY : domag
  USE timing_module,        ONLY : time_solver
  USE units_gw,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iunrec, iudvscf, iudwfm, iudwfp 
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin, current_k, et

  implicit none
  !
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes
  complex(DP) :: drhoscf (dfftp%nnr, nfs)
  ! output: the change of the scf charge
  complex(DP) :: dvbarein (dffts%nnr)
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable  ::  dvscfout (:,:)
  complex(DP)               :: alphabeta(2,nbnd,maxter_coul+1)
  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  complex(DP), allocatable  :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
                               dbecsum (:,:,:)
  complex(kind=DP)         :: ZDOTC
  complex(DP)              :: cw
  complex(DP), allocatable :: etc(:,:)

  ! HL careful now... complexifying preconditioner:
  real(DP) , allocatable :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian
  ! complex(DP) , allocatable :: h_diag (:,:)
  real(DP) :: thresh, anorm, averlt, dr2
  real(DP) :: x, a, b
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight, aux_avg (2)
  ! Misc variables for metals
  ! dos_ef: density of states at Ef
  real(DP), external :: w0gauss, wgauss
  ! For approx, mixing scheme.
  real(kind=DP) :: DZNRM2
  real(DP) :: tcpu, get_clock ! timing variables
  real(DP) :: meandvb
  !
  ! counter on frequencies.
  !
  integer :: iw, allocatestatus
  integer :: irr

  
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
             ig,         & ! counter on G vectors
             ndim,       & ! integer actual row dimension of dpsi
             is,         & ! counter on spin polarizations
             nt,         & ! counter on types
             na,         & ! counter on atoms
             nrec, nrec1,& ! the record number for dvpsi and dpsi
             ios,        & ! integer variable for I/O control
             mode,       & ! mode index
             igpert,     & ! bare perturbation g vector.
             lmres         ! number of gmres iterations to include when using bicgstabl.
  integer  :: gveccount
  integer  :: niters(nbnd)
  logical :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0      ! true if xq=(0,0,0) in a metal
  external ZDOTC, DZNRM2
  external cg_psi, ch_psi_all, h_psi_all, ch_psi_all_green
  !complex(DP)               :: dpsit(npwx, nbnd, nfs), dpsi(npwx,nbnd,nfs)
  !complex(DP), allocatable  :: dpsic(:,:,:)
  complex(DP)               :: dpsi(npwx,nbnd,nfs)
  complex(DP)               :: dpsipm(npwx, nbnd, 2*nfs-1)

  CALL start_clock (time_solver)

  if (rec_code_read > 20 ) RETURN
  irr    = 1
  ipert  = 1
  lter   = 0
  lmres  = 1
  allocate (dvscfout ( dfftp%nnr , nfs))    
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat, nspin_mag))    
!Complex eigenvalues:
  allocate (etc(nbnd, nkstot))
  allocate (h_diag ( npwx*npol, nbnd))    
  IF (.not.ASSOCIATED(igkq)) ALLOCATE(igkq(npwx))
! if(.not.prec_direct) ALLOCATE (dpsic(npwx,nbnd,maxter_coul+1))
  iter0 = 0
  convt =.FALSE.
  where_rec='no_recover'
  if (iter0==-1000) iter0=0
!No self-consistency:
  do kter = 1, 1
     iter = kter + iter0
     ltaver = 0
     lintercall = 0
     drhoscf(:,:)   = (0.d0, 0.d0)
     dbecsum(:,:,:) = (0.d0, 0.d0)
     do ik = 1, nksq
        if (lgamma)  npwq = npw
!       ikks(ik) = 2 * ik - 1
!       ikqs(ik) = 2 * ik
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        igkq= igk_k(:,ikq)

        if (lsda) current_spin = isk (ikk)
       !Calculates beta functions (Kleinman-Bylander projectors), with
       !structure factor, for all atoms, in reciprocal space
       !HL the beta functions (vkb) are being generated properly.  
        call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
       !Reads unperturbed wavefuctions psi(k) and psi(k+q)
        if (nksq.gt.1) then
           if (lgamma) then
               call get_buffer (evc, lrwfc, iuwfc, ikk)
           else
               call get_buffer (evc, lrwfc, iuwfc, ikk)
               call get_buffer (evq, lrwfc, iuwfc, ikq)
           endif
        endif
        CALL g2_kin(ikq)
        !
        ! compute preconditioning matrix h_diag used by cgsolve_all
        !
        IF (prec_direct) THEN
          CALL h_prec (ik, evq, h_diag)
        ELSE
          h_diag = 0.d0
          DO ibnd = 1, nbnd_occ(ikk)
            h_diag(:npwq,ibnd) = 1.0
          END DO
        END IF
        !
        mode = 1
        nrec = ik
        call dvqpsi_us (dvbarein, ik, .false.)
! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
! Apply -P_c^+.
! -P_c^ = - (1-P_v^):
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi(:,:,1), npwq, .false. )
!        if(.not.prec_direct) dpsic(:,:,:)     =  dcmplx(0.d0, 0.d0)
!        dpsit(:,:,:)     =  dcmplx(0.d0, 0.d0)
        dpsi(:,:,:)      =  dcmplx(0.d0, 0.d0)
        alphabeta(:,:,:) =  dcmplx(0.d0, 0.d0)
        niters(:)        = 0  
        thresh    = tr2_gw
        conv_root = .true.
        etc(:,:)  = CMPLX(et(:,:), 0.0d0 , kind=DP)
   if(prec_direct) then
     do iw = 1, nfs
        cw    = fiu(iw)
        !first frequency in list should be zero.
        if((real(cw).eq.0.0d0).and.(aimag(cw).eq.0.0d0)) then
                 CALL cgsolve_all (h_psi_all, cg_psi, et(1,ikk), dvpsi, dpsi(:,:,1), h_diag, &
                       npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol)
        else
                 CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsipm(:,:,1), h_diag, &
                      npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),          &
                      npol, cw,maxter_coul, .true.)
                 CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsipm(:,:,2), h_diag, &
                       npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol,   &
                      -cw,maxter_coul, .true.)
                 dpsi(:,:,iw) = dcmplx(0.5d0,0.0d0)*(dpsipm(:,:,1) + dpsipm(:,:,2))
        endif
     enddo
   else 
!       call cbcg_solve_coul(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsi, dpsic(1,1,1), h_diag, &
!                            npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), &
!                            npol, niters, alphabeta, .false.)
!       dpsi = dpsi^{+}
!       dpsi(:,:,:)    =  dcmplx(0.d0, 0.d0)
!       call coul_multishift(npwx, npwq, nfs, niters, dpsit, dpsic, alphabeta, fiu)
!       dpsi(:,:,:)    = dpsit(:,:,:)
!       call zcopy(npwx*nbnd*nfs, dpsit(1,1,1), 1, dpsi(1,1,1), 1)
!       dpsi = dpsi^{+} + dpsi^{-}
!       dpsit(:,:,:) = dcmplx(0.0d0, 0.0d0)
!       call coul_multishift(npwx, npwq, nfs, niters, dpsit, dpsic, alphabeta, ((-1.0d0,0.0d0)*fiu))
!       dpsi(:,:,:) = dcmplx(0.5d0,0.0d0)*(dpsi(:,:,:) + dpsit(:,:,:))
!       call zaxpy (npwx*npol*nbnd*nfs, dcmplx(0.5d0,0.0d0), dpsit(1,1,1), 1, dpsi(1,1,1), 1)
        dpsipm(:,:,:) = dcmplx(0.0d0,0.0d0)
        call coul_multi(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsipm, h_diag, fiu(1), 2*nfs-1, &
                        npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),          &
                        npol, niters, .false.)
        if(.not.conv_root)   WRITE(1000+mpime, '(5x,"kpoint NC increase maxiter!", i4)') ik
        dpsi(:,:,:)    = dcmplx(0.d0, 0.d0)
        dpsi(:,:,1)    = dpsipm(:,:,1)
        do iw = 2, nfs
           dpsi(:,:,iw) = dcmplx(0.5d0,0.0d0)*(dpsipm(:,:,iw) + dpsipm(:,:,iw+nfs-1))
        enddo
        do ibnd=1, nbnd 
           if (niters(ibnd).ge.maxter_coul) then
              !WRITE(1000+mpime, '(5x,"kpoint NC", i4)') ik
               dpsi(:,ibnd,:) = dcmplx(0.0d0,0.0d0)
           endif
        enddo
   endif
      ltaver = ltaver + lter
      lintercall = lintercall + 1
      nrec1 =  ik
      weight = wk (ikk)
      do iw = 1 , nfs
         call incdrhoscf_w (drhoscf(1, iw) , weight, ik, &
                            dbecsum(1,1,current_spin), dpsi(:,:,iw))
      enddo
     enddo !kpoints
     call mp_sum ( drhoscf(:,:), inter_pool_comm )
     do iw = 1, nfs
        call zcopy (dfftp%nnr*nspin_mag, drhoscf(1,iw),1, dvscfout(1,iw),1)
     enddo
!!!!!!!NEED THIS FOR ULTRASOFT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do iw = 1, nfs
!       call zcopy (nspin_mag*dfftp%nnr, drhoscf(1,iw), 1, drhoscfh(1,iw), 1)
!       call addusddens (drhoscfh(1,iw), dbecsum, imode0, npe, 0)
!       call zcopy (dfftp%nnr*nspin_mag, drhoscfh(1,iw),1, dvscfout(1,iw),1)
!    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SGW: here we enforce zero average variation of the charge density
! if the bare perturbation does not have a constant term
! (otherwise the numerical error, coupled with a small denominator
! in the Coulomb term, gives rise to a spurious dvscf response)
! One wing of the dielectric matrix is particularly badly behaved 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    meandvb = sqrt ((sum(dreal(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0) / float(dffts%nnr)
    do iw = 1, nfs 
       if (meandvb.lt.1.d-8) then 
          CALL fwfft  ('Dense', dvscfout(:,iw), dfftp)
          dvscfout  (nl(1), iw) = dcmplx(0.d0, 0.0d0)
          CALL invfft ('Dense', dvscfout(:,iw), dfftp)
       endif
    enddo
    !
    !   After the loop over the perturbations we have the linear change
    !   in the charge density for each mode of this representation.
    !   Here we symmetrize them ...
    !
    do iw = 1, nfs
       call dv_of_drho (1, dvscfout(1,iw), .true.)
    enddo
    averlt = DBLE (ltaver) / lintercall
    tcpu = get_clock ('GW')
    FLUSH( stdout )
    rec_code=10
  enddo !loop on kter (iterations)
!   after this point drhoscf is dv_hartree(RPA)
!   drhoscf(:,1) = dvscout(:,1)
!  -vc*\Chi
  drhoscf(:,:) = -dvscfout(:,:)
  if (convt) then
   if (fildvscf.ne.' ') then
    write(6, '("fildvscf")') 
   end if
  endif

!  if(.not.prec_direct) deallocate(dpsic)
  deallocate (h_diag)
  deallocate (dvscfout)
  deallocate (dbecsum)

  CALL stop_clock (time_solver)

END SUBROUTINE solve_lindir

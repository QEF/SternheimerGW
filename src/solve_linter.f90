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
SUBROUTINE solve_linter(dvbarein, iw, drhoscf)
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
! grep @10TION for all the instances where I've turned off the regular
! PH parallelism which always assumes we have split k-points across processors.
! Each processor has a full slice of k-points. 
!------------------------------------------------------------------------------


  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin
  USE cell_base,            ONLY : tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : lgauss, degauss, ngauss, xk, wk, nkstot
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,             ONLY : domag
  USE wvfct,                ONLY : nbnd, npw, npwx, igk,g2kin,  et
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE paw_variables,        ONLY : okpaw
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, flmixdpot, current_iq, &
                                   ext_recover, eta, high_io, maxter_coul
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iunrec, iudvscf, iudwfm, iudwfp 
  USE output_mod,           ONLY : fildrho, fildvscf
  USE eqv,                  ONLY : dvpsi, dpsi, evq, eprec
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE lr_symm_base,         ONLY : minus_q, irgq, nsymq, rtau 
  USE freq_gw,              ONLY : fpol, fiu, nfs, nfsmax
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE buffers,              ONLY : save_buffer, get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, ntask_groups, me_bgrp
  USE mp_world,             ONLY : mpime
  USE mp_images, ONLY : intra_image_comm
  USE gvect,           ONLY : ngm, g, nl
  USE gvecs,           ONLY : nls, doublegrid
  USE fft_base,        ONLY : dfftp, dffts
  USE fft_interfaces,  ONLY : invfft, fwfft
#ifdef __NAG
  USE f90_unix_io,     ONLY : flush
#endif

  implicit none

  complex(DP) :: drhoscf (dfftp%nnr, nspin_mag)
  ! output: the change of the scf charge
  complex(DP) :: dvbarein (dffts%nnr)
  complex(DP), allocatable, target :: dvscfin(:,:)
  ! change of the scf potential 
  complex(DP), pointer :: dvscfins (:,:)
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:), dvscfout (:,:)
  complex(DP) :: dpsip(npwx*npol, nbnd), dpsim(npwx*npol, nbnd)
  complex(DP), allocatable :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
                              dbecsum (:,:,:), dbecsum_nc(:,:,:,:), aux1 (:,:)
  complex(DP) :: cw
  complex(DP), allocatable :: etc(:,:)
  complex(kind=DP) :: ZDOTC
  ! Misc variables for metals
  ! dos_ef: density of states at Ef
  real(DP), external :: w0gauss, wgauss
  real(DP) , allocatable :: h_diag (:,:)
  real(DP) :: thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight, aux_avg (2)
  !HL dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  ! becsum1 PAW array.
  real(DP), allocatable :: becsum1(:,:,:)
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
  integer :: iw, ir 
  integer :: irr, imode0, npe

  external ch_psi_all, cg_psi, h_psi_all, ch_psi_all_nopv
  real(kind=DP)    :: DZNRM2
  external ZDOTC, DZNRM2

  logical :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0,     & ! true if xq=(0,0,0) in a metal
             cgsolver          

  allocate (dpsi(npwx*npol, nbnd))
 
  IF (rec_code_read > 20 ) RETURN

!HL- Allocate arrays for dV_scf (need to alter these from (dfftp%nnr, nspin_mag, npe) 
!to just (dfftp%nnr, nspin_mag).
  npe    = 1
  imode0 = 1
  irr    = 1
  ipert  = 1
  lter   = 0
  lmres  = 1

  call start_clock ('solve_linter')

  allocate (dvscfout ( dfftp%nnr , nspin_mag))    
  allocate (drhoscfh ( dfftp%nnr , nspin_mag))    
  allocate (etc(nbnd, nkstot))
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat, nspin_mag))

  allocate (dvscfin  ( dfftp%nnr , nspin_mag))    
  if (doublegrid) then
     allocate (dvscfins ( dffts%nnr , nspin_mag))    
  else
     dvscfins => dvscfin
  endif

!Complex eigenvalues
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin))
  allocate (aux1 ( dffts%nnr, npol))    
  allocate (h_diag ( npwx*npol, nbnd))    

  iter0 = 0
  convt =.FALSE.
  where_rec='no_recover'

! The outside loop is over the iterations.
! niter_gw := maximum number of iterations
  
  do kter = 1, niter_gw
     iter = kter + iter0
     ltaver = 0
     lintercall = 0
     drhoscf(:,:)   = (0.d0, 0.d0)
     drhoscfh(:,:)  = (0.d0, 0.d0)
     dbecsum(:,:,:) = (0.d0, 0.d0)

     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     if (nksq.gt.1) rewind (unit = iunigk)
!start kpoints loop
     do ik = 1, nksq
        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_linter', 'reading igk', abs (ios) )
        endif
! lgamma is a q=0 computation
        if (lgamma)  npwq = npw
! k and k+q mesh defined in initialize_gw:
!       ikks(ik) = 2 * ik - 1
!       ikqs(ik) = 2 * ik
        ikk = ikks(ik)
        ikq = ikqs(ik)

        if (lsda) current_spin = isk (ikk)
        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('solve_linter', 'reading igkq', abs (ios) )
        endif
       !Calculates beta functions (Kleinman-Bylander projectors), with
       !structure factor, for all atoms, in reciprocal space
       !HL the beta functions (vkb) are being generated properly.  
        call init_us_2 (npwq, igkq, xk (1, ikq), vkb)

       !Reads unperturbed wavefuctions psi(k) and psi(k+q)

        if (nksq.gt.1) then
           if (lgamma) then
              call get_buffer (evc, lrwfc, iuwfc, ikk)
           else
              call get_buffer (evc, lrwfc, iuwfc, ikk)
              call get_buffer (evq, lrwfc, iuwfc, ikq)
           endif
        endif

!IS TPA preconditioner better?
        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                          (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                          (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
        enddo
        h_diag = 0.d0
        do ibnd = 1, nbnd_occ (ikk)
           do ig = 1, npwq
              h_diag(ig,ibnd)= 1.d0/max(1.0d0, g2kin(ig)/eprec(ibnd,ik))
           enddo
           IF (noncolin) THEN
              do ig = 1, npwq
                 h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
              enddo
           END IF
        enddo
!HL indices freezing perturbations.
        nrec = ik
!and now adds the contribution of the self consistent term
        if (where_rec =='solve_lint'.or.iter>1) then
          ! After the first iteration dvbare_q*psi_kpoint is read from file
           call get_buffer (dvpsi, lrbar, iubar, nrec)

           call start_clock ('vpsifft')
           do ibnd = 1, nbnd_occ (ikk)
              call cft_wave (ik, evc (1, ibnd), aux1, +1)
              call apply_dpot(dffts%nnr, aux1, dvscfins(1,1), current_spin)
              call cft_wave (ik, dvpsi (1, ibnd), aux1, -1)
           enddo
           call stop_clock ('vpsifft')
           !  In the case of US pseudopotentials there is an additional
           !  selfconsist term which comes from the dependence of D on
           !  V_{eff} on the bare change of the potential
           !
           !Need to check this for ultrasoft              
           !HL THIS TERM PROBABLY NEEDS TO BE INCLUDED.
           !KC: This term needs to be included for USPP.
           !KC: add the augmentation charge term for dvscf
         !!  call adddvscf (1, ik)
        else
            call dvqpsi_us (dvbarein, ik, .false.)
         ! USPP
         ! add the augmentation charge term for dvext and dbext
         ! call adddvscf (1, ik)
           call save_buffer (dvpsi, lrbar, iubar, nrec)
        endif
        ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
        ! Apply -P_c^+.
        !-P_c^ = - (1-P_v^):
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .false.)
        
        if(where_rec=='solve_lint'.or.iter > 1) then
           if(high_io) then
              call get_buffer( dpsip, lrdwf, iudwfp, ik)
              if(iw.gt.1) call get_buffer( dpsim, lrdwf, iudwfm, ik)
           else
              dpsim(:,:)     = (0.d0, 0.d0) 
              dpsip(:,:)     = (0.d0, 0.d0) 
           endif
           thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
        else
         !
         ! At the first iteration dpsi and dvscfin are set to zero
         !
          dpsi(:,:)      = (0.d0, 0.d0) 
          dpsim(:,:)     = (0.d0, 0.d0) 
          dpsip(:,:)     = (0.d0, 0.d0) 
          dvscfin(:, :)  = (0.d0, 0.d0)
          dvscfout(:, :) = (0.d0, 0.d0)
          thresh         =  1.0d-2
        endif

        conv_root = .true.
        etc(:,:)  = CMPLX( et(:,:), 0.0d0 , kind=DP)
        cw        = fiu(iw) 

        if (iw.eq.1) then
               call cgsolve_all (h_psi_all, cg_psi, et(1,ikk), dvpsi, dpsip, h_diag, & 
                      npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), &
                      npol)
               do ibnd = 1, nbnd_occ(ikk)
                  call ZCOPY (npwx*npol, dpsip (1, ibnd), 1, dpsi(1, ibnd), 1)
               enddo
        else 

              conv_root = .true.
              call cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsip, h_diag, &
                   npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk),   &
                   npol, cw, maxter_coul, .true.)

              conv_root = .true.
              call cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsim, h_diag,     &
                   npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, &
                  -cw, maxter_coul, .true.)

              dpsi(:,:) = dcmplx(0.0d0, 0.0d0)
              do ibnd =1 , nbnd_occ(ikk)
                 call ZAXPY (npwx*npol, dcmplx(0.5d0,0.0), dpsim(1,ibnd), 1, dpsi(1,ibnd), 1)
                 call ZAXPY (npwx*npol, dcmplx(0.5d0,0.0), dpsip(1,ibnd), 1, dpsi(1,ibnd), 1)
              enddo
        endif

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

     meandvb = sqrt ((sum(dreal(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0 )/float(dffts%nnr)
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
        call mix_potential_c(dfftp%nnr*nspin_mag, dvscfout(1,1), dvscfin(1,1), alpha_mix(kter),& 
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
              call cinterpolate (dvscfin(1,is), dvscfins(1,is), -1)
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
     CALL FLUSH( stdout )
     rec_code=10

!     WRITE(1000+mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!           "secs   av.it.: ",f5.1)') iter, tcpu, averlt
!     WRITE(1000+mpime, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
!           &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
     if (convt) goto 155

  enddo !loop on kter (iterations)

155 iter0=0

   WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
        & "secs av.it.:",f5.1)') iter, tcpu, averlt
!   WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
!          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
!   WRITE(1000+mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!        "secs   av.it.: ",f5.1)') iter, tcpu, averlt

  drhoscf(:,:) = dvscfin(:,:)
  
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
  call stop_clock ('solve_linter')
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
  RETURN
  !
END SUBROUTINE

!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! grep @10TION for all the instances where I've turned off the regular
! PH parallelism which always assumes we have split k-points across processors.
! Each processor has a full slice of k-points. 


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
                                   ext_recover, eta
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iunrec, iudvscf, iudwfm, iudwfp 
  USE output,               ONLY : fildrho, fildvscf
  USE eqv,                  ONLY : dvpsi, dpsi, evq, eprec
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE modes,                ONLY : npertx, npert, u, t, irotmq, tmq, &
                                   minus_q, irgq, nsymq, rtau 
  USE freq_gw,              ONLY : fpol, fiu, nfs, nfsmax
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE buffers,              ONLY : save_buffer, get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, ntask_groups, me_bgrp
  USE mp_world,             ONLY : mpime
  USE gvect,           ONLY : ngm, g, nl
  USE gvecs,           ONLY : nls, doublegrid
  USE fft_base,        ONLY : dfftp, dffts
  USE fft_interfaces,  ONLY : invfft, fwfft

  implicit none

  ! counter on frequencies.

  integer :: iw, ir 
  integer :: irr, imode0, npe

  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  complex(DP) :: drhoscf (dfftp%nnr, nfs)
  ! output: the change of the scf charge
  complex(DP) :: dvbarein (dffts%nnr)

! HL prec
! HL careful now... complexifying preconditioner:
  real(DP) , allocatable :: h_diag (:,:)
! h_diag: diagonal part of the Hamiltonian
  real(DP) :: thresh, anorm, averlt, dr2

  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight, aux_avg (2)

  ! Misc variables for metals
  ! dos_ef: density of states at Ef
  real(DP), external :: w0gauss, wgauss

  ! functions computing the delta and theta function
  complex(DP), allocatable, target :: dvscfin(:,:)

  ! change of the scf potential 
  complex(DP), pointer :: dvscfins (:,:)

  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:), dvscfout (:,:)

  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  complex(DP), allocatable :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
      dbecsum (:,:,:), dbecsum_nc(:,:,:,:), aux1 (:,:)
  complex(DP) :: cw
  complex(DP), allocatable :: etc(:,:)


  !HL dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  ! becsum1 PAW array.
  REAL(DP), allocatable :: becsum1(:,:,:)

  !For approx, mixing scheme.
  real(kind=DP) :: DZNRM2
  complex(kind=DP) :: ZDOTC
  external ZDOTC, DZNRM2

  logical :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0,     & ! true if xq=(0,0,0) in a metal
             cgsolver          
  
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

  real(DP) :: tcpu, get_clock ! timing variables
  real(DP) :: meandvb
 
  external ch_psi_all, cg_psi, h_psi_all
  COMPLEX(DP) :: dpsip(npwx*npol, nbnd), dpsim(npwx*npol, nbnd)

  allocate (dpsi(npwx*npol, nbnd))
 
  IF (rec_code_read > 20 ) RETURN

  !HL- Allocate arrays for dV_scf (need to alter these from (dfftp%nnr, nspin_mag, npe) to just (dfftp%nnr, nspin_mag).
  npe    = 1
  imode0 = 1
  irr    = 1
  ipert  = 1
  lter   = 0
  lmres  = 1

!HLallocate (hpsi(npwx*npol, 4)) ! Test array for whether linear system is being properly solved
  call start_clock ('solve_linter')

  allocate (dvscfout ( dfftp%nnr , nfs))    
  allocate (drhoscfh ( dfftp%nnr , nfs))    
  allocate (dvscfin  ( dfftp%nnr , nfs))    
  allocate (etc(nbnd, nkstot))
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat, nspin_mag))    

  if (doublegrid) then
     allocate (dvscfins ( dffts%nnr , nfs))    
  else
     dvscfins => dvscfin
  endif

!Complex eigenvalues
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat , nspin))
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
              h_diag(ig,ibnd)= 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
           enddo
           IF (noncolin) THEN
              do ig = 1, npwq
                 h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
              enddo
           END IF
        enddo
!HL indices freezing perturbations.
           mode = 1
           nrec = ik
!and now adds the contribution of the self consistent term
           if (where_rec =='solve_lint'.or.iter>1) then
             ! After the first iteration dvbare_q*psi_kpoint is read from file
              call get_buffer (dvpsi, lrbar, iubar, nrec)
             ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
             ! dvscf_q from previous iteration (mix_potential)
              call start_clock ('vpsifft')
              do ibnd = 1, nbnd_occ (ikk)
                 call cft_wave (evc (1, ibnd), aux1, +1)
                 call apply_dpot(aux1, dvscfins(1,iw), current_spin)
                 call cft_wave (dvpsi (1, ibnd), aux1, -1)
              enddo
              call stop_clock ('vpsifft')
           else
               call dvqpsi_us (dvbarein, ik, .false.)
               call save_buffer (dvpsi, lrbar, iubar, nrec)
           endif
        ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
        ! Apply -P_c^+.
        ! -P_c^ = - (1-P_v^):
        !
           CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi)
        !
           if (where_rec=='solve_lint'.or.iter > 1) then
              call get_buffer( dpsip, lrdwf, iudwfp, ik)
              call get_buffer( dpsim, lrdwf, iudwfm, ik)

             !  dpsi(:,:)  = (0.d0, 0.d0) 
             !  dpsim(:,:) = (0.d0, 0.d0) 
             !  dpsip(:,:) = (0.d0, 0.d0) 
             !threshold for iterative solution of the linear system
             !write(6,*)1.d-1*sqrt(dr2), 1.d-4
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
              !
              ! starting threshold for iterative solution of the linear system
              !
              thresh = 1.0d-2
           endif

       etc(:,:) = CMPLX( et(:,:), 0.0d0 , kind=DP)
       cw       = fiu(iw) 

!HL should just use cgsolve_all when fiu(iw) = 0.0d0! 
!probably gain at least factor of two for static case.

       conv_root = .true.

       IF (iw.eq.1) THEN
               CALL cgsolve_all (h_psi_all, cg_psi, et(1,ikk), dvpsi, dpsip, h_diag, & 
                      npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol)
               dpsim(:,:) = dpsip(:,:)
               dpsi(:,:) = dcmplx(0.5d0,0.0d0)*(dpsim(:,:) + dpsip(:,:) ) 
       ELSE
               CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsip, h_diag, &
                     npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, cw, .true.)
               CALL cbcg_solve(ch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsim, h_diag, &
                     npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, -cw, .true.)
               dpsi(:,:) = dcmplx(0.5d0,0.0d0)*(dpsim(:,:) + dpsip(:,:) ) 
       ENDIF

       ltaver = ltaver + lter
       lintercall = lintercall + 1

       IF (.NOT.conv_root) WRITE(1000+mpime, '(5x,"kpoint",i4," ibnd",i4, &
                  &              " solve_linter: root not converged ",e10.3)')  &
                  &                ik , ibnd, anorm
           nrec1 =  ik
         ! calculates dvscf, sum over k => dvscf_q_ipert
         ! incdrhoscf:  This routine computes the change of the charge density due to the
           call save_buffer (dpsim, lrdwf, iudwfp, ik)
           call save_buffer (dpsip, lrdwf, iudwfm, ik)
         ! perturbation. It is called at the end of the computation of the
         ! change of the wavefunction for a given k point.
           weight = wk (ikk)
           call incdrhoscf_w (drhoscf(1, iw) , weight, ik, &
                              dbecsum(1,1,current_spin), dpsi(1,1))
     enddo 

!HLPARA
!     call mp_sum ( dbecsum, intra_bgrp_comm )

     if (doublegrid) then
         do is = 1, nspin_mag
            call cinterpolate (drhoscfh(1,iw), drhoscf(1,iw), 1)
         enddo
     else
            call zcopy (nspin_mag*dfftp%nnr, drhoscf(1,iw), 1, drhoscfh(1,iw), 1)
     endif

!     call mp_sum ( dvscfout, inter_pool_comm )
!    if (.not.lgamma_gamma) then
!         call psyme (dvscfout)
!         IF ( noncolin.and.domag ) CALL psym_dmage(dvscfout)
!    endif
     call zcopy (dfftp%nnr*nspin_mag, drhoscfh(1,iw), 1, dvscfout(1,iw),1)

     meandvb = sqrt ((sum(dreal(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0 )/float(dffts%nnr)

     if (meandvb.lt.1.d-8) then 
         CALL fwfft ('Dense', dvscfout(:,iw), dfftp)
         dvscfout ( nl(1), current_spin ) = (0.d0, 0.0d0)
         CALL invfft ('Dense', dvscfout(:,iw), dfftp)
     endif

     call mp_sum ( drhoscf, inter_pool_comm )
     call mp_sum ( drhoscfh, inter_pool_comm )
!     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )


! for q->0 the Fermi level can shift.
! IF (lmetq0) call ef_shift(drhoscfh,ldos,ldoss,dos_ef,irr,npe,.false.)

     call dv_of_drho (1, dvscfout(1,iw), .true.)

     nmix_gw = 5

     if (iw.eq.1) then
!Density reponse in real space should be real at zero freq no matter what!
!just using standard broyden for the zero freq. case.
        call mix_potential(2*dfftp%nnr*nspin_mag, dvscfout(1,iw), dvscfin(1,iw), alpha_mix(kter), &
                           dr2, tr2_gw, iter, nmix_gw, flmixdpot, convt)
     else
        call mix_potential_c(dfftp%nnr, dvscfout(1,iw), dvscfin(1,iw), &
                             alpha_mix(kter), dr2, tr2_gw, iter, &
                             nmix_gw, convt)
     endif

     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin_mag
              call cinterpolate (dvscfin(1,iw), dvscfins(1,iw), -1)
           enddo
        enddo
     endif

#ifdef __PARA
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif
     tcpu = get_clock ('SGW')
     dr2 = dr2 
     CALL flush_unit( stdout )
     rec_code=10
     if (convt) goto 155
  enddo !loop on kter (iterations)

155 iter0=0

   WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
         "secs av.it.:",f5.1)') iter, tcpu, averlt
   WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
!   WRITE(1000+mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!        "secs   av.it.: ",f5.1)') iter, tcpu, averlt

  drhoscf(:,iw) = dvscfin(:,iw)

  if (convt) then
   if (fildvscf.ne.' ') then
       write(6, '("fildvscf")') 
   end if
  endif
 
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

SUBROUTINE setmixout(in1, in2, mix, dvscfout, dbecsum, ndim, flag )
USE kinds, ONLY : DP
USE mp_global, ONLY : intra_pool_comm
USE mp, ONLY : mp_sum
IMPLICIT NONE
INTEGER :: in1, in2, flag, ndim, startb, lastb
COMPLEX(DP) :: mix(in1+in2), dvscfout(in1), dbecsum(in2)

CALL divide (in2, startb, lastb)
ndim=lastb-startb+1

IF (flag==-1) THEN
   mix(1:in1)=dvscfout(1:in1)
   mix(in1+1:in1+ndim)=dbecsum(startb:lastb)
ELSE
   dvscfout(1:in1)=mix(1:in1)
   dbecsum=(0.0_DP,0.0_DP)
   dbecsum(startb:lastb)=mix(in1+1:in1+ndim)
#ifdef __PARA
   CALL mp_sum(dbecsum, intra_pool_comm)
#endif
ENDIF
END SUBROUTINE setmixout


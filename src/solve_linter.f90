!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
!SUBROUTINE solve_linter(igpert, iw, drhoscf)
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
  USE gvect,                ONLY : nrxx, g, nl, nr1, nr2, nr3, nrx1, nrx2, nrx3
  USE gsmooth,              ONLY : doublegrid, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,             ONLY : domag
  USE wvfct,                ONLY : nbnd, npw, npwx, igk,g2kin,  et
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential, paw_dusymmetrize, &
                                   paw_dumqsymmetrize
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, flmixdpot, current_iq, &
                                   ext_recover, eta
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iunrec, iudvscf, iudwfm, iudwfp 
  USE output,               ONLY : fildrho, fildvscf
  USE gwus,                 ONLY : int3_paw, becsumort
  USE eqv,                  ONLY : dvpsi, dpsi, evq, eprec, dpsim, dpsip
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE modes,                ONLY : npertx, npert, u, t, irotmq, tmq, &
                                   minus_q, irgq, nsymq, rtau 
  USE recover_mod,          ONLY : read_rec, write_rec

  ! used oly to write the restart file
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm, mpime, mp_global_end,&
                                   intra_image_comm

  USE mp,                   ONLY : mp_sum, mp_barrier
  !  
  USE freq_gw,              ONLY : fpol, fiu, nfs, nfsmax

  implicit none

  ! counter on frequencies.

  integer :: iw 
  integer :: irr, imode0, npe

  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  complex(DP) :: drhoscf (nrxx, nspin_mag)
  ! output: the change of the scf charge
  complex(DP) :: dvbarein (nrxxs)

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
      dbecsum (:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:)
  complex(DP) :: cw
  complex(DP), allocatable :: etc(:,:)


  !HL dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  ! becsum1 PAW array.
  REAL(DP), allocatable :: becsum1(:,:,:)
  COMPLEX (DP), ALLOCATABLE :: hpsi(:,:)
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
             igpert        ! bare perturbation g vector.

  real(DP) :: tcpu, get_clock ! timing variables
  real(DP) :: meandvb
 
  !external ch_psi_all, cg_psi, ccg_psi, cch_psi_all_fix
   external ch_psi_all, cg_psi, cch_psi_all_fix
  
  IF (rec_code_read > 20 ) RETURN

  !HL- Allocate arrays for dV_scf (need to alter these from (nrxx, nspin_mag, npe) to just (nrxx, nspin_mag).
  npe    = 1
  imode0 = 1
  irr    = 1
  ipert  = 1
  lter   = 0

!HLallocate (hpsi(npwx*npol, 4)) ! Test array for whether linear system is being properly solved
  call start_clock ('solve_linter')

  !HL allocate (dvscfin ( nrxx , nspin_mag , npe))
  allocate (dvscfout ( nrxx , nspin_mag))    
  allocate (dvscfin ( nrxx , nspin_mag))    

  if (doublegrid) then
  !HL allocate (dvscfins ( nrxxs , nspin_mag , npe))
     allocate (dvscfins ( nrxxs , nspin_mag))    
  else
     dvscfins => dvscfin
  endif

  allocate (drhoscfh ( nrxx , nspin_mag))    
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat, nspin_mag))    

!Complex eigenvalues
  allocate (etc(nbnd, nkstot))

  IF (okpaw) THEN
     allocate (mixin(nrxx*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     allocate (mixout(nrxx*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     mixin=(0.0_DP,0.0_DP)
  ENDIF
  IF (noncolin) allocate (dbecsum_nc (nhm,nhm, nat , nspin, npe))
  allocate (aux1 ( nrxxs, npol))    
  allocate (h_diag ( npwx*npol, nbnd))    

  if (rec_code_read == 10.AND.ext_recover) then
    ! restart from GW calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh, dbecsum)
        CALL setmixout(npe*nrxx*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*npe)/2, &
                    mixin, dvscfin, dbecsum, ndim, -1 )
     ELSE
     ENDIF
     rec_code=0
  else
    iter0 = 0
    convt =.FALSE.
    where_rec='no_recover'
  endif

  IF (convt) GOTO 155
! In this case it has recovered after computing the contribution
! to the dynamical matrix. This is a new iteration that has to 
! start from the beginning.

  IF (iter0==-1000) iter0=0

! The outside loop is over the iterations.
! niter_gw := maximum number of iterations

!HL all fine...
!write(6,*) nbnd_occ(:)

  do kter = 1, niter_gw
     iter = kter + iter0
     ltaver = 0

     lintercall = 0
     drhoscf(:,:) = (0.d0, 0.d0)
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
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
           else
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
              call davcio (evq, lrwfc, iuwfc, ikq, - 1)
           endif
        endif

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
             !HL is this necessary?
             ! After the first iteration dvbare_q*psi_kpoint is read from file
              call davcio (dvpsi, lrbar, iubar, nrec, - 1)
             ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
             ! dvscf_q from previous iteration (mix_potential)
              call start_clock ('vpsifft')
              do ibnd = 1, nbnd_occ (ikk)
                !FFT translated according to igk
                 call cft_wave (evc (1, ibnd), aux1, +1) 
                 call apply_dpot(aux1, dvscfins(1,1), current_spin)
                !FFT translated according to igkq: DeltaV(q)psi(k).
                 call cft_wave (dvpsi (1, ibnd), aux1, -1)
              enddo
              call stop_clock ('vpsifft')
           else
        !  At the first iteration dvbare_q*psi_kpoint is calculated
        !  and written to file
            call dvqpsi_us (dvbarein, ik, 1, .false.)
            call davcio (dvpsi, lrbar, iubar, nrec, +1)
           endif

        ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
        ! Apply -P_c^+.
        ! -P_c^ = - (1-P_v^):
        ! SGW := call emptyproj ( evq, dvpsi)
           CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi)
           if (where_rec=='solve_lint'.or.iter > 1) then
              !starting value for delta_psi is read from iudwf
               nrec1 = ik
              !HL Don't need to read/write the full wave fxn at each iteration...
             ! call davcio ( dpsi, lrdwf, iudwf, nrec1, -1)
             ! call davcio ( dpsip, lrdwf, iudwfp, nrec1, -1)
             ! call davcio ( dpsim, lrdwf, iudwfm, nrec1, -1)
               dpsi(:,:)  = (0.d0, 0.d0) 
               dpsim(:,:) = (0.d0, 0.d0) 
               dpsip(:,:) = (0.d0, 0.d0) 
              !threshold for iterative solution of the linear system
              !write(6,*)1.d-1*sqrt(dr2), 1.d-4
               thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
           else
            !
            ! At the first iteration dpsi and dvscfin are set to zero
            !
              dpsi(:,:) = (0.d0, 0.d0) 
              dpsim(:,:) = (0.d0, 0.d0) 
              dpsip(:,:) = (0.d0, 0.d0) 
              dvscfin(:, :) = (0.d0, 0.d0)
              !
              ! starting threshold for iterative solution of the linear system
              !
              thresh = 1.0d-2
           endif

       etc(:,:) = CMPLX( et(:,:), 0.0d0 , kind=DP)
       cw       = CMPLX(0.0d0,  fiu(iw), kind=DP) 

!HL should just use cgsolve_all when fiu(iw) = 0.0d0! probably gain at least factor of 
!two for static case...

       if(iw.eq.1) then
           call cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, dpsip, h_diag, & 
                      npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol)
                      dpsim(:,:) = dpsip(:,:)
       else
            call cbcg_solve_fix(cch_psi_all_fix, cg_psi, etc(1,ikk), dvpsi, dpsip, h_diag, &
                       npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, cw, .true.)

            call cbcg_solve_fix(cch_psi_all_fix, cg_psi, etc(1,ikk), dvpsi, dpsim, h_diag, &
                      npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, -cw, .true.)
       endif

           ltaver = ltaver + lter
           lintercall = lintercall + 1

         if (.not.conv_root) WRITE(1000+mpime, '(5x,"kpoint",i4," ibnd",i4, &
              &              " solve_linter: root not converged ",e10.3)')  &
              &                ik , ibnd, anorm

           nrec1 =  ik
           dpsi(:,:) = (0.5d0,0.0d0) * (dpsim(:,:) + dpsip(:,:) ) 

          ! call davcio (dpsi, lrdwf, iudwf, nrec1, + 1)
          ! call davcio (dpsim, lrdwf, iudwfm, nrec1, + 1)
          ! call davcio (dpsip, lrdwf, iudwfp, nrec1, + 1)

          ! calculates dvscf, sum over k => dvscf_q_ipert
          ! incdrhoscf:  This routine computes the change of the charge density due to the
          ! perturbation. It is called at the end of the computation of the
          ! change of the wavefunction for a given k point.
           weight = wk (ikk)
           IF (noncolin) THEN
              call incdrhoscf_nc(drhoscf(1,1),weight,ik, &
                                       dbecsum_nc(1,1,1,1,ipert))
           ELSE
              call incdrhoscf ( drhoscf(1,current_spin) , weight, ik, &
                                dbecsum(1,1,current_spin))
           END IF
     enddo 

        if (doublegrid) then
             do is = 1, nspin_mag
                call cinterpolate (drhoscfh(1,is), drhoscf(1,is), 1)
             enddo
        else
                call zcopy (nspin_mag*nrxx, drhoscf, 1, drhoscfh, 1)
        endif

        call addusddens (drhoscfh, dbecsum, imode0, npe, 0)

!!!! HL @10TION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef __PARA
!   Reduce the delta rho across pools
!    call mp_sum ( drhoscf, inter_pool_comm )
!    call mp_sum ( drhoscfh, inter_pool_comm )
!    IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
!#endif

!       if (fildrho.ne.' ') call davcio_drho (drhoscfh(1,1), lrdrho, &
!                                             iudrho, imode0+ipert, +1)

        call zcopy (nrxx*nspin_mag,drhoscfh(1,1),1,dvscfout(1,1),1)

     ! SGW: here we enforce zero average variation of the charge density
     ! if the bare perturbation does not have a constant term
     ! (otherwise the numerical error, coupled with a small denominator
     ! in the coulomb term, gives rise to a spurious dvscf response)
     ! One wing of the dielectric matrix is particularly badly behaved 

          meandvb = sqrt ( (sum(dreal(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0 ) / float(nrxxs)
          if (meandvb.lt.1.d-8) then 
     !        WRITE(1000+mpime,'("Zeroing dvscf")')  
     !        WRITE(1000+mpime,*) meandvb
     !        WRITE(1000+mpime,*) dvscfout(nl(1), current_spin) 
             call cft3 (dvscfout, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
             dvscfout ( nl(1),current_spin ) = (0.d0, 0.0d0)
             call cft3 (dvscfout, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     !     else
     !        WRITE(1000+mpime,'("NOT Zeroing dvscf")')  
     !        WRITE(1000+mpime,*) meandvb
     !        WRITE(1000+mpime,*) dvscfout(nl(1), current_spin) 
          endif
        call dv_of_drho (1, dvscfout(1,1), .true.)

        call mix_potential_c(nrxx, dvscfout, dvscfin, &
                             alpha_mix(kter), dr2, tr2_gw, iter, &
                             nmix_gw, convt)
     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin_mag
              call cinterpolate (dvscfin(1,is), dvscfins(1,is), -1)
           enddo
        enddo
     endif

     ! with the new change of the potential we compute the integrals
     ! of the change of potential and Q 
     ! HL-Q denotes the augmentation charge. Look at Vanderbilt PRB 41 7892
     ! Q_{ij} = <psi_{j}|psi_{i}> - <phi_{i}|phi_{j}> 
     ! USPP valence charge density is described 
     ! n_v(r) =  \sum_{n,k} \phi^{*}_{nk} (r) \phi_{nk}(r) + \sum_{i,j} p_{i,j}Q_{j,i}
     ! p_{i,j} = \sum_{n,k} <beta_{i}|\phi_{nk}><phi_{nk}|beta_{j}> 

#ifdef __PARA
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)

!@10TION
!     call mp_sum ( aux_avg, inter_pool_comm )

     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif
     tcpu = get_clock ('GW')
     dr2 = dr2 / DBLE(npe)

!   WRITE( 1000 + mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!   " secs   av.it.: ",f5.1)') iter, tcpu, averlt
!   WRITE( 1000+mpime, '(5x," thresh=",e10.3, " alpha_mix = ",f6.3, &
!   " |ddv_scf|^2 = ",e10.3 )') thresh, alpha_mix (kter) , dr2
!   Here we save the information for recovering the run from this point
!   HL-  CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, dvscfin, drhoscfh) 
!   write_rec('solve_linte, k, dr2, iter, cont, q, dvscfin, drhoscfh) 
     
     CALL flush_unit( stdout )

     rec_code=10

!@10TION more hanging, more clumsy parallelisation.
!     IF (okpaw) THEN
!        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
!                                               dvscfin, drhoscfh, dbecsum)
!     ELSE
!        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
!                                               dvscfin, drhoscfh)
!     ENDIF
!     if (check_stop_now()) call stop_smoothly_gw (.false.)

     if (convt) goto 155
  enddo !loop on kter (iterations)

155 iter0=0

!   WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!         "secs av.it.:",f5.1)') iter, tcpu, averlt
!   WRITE(1000+mpime, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!        "secs   av.it.: ",f5.1)') iter, tcpu, averlt

!   HL setting drhoscf to dvscfin here this is a temporary hack. 
!   need to understand why drhoscf is zeroed in PH code...
!   possibly because they write to disc davcio_drho ?
!   drhoscf(igpert, 1) = dvscfin + dvbare
!   after this point drhoscf is dv_hartree(RPA)
    drhoscf(:,1) = dvscfin(:,1)
    if (convt) then
    if (fildvscf.ne.' ') then
    write(6, '("fildvscf")') 
  !HL 
  ! do ipert = 1, npe
  ! call davcio_drho ( dvscfin(1,1),  lrdrho, iudvscf, 
  ! imode0 + ipert+(current_iq-1)*3*nat, +1 )
  ! enddo
     end if
  endif

  deallocate (h_diag)
  deallocate (aux1)
  deallocate (dbecsum)

  IF (okpaw) THEN
     if (lmetq0.and.allocated(becsum1)) deallocate (becsum1)
     deallocate (mixin)
     deallocate (mixout)
  ENDIF
  IF (noncolin) deallocate (dbecsum_nc)
  deallocate (dvscfout)
  deallocate (drhoscfh)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
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


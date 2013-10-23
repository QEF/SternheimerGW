!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
SUBROUTINE solve_direct_shift(dvbarein, iw, drhoscf)
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

  complex(DP) :: drhoscf (nrxx, nspin_mag, nfs)
  ! output: the change of the scf charge
  complex(DP) :: dvbarein (nrxxs)

  complex(DP), allocatable :: dpsiwp(:,:), dpsiwm(:,:)
! HL prec
! HL careful now... complexifying preconditioner:
  real(DP) , allocatable :: h_diag (:,:)
! h_diag: diagonal part of the Hamiltonian
  real(DP) :: thresh, anorm, averlt, dr2
  real(DP) :: x
!keeps track of iterations for seed system
  INTEGER, ALLOCATABLE  :: niters(:)
  INTEGER               :: gveccount , ngvecs

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
  complex(DP), allocatable :: drhoscfh (:,:,:), dvscfout (:,:,:)

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
             igpert,     & ! bare perturbation g vector.
             lmres         ! number of gmres iterations to include when using bicgstabl.

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
  lmres  = 1

  call start_clock ('solve_linter')

  allocate (dvscfout ( nrxx , nspin_mag, nfs))
  allocate (dvscfin  ( nrxx , nspin_mag))


  allocate (dpsiwm  (npwx, nfs))
  allocate (dpsiwp  (npwx, nfs))


  if (doublegrid) then
     allocate (dvscfins ( nrxxs , nspin_mag))    
  else
     dvscfins => dvscfin
  endif

  allocate (drhoscfh ( nrxx, nspin_mag, nfs))

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

  iter0 = 0
  convt =.FALSE.
  where_rec='no_recover'

  IF (convt) GOTO 155
  IF (iter0==-1000) iter0=0

! No self-consistency:
  do kter = 1, 1
     iter = kter + iter0
     ltaver = 0

     lintercall = 0
     drhoscf(:,:,:) = (0.d0, 0.d0)
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

!No preconditioning in multishift.
        h_diag = 0.d0
        do ibnd = 1, nbnd_occ (ikk)
           do ig = 1, npwq
                  h_diag(ig,ibnd) = 1.0d0
           enddo
        enddo

!HL indices freezing perturbations.
           mode = 1
           nrec = ik

!and now adds the contribution of the self consistent term
           if (where_rec =='solve_lint'.or.iter>1) then
             !After the first iteration dvbare_q*psi_kpoint is read from file
              call davcio (dvpsi, lrbar, iubar, nrec, - 1)
             !calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
             !dvscf_q from previous iteration (mix_potential)
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
         ! At the first iteration dvbare_q*psi_kpoint is calculated
         ! and written to file
            call dvqpsi_us (dvbarein, ik, 1, .false.)
         ! call davcio (dvpsi, lrbar, iubar, nrec, +1)
           endif
         ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
         ! Apply -P_c^+.
         ! -P_c^ = - (1-P_v^):
           CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi)
         !
         ! At the first iteration dpsi and dvscfin are set to zero
         !
           dpsi(:,:)     = (0.d0, 0.d0) 
           dpsim(:,:)    = (0.d0, 0.d0) 
           dpsip(:,:)    = (0.d0, 0.d0) 
           dvscfin(:, :) = (0.d0, 0.d0)
         ! starting threshold for iterative solution of the linear system
           thresh        = tr2_gw
           etc(:,:)      = CMPLX(et(:,:), 0.0d0 , kind=DP)

!want to construct full frequency density response, one band at a time, 
!for each k.
        gveccount = 1 
        ngvecs = 1
        if(.not.allocated(niters)) ALLOCATE(niters(ngvecs))
        niters = 0 
        do ibnd = 1, nbnd_occ(ikk)
           call cbcg_solve_coul(cch_psi_all_fix, cg_psi, etc(1,ikk), dvpsi(:,ibnd), dpsip, h_diag,  &
                                  npwx, npwq, thresh, ik, lter, conv_root, anorm, 1, npol, &
                                  cw, niters(gveccount))

           if (.not.conv_root) WRITE(1000+mpime, '(5x,"kpoint",i4," ibnd",i4,    &
                    &               "solve_linter: root not converged ",e10.3)') &
                    &                ik , ibnd, anorm

           call coul_multishift(npwx, npwq, nfs, niters(gveccount), -fiu(iw), 1, dpsiwm(:,:))
           call coul_multishift(npwx, npwq, nfs, niters(gveccount), fiu(iw), 1, dpsiwp(:,:))

           ltaver = ltaver + lter
           lintercall = lintercall + 1

           nrec1 =  ik
           dpsi(:,:) = (0.5d0,0.0d0) * (dpsiwm(:,:) + dpsiwp(:,:)) 

          !calculates dvscf, sum over k => dvscf_q_ipert
          !incdrhoscf:  This routine computes the change of the charge 
          !density due to the perturbation. It is called at the end of 
          !the computation of the change of the wavefunction for a given
          !k point.

           weight = wk (ikk)
           do iw = 1, nfs
              call incdrhoscf(drhoscf(1, current_spin, iw), weight, ik, &
                              dbecsum(1,1,current_spin))
           enddo!iw
        enddo !ibnd
     enddo 

        if (doublegrid) then
             do is = 1, nspin_mag
                call cinterpolate (drhoscfh(1,is,1), drhoscf(1,is,1), 1)
             enddo
        else
                call zcopy (nspin_mag*nrxx, drhoscf, 1, drhoscfh, 1)
        endif

!No USPP:
!    call addusddens (drhoscfh, dbecsum, imode0, npe, 0)
     call zcopy (nrxx*nspin_mag, drhoscfh(1,1,1), 1, dvscfout(1,1,1), 1)

     ! SGW: here we enforce zero average variation of the charge density
     ! if the bare perturbation does not have a constant term
     ! (otherwise the numerical error, coupled with a small denominator
     ! in the coulomb term, gives rise to a spurious dvscf response)
     ! One wing of the dielectric matrix is particularly badly behaved 

     meandvb = sqrt ( (sum(dreal(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0 ) / float(nrxxs)
     do iw = 1, nfs
        if (meandvb.lt.1.d-8) then 
            call cft3 (dvscfout(1,1,iw), nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
            dvscfout ( nl(1), current_spin, iw) = dcmplx(0.d0, 0.0d0)
            call cft3 (dvscfout(1,1,iw), nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        endif 
     enddo

     do iw = 1, nfs
        call dv_of_drho (1, dvscfout(1,1,iw), .true.)
     enddo

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
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif
     tcpu = get_clock ('GW')
     dr2 = dr2 / DBLE(npe)
     CALL flush_unit( stdout )
     rec_code=10
  enddo !loop on kter (iterations)
155 iter0=0
!   -vc*\Chi
    do iw = 1 , nfs
       drhoscf(:,1,iw) = -dvscfout(:,1,iw)
    enddo

    if (convt) then
    if (fildvscf.ne.' ') then
    write(6, '("fildvscf")') 
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
END SUBROUTINE solve_direct_shift

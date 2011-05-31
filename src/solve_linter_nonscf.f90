!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter_nonscf(igpert, drhoscf)
!-----------------------------------------------------------------------
!
!HL- This non-scf version just does a single iteration of the linear
! system solver.
!
!  Driver routine for the solution of the linear system which
!  defines the change of the wavefunction due to a perturbing potential
!  parameterized in r and w. i.e. v_{r,w} (r')
!  It performs the following tasks:
!   a) computes the bare potential term Delta V | psi > 
!   b) adds to it the screening term Delta V_{SCF} | psi >
!   c) applies P_c^+ (orthogonalization to valence states)
!   d) calls cgsolve_all to solve the linear system
!   e) computes Delta rho, Delta V_{SCF} and symmetrizes them... 
!   Currently symmetrized in terms of mode etc. Might need to strip this out
!   and check PW for how it stores/symmetrizes charge densities.
!

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin
  USE cell_base,            ONLY : tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : lgauss, degauss, ngauss, xk, wk
  USE gvect,                ONLY : nrxx, g, nl
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
                                   ext_recover
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iunrec, iudvscf, &
                                   this_pcxpsi_is_on_file
  USE output,               ONLY : fildrho, fildvscf
  USE gwus,                 ONLY : int3_paw, becsumort
  USE eqv,                  ONLY : dvpsi, dpsi, evq, eprec
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE modes,                ONLY : npertx, npert, u, t, irotmq, tmq, &
                                   minus_q, irgq, nsymq, rtau 
  USE recover_mod,          ONLY : read_rec, write_rec

  ! used oly to write the restart file
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none

  integer :: irr, imode0, npe
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  complex(DP) :: drhoscf (nrxx, nspin_mag)
  ! output: the change of the scf charge

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
       dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  REAL(DP), allocatable :: becsum1(:,:,:)
  
!HL temp array so I can look at fourier components of drho.
  complex(DP), allocatable :: drhoaux (:,:)
  COMPLEX (DP), ALLOCATABLE :: hpsi(:,:)

  logical :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0        ! true if xq=(0,0,0) in a metal

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

  external ch_psi_all, cg_psi
  
  IF (rec_code_read > 20 ) RETURN

!HL- Allocate arrays for dV_scf (need to alter these from (nrxx, nspin_mag, npe) to just (nrxx, nspin_mag).

  npe = 1
  imode0 = 1
  irr = 1
  ipert = 1

  allocate (hpsi(npwx*npol, 4)) ! Test array for whether linear system is being properly solved.
  call start_clock ('solve_linter')
  allocate (dvscfin ( nrxx , nspin_mag))    
  if (doublegrid) then
     allocate (dvscfins ( nrxxs , nspin_mag))    
  else
     dvscfins => dvscfin
  endif

  allocate (drhoscfh ( nrxx , nspin_mag))    
  allocate (dvscfout ( nrxx , nspin_mag))    
  allocate (drhoaux ( nrxx , nspin_mag))
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag, npe))    
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
        CALL read_rec(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh)
     ENDIF
     rec_code=0
  else
    iter0 = 0
    convt =.FALSE.
    where_rec='no_recover'
  endif

  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL DIROPN (iudrho, TRIM(fildrho)//'.u', lrdrho, exst)
  END IF

  IF (convt) GOTO 155

  ! In this case it has recovered after computing the contribution
  ! to the dynamical matrix. This is a new iteration that has to 
  ! start from the beginning.

  IF (iter0==-1000) iter0=0

  !
  ! The outside loop is over the iterations niter_gw maximum number of iterations
  !   

  !if(okvan) WRITE(6,'("OKVAN")')

!  do kter = 1, niter_gw
     kter = 1
     iter = kter + iter0
     ltaver = 0
     lintercall = 0
     drhoscf(:,:) = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     if (nksq.gt.1) rewind (unit = iunigk)

!START LOOP OVER KPOINTS
     do ik = 1, nksq
        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_linter', 'reading igk', abs (ios) )
        endif

! lgamma is a q=0 computation
        if (lgamma)  npwq = npw
! k and k+q mesh defined in initialize_gw
        ikk = ikks(ik)
        ikq = ikqs(ik)

!      WRITE(stdout, '("ikk  ", i4, " ikq ", i4, " ik ", i4, " nksq ", i4)')ikk, ikq, ik, nksq

        if (lsda) current_spin = isk (ikk)
        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('solve_linter', 'reading igkq', abs (ios) )
        endif

        !   Calculates beta functions (Kleinman-Bylander projectors), with
        !   structure factor, for all atoms, in reciprocal space
        !   HL the beta functions (vkb) are being generated properly.  

        call init_us_2 (npwq, igkq, xk (1, ikq), vkb)

        !
        ! reads unperturbed wavefuctions psi(k) and psi(k+q)
        !
        if (nksq.gt.1) then
           if (lgamma) then
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
           else
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
              call davcio (evq, lrwfc, iuwfc, ikq, - 1)
           endif
        endif

!DEBUG the eigenvectors at k and k+q should be exactly the same as in the PH sample code here
!        do ig = 1, ngm
!        WRITE (stdout, '("psi_k   ",  3f7.4)')evc(ig, 4)
!        WRITE (stdout, '("psi_q   ",  3f7.4)')evq(ig, 4)
!        end do
!        stop
!        compute the kinetic energy
!HL- I checked this. It is reading in the eigenfunctions at this point just right. Through all iterations

        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                          (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                          (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
        enddo

!HL gvectors/kvectors are also being read properly...

        h_diag=0.d0
        do ibnd = 1, nbnd_occ (ikk)
           do ig = 1, npwq
              h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
           enddo
           IF (noncolin) THEN
              do ig = 1, npwq
                 h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
              enddo
           END IF
        enddo

        !
        ! diagonal elements of the unperturbed hamiltonian
        ! do ipert = 1, npe
        !   mode = imode0 + ipert
        !   nrec = (ipert - 1) * nksq + ik

            mode = 1
            nrec = ik
           !
           !  and now adds the contribution of the self consistent term
           !

           if (where_rec =='solve_lint'.or.iter>1) then
              !
              ! After the first iteration dvbare_q*psi_kpoint is read from file

              call davcio (dvpsi, lrbar, iubar, nrec, - 1)

              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !

              call start_clock ('vpsifft')
              do ibnd = 1, nbnd_occ (ikk)
                 call cft_wave (evc (1, ibnd), aux1, +1) 
                 call apply_dpot(aux1, dvscfins(1,1), current_spin)
                 call cft_wave (dvpsi (1, ibnd), aux1, -1)
              enddo

              call stop_clock ('vpsifft')

        !
        !  In the case of US pseudopotentials there is an additional 
        !  selfconsist term which comes from the dependence of D on 
        !  V_{eff} on the bare change of the potential
        !
        ! HL- adddvscf now needed to include the augmentation charges and core charges.
        ! This routine computes the contribution of the selfconsistent
        ! change of the potential (i.e. of the augmentation charge  to the known part of the linear
        ! system and adds it to dvpsi.
        ! It implements the second term in Eq. B30 of PRB 64, 235118 (2001). 

           call adddvscf (ipert,ik)

           else
              !  At the first iteration dvbare_q*psi_kpoint is calculated
              !  and written to file
              !  call dvqpsi_us (ik, 1, .false.)

            call dvqpsi_us (igpert, ik, 1, .false.)
            call davcio (dvpsi, lrbar, iubar, nrec, 1)

           endif

           !
           ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
           ! Apply -P_c^+.
           ! -P_c^ = - (1-P_v^)  need to check consistency of GW routine and this one:
           ! SGW := call emptyproj ( evq, dvpsi)
           ! 

           CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi)
           !
           if (where_rec=='solve_lint'.or.iter > 1) then

              !
              ! starting value for delta_psi is read from iudwf
              !

              nrec1 = ik
              call davcio ( dpsi, lrdwf, iudwf, nrec1, -1)
              !
              ! threshold for iterative solution of the linear system
              !

              thresh = min (1.d-1 * sqrt (dr2), 1.d-2)

           else

              !
              !  At the first iteration dpsi and dvscfin are set to zero
              !

              dpsi(:,:) = (0.d0, 0.d0) 
              dvscfin (:, :) = (0.d0, 0.d0)

              !
              ! starting threshold for iterative solution of the linear system
              !

              thresh = 1.0d-2
           endif

           conv_root = .true.

!        call bcgsolve_all   (ch_psi_all_eta, et(:,ikk), dvpsi, dpsip, h_diag, &
!             ngm, ngm, tr_cgsolve_now, ik, lter, conv_root, anorm, nbnd_occ, &
!             g2kin, vr, evq,  cw )
!          iterative solution of the linear system (H-eS)*dpsi=dvpsi,
!          dvpsi = -P_c^+ (dvbare+dvscf)*psi , dvscf fixed.

           call cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, dpsi, &
                             h_diag, npwx, npwq, thresh, ik, lter, conv_root, &
                             anorm, nbnd_occ(ikk), npol )

! HL - DEBUG same as in SGW. Going to check that (H-e_{vk}+alphaPv*)dpsi + Pcdvpsi = 0
! For some reason in SGW they only check (H-et+alpha*Pv)*dpsi-dvpsi = 0 
! Reason = dvpsi already projected onto valence states?
!          WRITE(stdout, '("(H-et+alpha*Pv)*dpsi-dvpsi=0")')
!          do ibnd = 1, nbnd_occ(ikk)
!             call ch_psi_all(ndim, dpsi, hpsi, et(1,ikk), ikk, nbnd_occ(ikk))
!             hpsi(:,1)  = hpsi(:,1) - dvpsi(:,1) 
!             do ig = 1, npwq
!              WRITE(stdout, '("LHS - RHS ", 3f15.5 )'), hpsi(ig,1) 
!             enddo
!          enddo
!-HL
           ltaver = ltaver + lter
           lintercall = lintercall + 1
           if (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4,  &
                &              " solve_linter: root not converged ",e10.3)') &
                &              ik , ibnd, anorm
           ! writes delta_psi on iunit iudwf, k=kpoint,
           ! nrec1 = (ipert - 1) * nksq + ik
           nrec1 =  ik
           call davcio (dpsi, lrdwf, iudwf, nrec1, + 1)
           !
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
                            dbecsum(1,1,current_spin,ipert))
           END IF
     enddo ! END LOOP ON K-POINTS

!HL debug by printing out fourier components of drhoscf.
!         drhoaux(:,1) = drhoscf(:,1)  
!         call cft3 (drhoaux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)
!         write(6, '("")')
!         do ig = 1, 3
!         write (6, '("drhoaux(nl(ig)) ", 3f7.4, 3f7.3)' )drhoaux(nl(ig),1)
!         write (6, '("drhoscf(nl(ig)) ", 3f7.4, 3f7.3)' )drhoscf(nl(ig),1)
!         write (6, '("drhscfh(nl(ig)) ", 3f7.4, 3f7.3)' )drhoscfh(nl(ig),1)
!         end do
!
!end drho debug

! call addusddens (drhoscfh, dbecsum, imode0, npe, 0)

#ifdef __PARA
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_pool_comm )
     ELSE
        call mp_sum ( dbecsum, intra_pool_comm )
     ENDIF
#endif

     if (doublegrid) then
        do is = 1, nspin_mag
              call cinterpolate (drhoscfh(1,is), drhoscf(1,is), 1)
        enddo
     else
        call zcopy (nspin_mag*nrxx, drhoscf, 1, drhoscfh, 1)
     endif
     !
     !  In the noncolinear, spin-orbit case rotate dbecsum
     !
#ifdef __PARA
     !
     !   Reduce the delta rho across pools
     !
     call mp_sum ( drhoscf, inter_pool_comm )
     call mp_sum ( drhoscfh, inter_pool_comm )
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
#endif
     !   q=0 in metallic case deserve special care (e_Fermi can shift)
     !   After the loop over the perturbations we have the linear change 
     !   in the charge density for each mode of this representation. 
     !   Here we symmetrize them ...
     IF (.not.lgamma_gamma) THEN
#ifdef __PARA

!HL Turning symmetry off for the moment
!       call psymdvscf (npe, irr, drhoscfh)

        IF ( noncolin.and.domag ) &
           CALL psym_dmag( npe, irr, drhoscfh)
#else
        call symdvscf (npe, irr, drhoscfh)
        IF ( noncolin.and.domag ) CALL sym_dmag( npe, irr, drhoscfh)
#endif
        IF (okpaw) THEN
           IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,npe,irr, &
                             npertx,irotmq,rtau,xq,tmq)
           CALL  &
                PAW_dusymmetrize(dbecsum,npe,irr,npertx,nsymq,irgq,rtau,xq,t)
        END IF
     ENDIF

     ! 
     !   ... save them on disk and 
     !   compute the corresponding change in scf potential 
     !
   

      if (fildrho.ne.' ') call davcio_drho (drhoscfh(1,1), lrdrho, &
                                            iudrho, imode0+ipert, +1)
      call zcopy (nrxx*nspin_mag,drhoscfh(1,1),1,dvscfout(1,1),1)

!  HL - old version with modes and perts... call dv_of_drho (imode0+ipert, dvscfout(1,1), .true.)

      call dv_of_drho (1, dvscfout(1,1), .true.)

!   do ig= 1, 20
!      WRITE(stdout,'(4x,4x,"dvdrho = ", 2f9.5)') dvscfout(nl(ig), nspin_mag)
!   enddo
!   And we mix with the old potential

     call mix_potential (2*npe*nrxx*nspin_mag, dvscfout, dvscfin, &
                         alpha_mix(kter), dr2, npe*tr2_gw/npol, iter, &
                         nmix_gw, flmixdpot, convt)

!HL dvscfin is the quantity multiplied by the appropriate prefactor fpi e2 ... 
!In the non scf case is think this number should be 12.
!do ig = 1, 50
!    WRITE (stdout,'("dvscfin(R(1)) =  ", 3f10.5)' ) -real(dvscfin( nl(1),1))

     call cft3s( dvscfin, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
     WRITE (stdout,'("dvscfin (G(1)) =  ", 3f10.5)' ) -real(dvscfin( nl(1),1))

!    call cft3s( dvscfin, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!    WRITE (stdout,'("dvscfin(R(1)) =  ", 3f10.5)' )  -real(dvscfin( nl(1),1))
!enddo

     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin_mag
              call cinterpolate (dvscfin(1,is), dvscfins(1,is), -1)
           enddo
        enddo
     endif

!   calculate here the change of the D1-~D1 coefficients due to the PH
!   perturbation

     IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
     !
     ! with the new change of the potential we compute the integrals
     ! of the change of potential and Q 
     !
     ! HL-Q denotes the augmentation charge. Look at Vanderbilt PRB 41 7892
     ! Q_{ij} = <psi_{j}|psi_{i}> - <phi_{i}|phi_{j}> 
     ! USPP valence charge density is described 
     ! n_v(r) =  \sum_{n,k} \phi^{*}_{nk} (r) \phi_{nk}(r) + \sum_{i,j} p_{i,j}Q_{j,i}
     ! p_{i,j} = \sum_{n,k} <beta_{i}|\phi_{nk}><phi_{nk}|beta_{j}> 
     ! HL newdq ULTRASOFT routine  
     !     call newdq (dvscfin)

#ifdef __PARA
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)
     call mp_sum ( aux_avg, inter_pool_comm )
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif

     tcpu = get_clock ('GW')
     dr2 = dr2 / npe

!   WRITE (stdout, '(/,5x," iter ",i3," cpu ", f8.1, " ndim ", i4)')iter, tcpu, ndim
!   WRITE (stdout, '(5x," kter", i4)')kter 
!   WRITE (stdout, '(5x," alpha mix", f6.3)')alpha_mix (kter) 
!   WRITE (stdout, '(5x,"dr2 ",e10.3," npe ", i4)')dr2, npe
!   WRITE (stdout, '(5x," thresh ", e10.3)')thresh  
!   WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
!   "secs   av.it.: ",f5.1)') iter, tcpu, averlt
!   WRITE( stdout, '(5x," thresh=",e10.3, " alpha_mix = ",f6.3, &
!   "|ddv_scf|^2 = ",e10.3 )') thresh, alpha_mix (kter) , dr2
!   if (convt)  write(stdout, '("CONVERGED")')
!   Here we save the information for recovering the run from this point
!   HL-  CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, dvscfin, drhoscfh) 
!   To implement rewrite function for GW might want something like:
!   write_rec('solve_linte, k, dr2, iter, cont, q, dvscfin, drhoscfh) 

     CALL flush_unit( stdout )

     rec_code=10
     IF (okpaw) THEN
        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                               dvscfin, drhoscfh, dbecsum)
     ELSE

        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                               dvscfin, drhoscfh)
     ENDIF

     if (check_stop_now()) call stop_smoothly_gw (.false.)

!DEBUG
!-HL WRITE(stdout, '("(H-et+alpha*Pv)*dpsi-dvpsi=0")')
!          if (convt) then 
!          do ibnd = 1, nbnd_occ(1)
!             call ch_psi_all(ndim, dpsi, hpsi, et(1,ikk), ikk, nbnd_occ(ikk))
!             call ch_psi_all(ndim, dpsi, hpsi, et(1,1), 1, nbnd_occ(1))
!             hpsi(:,1)  = hpsi(:,1) - dvpsi(:,1) 
!             do ig = 1, npwq
!              WRITE(stdout, '("LHS - RHS ", 3f15.5 )'), hpsi(ig,1) 
!             enddo
!          enddo
!          goto 155
!          endif
!-HL moved if cont loop into above if block
!     if (convt) goto 155
!  enddo !loop on kter (iterations)
!END DEBUG


155 iter0=0

  !    A part of the dynamical matrix requires the integral of
  !    the self consistent change of the potential and the variation of 
  !    the charge due to the displacement of the atoms. 
  !    We compute it here. 

  if (convt) then
   !HL add the contribution drhodvus to dynamical matrix  call drhodvus (irr, imode0, dvscfin, npe)
     if (fildvscf.ne.' ') then
     write(6, '("fildvscf")') 
    !call davcio_drho ( dvscfin(1,1),  lrdrho, iudvscf, imode0 + ipert+(current_iq-1)*3*nat, +1 )
     call davcio_drho ( dvscfin(1,1),  lrdrho, iudvscf, 1, +1 )
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
END SUBROUTINE solve_linter_nonscf

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


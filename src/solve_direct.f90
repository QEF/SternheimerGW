!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
SUBROUTINE solve_lindir(dvbarein, iw, drhoscf)
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
  USE cell_base,            ONLY : tpiba2,at,bg
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
                                   ext_recover, eta, maxter_green
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iunrec, iudvscf, iudwfm, iudwfp 
  USE output,               ONLY : fildrho, fildvscf
  USE gwus,                 ONLY : int3_paw, becsumort
! USE eqv,                  ONLY : dvpsi, dpsi, evq, eprec, dpsim, dpsip
  USE eqv,                  ONLY : dvpsi, evq, eprec
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

  integer :: iw, allocatestatus
  integer :: irr, imode0, npe

  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes
  complex(DP) :: drhoscf (nrxx, nfs)
  ! output: the change of the scf charge
  complex(DP) :: dvbarein (nrxxs)

! HL prec
! HL careful now... complexifying preconditioner:
  real(DP) , allocatable :: h_diag (:,:)
! h_diag: diagonal part of the Hamiltonian
! complex(DP) , allocatable :: h_diag (:,:)
  real(DP) :: thresh, anorm, averlt, dr2
  real(DP) :: x, a, b, norm
  real(DP) :: xkloc(3)

  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight, aux_avg (2)

  ! Misc variables for metals
  ! dos_ef: density of states at Ef
  real(DP), external :: w0gauss, wgauss
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:), dvscfout (:,:)

  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  complex(DP), allocatable :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
                              dbecsum (:,:,:), dbecsum_nc(:,:,:,:,:)
  complex(DP) :: cw
  complex(DP), allocatable :: etc(:,:)


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

  INTEGER                   :: gveccount

! COMPLEX(DP), ALLOCATABLE :: dpsic(:,:,:), dpsit(:,:,:), dpsi(:,:,:)
! complex(DP), allocatable :: alphabeta(:,:,:)
! INTEGER, ALLOCATABLE      :: niters(:)

  INTEGER     :: niters(nbnd)
  COMPLEX(DP) :: dpsic(npwx,nbnd,maxter_green+1), dpsit(npwx, nbnd, nfs), dpsi(npwx,nbnd,nfs)
  COMPLEX(DP) :: alphabeta(npwx,nbnd,maxter_green+1)
 
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
  allocate (dvscfout ( nrxx , nfs))    
  allocate (drhoscfh ( nrxx , nfs))    
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat, nspin_mag))    
!Complex eigenvalues:
  allocate (etc(nbnd, nkstot))
  allocate (h_diag ( npwx*npol, nbnd))    

!Multishift arrays.
!  allocate (niters(nbnd))
!  allocate (dpsit ( npwx, nbnd, nfs))
!  allocate (dpsi  ( npwx, nbnd, nfs))
!  allocate (dpsic ( npwx, nbnd, maxter_green+1))
!  allocate (alphabeta ( 2, nbnd, maxter_green+1))

  iter0 = 0
  convt =.FALSE.
  where_rec='no_recover'

  IF (iter0==-1000) iter0=0

!No self-consistency:
  do kter = 1, 1
     iter = kter + iter0
     ltaver = 0

     lintercall = 0
     drhoscf(:,:)   = (0.d0, 0.d0)
     drhoscfh(:,:)  = (0.d0, 0.d0)
     dbecsum(:,:,:) = (0.d0, 0.d0)

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
!MULTISHIFT No Preconditioning.
        h_diag = 0.d0
        do ibnd = 1, nbnd_occ (ikk)
           do ig = 1, npwq
              h_diag(ig,ibnd) =  (1.0d0, 0.0d0)
           enddo
        enddo
!mode relic from phonon days:
        mode = 1
        nrec = ik

             call dvqpsi_us (dvbarein, ik, 1, .false.)

! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
! Apply -P_c^+.
! -P_c^ = - (1-P_v^):

             CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi(:,:,iw))

             dpsic(:,:,:)   =  dcmplx(0.d0, 0.d0)
             dpsit(:,:,:)   =  dcmplx(0.d0, 0.d0)
             dpsi(:,:,:)    =  dcmplx(0.d0, 0.d0)

             thresh    = tr2_gw
             conv_root = .true.

             etc(:,:)  = CMPLX(et(:,:), 0.0d0 , kind=DP)


             call cbcg_solve_coul(cch_psi_all_fix, cg_psi, etc(1,ikk), dvpsi, dpsi(:,:,:), dpsic, h_diag, &
                                  npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), &
                                  npol, niters(:), alphabeta(:,:,:))
!            dpsi = dpsi^{+}

            call coul_multishift(npwx, npwq, nfs, niters, dpsit, dpsic, alphabeta, fiu)
            dpsi(:,:,:) = dpsit(:,:,:) 

!            dpsi = dpsi^{+} + dpsi^{-}
             dpsit(:,:,:) = dcmplx(0.0d0, 0.0d0)

            call coul_multishift(npwx, npwq, nfs, niters, dpsit, dpsic, alphabeta, -fiu)
            dpsi(:,:,:) = dpsi(:,:,:) + dpsit(:,:,:)


             ltaver = ltaver + lter
             lintercall = lintercall + 1
             nrec1 =  ik
             weight = wk (ikk)

           do iw = 1 , nfs
                 call incdrhoscf_w (drhoscf(1, iw) , weight, ik, &
                                    dbecsum(1,1,current_spin), dpsi(1,1,iw))
           enddo
     enddo !kpoints

     print*, "finished kpoints"

!     deallocate (niters)
!     deallocate (alphabeta)
!     deallocate (dpsi)
!     deallocate (dpsic)
!     deallocate (dpsit)


     call zcopy (nspin_mag*nrxx, drhoscf, 1, drhoscfh, 1)
     call addusddens (drhoscfh, dbecsum, imode0, npe, 0)
     call zcopy (nrxx*nspin_mag, drhoscfh(1,1),1, dvscfout(1,1),1)

     ! SGW: here we enforce zero average variation of the charge density
     ! if the bare perturbation does not have a constant term
     ! (otherwise the numerical error, coupled with a small denominator
     ! in the coulomb term, gives rise to a spurious dvscf response)
     ! One wing of the dielectric matrix is particularly badly behaved 

     meandvb = sqrt ( (sum(dreal(dvbarein)))**2.d0 + (sum(aimag(dvbarein)))**2.d0 ) / float(nrxxs)
     do iw = 1, nfs 
         if (meandvb.lt.1.d-8) then 
             call cft3 (dvscfout, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
             dvscfout  (nl(1), iw) = dcmplx(0.d0, 0.0d0)
             call cft3 (dvscfout, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
         endif
    enddo

     call dv_of_drho (1, dvscfout(1,1), .true.)

     averlt = DBLE (ltaver) / lintercall
     tcpu = get_clock ('GW')
     dr2 = dr2 / DBLE(npe)
     CALL flush_unit( stdout )
     rec_code=10
  enddo !loop on kter (iterations)

!   after this point drhoscf is dv_hartree(RPA)
!   drhoscf(:,1) = dvscout(:,1)
!  -vc*/Chi
    drhoscf(:,:) = -dvscfout(:,:)

  if (convt) then
   if (fildvscf.ne.' ') then
    write(6, '("fildvscf")') 
   end if
  endif

  deallocate (h_diag)
  deallocate (dvscfout)
  deallocate (drhoscfh)
  deallocate (dbecsum)

  call stop_clock ('solve_linter')
END SUBROUTINE solve_lindir


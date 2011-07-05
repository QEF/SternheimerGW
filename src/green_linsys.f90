SUBROUTINE green_linsys(ik0, iq)
! SGW: call green_linsys ( vr, g2kin, k0mq, nwgreen, wgreen, igstart, igstop, ik0, iq )
! "In order to calculate the analytic component G^{A} we consider G^{A}(r,r',w) as a
! parametric function of the first space variable and of the frequency: G^{A}_[r,w](r') 
! G^{A} = \sum_{n}\psi_n(r)\psi_n(r')/(w - ev + idelta)    (16)
! (H - w^{+})G^{A}_[r,w] = -\delta_[r](r')                 (18)
! let n -> nk so we are dealing with bloch states and expand the functions in terms of
! the plane waves: 
! G_[r,w]^{A}(r') = \frac{1}{(N_{k}\Omega)}\sum_{kg}g^{A}_[k,\omega,G](\r')e^{-i(k + G)\r} e^{ik\r'}
! The equation of motion then becomes:
! (H_{k} - \omega^{+})g_{[k,G,\omega]}(G') = -\delta_{GG'}

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin, pi, tpi, RYTOEV
  USE cell_base,            ONLY : tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : lgauss, degauss, ngauss, xk, wk, nkstot
  USE gvect,                ONLY : nrxx, g, nl, ngm
  USE gsmooth,              ONLY : doublegrid, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, ngms
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,             ONLY : domag
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et
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
  USE units_gw,             ONLY : iuwfc, lrwfc, iuwfcna, iungreen
  USE eqv,                  ONLY : evq, eprec
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE recover_mod,          ONLY : read_rec, write_rec
! used oly to write the restart file
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
!
  USE disp,        ONLY : nqs
  USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax, nwgreen, wgreen
  USE gwsigma,     ONLY : ngmsig

  IMPLICIT NONE 

  real(DP) :: thresh, anorm, averlt, dr2
  logical :: conv_root

  !HL Again need to think about the dimensions of the green's fxn.
  !It is a matrix of the dimension of the plane wave description before
  !truncation. The second index is a dummy band index although it could come in handy
  !if we started to look at spin systems.  

  COMPLEX(DP) :: gr_A(npwx,1), rhs(npwx, 1)
  COMPLEX(DP) :: gr_N(npwx,1), gr(npwx,1), ci, cw, green(ngmsig,ngmsig)
  COMPLEX(DP), ALLOCATABLE :: etc(:,:)

  INTEGER :: iw, igp
  INTEGER :: iq, ik0
  INTEGER :: rec0, n1
  REAL(DP) :: dirac, x, delta
  real(DP) :: k0mq(3) 
  !HL eta should probably be a constant user defined parameter as well...
  real(DP) :: eta 
  real(DP), allocatable :: et2(:,:)

  real(DP) :: w_ryd(nwgreen)
  external cch_psi_all, ch_psi_all, cg_psi, cch_psi_all_fix, cch_psi_all_green
  real(DP) , allocatable :: h_diag (:,:)
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
             mode          ! mode index


!Parameter for the dirac delta function
!real(DP), parameter :: eta = 0.04
!@10TION
!HL need a threshold here for the linear system solver. This could also go in the punch card
!with some default at a later date. 

        REAL(DP) :: tr_cgsolve = 1.0d-5
        allocate (h_diag ( npwx*npol, nbnd))
        allocate (etc(nbnd, nkstot))
  
!dummy array
  
        ci = (0.0d0, 1.0d0)
        eta = 0.04d0

!convert freq array generated in freqbins into rydbergs.
        w_ryd(:) = wgreen(:)/RYTOEV
        CALL start_clock('greenlinsys')

! WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,ikq)*RYTOEV
! generally the eigenvalues/for each kpoint we want will be stored somewhere in the array xk
! i.e. if we want to look at a specific k-point we need to know where it sits in that array.

! for now I am only interested in k = Gamma, and therefore k-q will always be in the second spot in all these arrays.
! When I want to start looking at a few k-points in the BZ then things might become a little bit trickier and 
! could try to foil me. 

! ikks(ik) = 2*ik. ikqs(ik) = 2*ik. ik=1=Gamma.

!Again
!        ikk = ikks(ik)
!        ikq = ikqs(ik)

        ikq = 2
        call davcio (evq, lrwfc, iuwfc, ikq, - 1)

! Various plane wave variables.
! I never checked to make sure that the eigenfxns, kinetic energy, pseudo-potentials, 
! form factors etc. were being generated properly. This is a very silly state of affairs.
! Tomorrow need to check all these things carefully so that green_linsys is a bit more like
! solve_linter. But it certainly appears necessary to have these variables set properly for
! when h_psiq is called and the linear system is run. 

        if (nksq.gt.1) rewind (unit = iunigk)

        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('green_linsys', 'reading igk', abs (ios) )
        endif

        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('green_linsys', 'reading igkq', abs (ios) )
        endif

        !   Calculates beta functions (Kleinman-Bylander projectors), with
        !   structure factor, for all atoms, in reciprocal space

       call init_us_2 (npwq, igkq, xk (1, ikq), vkb)

       do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                          (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                          (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
       enddo

       h_diag = 0.d0
       do ibnd = 1, nbnd_occ (ikq)
          do ig = 1, npwq
             h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ikq))
          enddo
          IF (noncolin) THEN
             do ig = 1, npwq
                h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ikq))
             enddo
          END IF
       enddo

!@10TION
! Here it becomes important that all the valence band eigenvalues are < 0.
!  
! et(:,:) = et(:,:) - (6.58D0/13.605)

WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,2)*RYTOEV
DO iw = 1, nwgreen
    green =(0.0d0, 0.0d0)
!    do ig = 1, ngmsig !zero after this point
     do ig = 1, 1 
       rhs(:,:)  = (0.0d0, 0.0d0)
       rhs(ig,1) = -(1.0d0, 0.0d0)
       gr_A(:,:) = (0.0d0, 0.0d0) 
       conv_root = .true.
       lter = 0

!      anorm = 0.0d0
       etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
       cw = CMPLX( w_ryd(iw), eta, kind=DP)
!      cw = CMPLX( 0.0d0, w_ryd(iw), kind=DP)

!      cw = CMPLX( 0.0d0, w_ryd(iw), kind=DP)
!      call cbcg_solve(cch_psi_all, cg_psi, etc(:,ikq), rhs, gr_A, h_diag, &
!           npwx, npwq, tr_cgsolve, ikq, lter, conv_root, anorm, 1, npol, .false.)
!      call cbcg_solve_fix(cch_psi_all_fix, cg_psi, etc, rhs, gr_A, h_diag, &
!            npwx, npwq, 1.0d-2, 1, lter, conv_root, anorm, 1, npol, cw, .true.)
!      call cgsolve_all (ch_psi_all, cg_psi, w_ryd(iw), rhs, gr_A, &
!                         h_diag, npwx, npwq, 1.0d-2, ikq, lter, conv_root, &
!                         anorm, nbnd_occ(ikq), npol )

       call cbcg_solve_fix(cch_psi_all_green, cg_psi, etc(1,ikq), rhs, gr_A, h_diag, &
             npwx, npwq, tr_cgsolve, ikq, lter, conv_root, anorm, 1, npol, cw, .false.)

       IF (.not.conv_root) then
        WRITE(6,'(" Root Not Converged ")')
!       WRITE( stdout, '(/,5x," iter # ",i3," anorm :",f8.3, &
!      "gr_A(1): ",f8.3)') lter, anorm, real(gr_A(ig))
!       gr_A = (0.d0, 0.0d0)
       ENDIF

      gr_N(:,1) = (0.0d0, 0.0d0)

      do ibnd = 1, 4
!       x = w_ryd(iw) - eval(ibnd)
        x = w_ryd(iw) - et(ibnd, ikq)
!        x = (0.0d0, 1.0d0)*w_ryd(iw) - et(ibnd, ikq)
        dirac = eta / pi / (x**2.d0 + eta**2.d0)

!      no spin factor (see above) 
! HL   I might need to divide this by two since q-espresso has the spin factor in the kweights.
!      gr_N = gr_N + 2.d0 * twopi * ci * conjg ( psi(ig,ibnd) ) * psi(1:ngm,ibnd) * dirac
!      gr_N = gr_N + tpi * ci * conjg( evq(ig,ibnd) ) * evq(1:ngm,ibnd) * dirac
!      do igp = 1, npwq

       do igp = 1, ngmsig
!       do igp = 1, 1
!       SGW  gr_N(ig) = gr_N(ig) + tpi * ci * conjg( evq(ig,ibnd) ) * evq(igp,ibnd) * dirac
             gr_N(ig,1) = gr_N(ig,1) + tpi * ci * conjg(evq(ig,ibnd)) *  evq(igp,ibnd) * dirac
       enddo
      enddo

       gr = gr_A + gr_N
!       gr = gr_A

      do igp = 1, ngmsig
         green (ig,igp) = gr (igp,1)
         if (ig.eq.1.and.igp.eq.1) write(6,'(3f15.10)') wgreen(iw), green(ig,igp)
      enddo
   enddo

!Collect G vectors across processors and then write the full green's function to file. 
!#ifdef __PARA
!    use poolreduce to bring together the results from each pool
!    call poolreduce ( 2 * ngms * ngms, green)
!    if (me.eq.1.and.mypool.eq.1) then
!#endif
!  rec0 = (iw-1) * nk0 * nq + (ik0-1) * nq + (iq-1) + 1
   rec0 = (iw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
   write ( iungreen, rec = rec0, iostat = ios) green
!#ifdef __PARA
!    endif
!#endif
ENDDO !enddo on iw
STOP
CALL stop_clock('greenlinsys')
RETURN
END SUBROUTINE green_linsys

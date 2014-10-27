SUBROUTINE green_linsys_test(ik0, iq)

!SGW: call green_linsys ( vr, g2kin, k0mq, nwgreen, wgreen, igstart, igstop, ik0, iq )
!"In order to calculate the analytic component G^{A} we consider G^{A}(r,r',w) as a
!parametric function of the first space variable and of the frequency: G^{A}_[r,w](r')
!G^{A} = \sum_{n}\psi_n(r)\psi_n(r')/(w - ev + idelta)    (16)
!(H - w^{+})G^{A}_[r,w] = -\delta_[r](r')                 (18)
!let n -> nk so we are dealing with bloch states and expand the functions in terms of
!the plane waves: 
!G_[r,w]^{A}(r') = \frac{1}{(N_{k}\Omega)}\sum_{kg}g^{A}_[k,\omega,G](\r')e^{-i(k + G)\r} e^{ik\r'}
!The equation of motion then becomes:
!(H_{k} - \omega^{+})g_{[k,G,\omega]}(G') = -\delta_{GG'}

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
  ! used only to write the restart file
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  USE disp,        ONLY : nqs
  USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax, nwgreen, wgreen
  USE gwsigma,     ONLY : ngmsig

  IMPLICIT NONE 

  real(DP) :: thresh, anorm, averlt, dr2
  logical :: conv_root

  COMPLEX(DP) :: gr_A(ngm), gr_N(ngm), rhs(ngm), gr(ngm), gr_exp(nwgreen), ci, cw, green(ngmsig,ngmsig)
  COMPLEX(DP), ALLOCATABLE :: etc(:,:)

  INTEGER :: iw, igp, ig1, ig2
  INTEGER :: iq, ik0
  INTEGER :: rec0, n1
  REAL(DP) :: dirac, x, delta
  real(DP) :: k0mq(3) 
  !HL eta should probably be a constant user defined parameter as well...
  real(DP) :: eta 
  real(DP) :: w_ryd(nwgreen)
  external cch_psi_all, ch_psi_all, cg_psi

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

!HL need a threshold here for the linear system solver. This could also go in the punch card
!with some default at a later date. 

        REAL(DP) :: tr_cgsolve = 1.0d-10

        allocate (h_diag ( npwx*npol, nbnd))
        allocate (etc(nbnd, nkstot))


        ci = (0.0d0, 1.0d0)
        eta = 0.04d0

!Convert Freq array generated in freqbins into rydbergs.
        w_ryd(:) = wgreen(:)/RYTOEV

!       write(6,*) nwgreen
!       do iw = 1, nwgreen
!       write(6,*) w_ryd(iw), wgreen(iw)
!       enddo
!       this is now set in prepare_kmq
!       k0mq(:) = xk(:,ik0) - xq(:) the wave functions at k0mq should be stored in:
!       if (.not.lgamma.and.nksq.gt.1) then
!          read (iunigk, err = 200, iostat = ios) npwq, igkq
!200       call errore ('solve_linter', 'reading igkq', abs (ios) )
!       endif
!  Write this wave function \psi_{k-q} to file to be used again later when I need to form the 
!  G^{NA}v(r,r') component of the self-energy. There should be a more elegant way of doing this. 
!  I should probably sit down and make some kind of control plot 
!  for where I need to calculate the eigenvectors
!  and where I need to do fourier transforms; etc. etc. 
!  call davcio (evq, lrwfc, iuwfc, ikq, -1)
!     do ibnd = 1,1
!      do ig=1,npwx
!          write(6,*)evq(ig, ibnd)
!      enddo
!     enddo

CALL start_clock('greenlinsys')

!@10TION
!WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,ikq)*RYTOEV
! generally the eigenvalues/for each kpoint we want will be stored somewhere in the array xk
! i.e. if we want to look at a specific k-point we need to know where it sits in that array.
! for now I am only interested in k = Gamma, and therefore k-q will always be in the second spot in all these arrays.
! ikks(ik) = 2*ik. ikqs(ik) = 2*ik. ik=1=Gamma.

       ikq = 2
       call davcio (evq, lrwfc, iuwfc, ikq, - 1)

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

WRITE(6, '("nbnd ", i4)')nbnd
WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,ikq)*RYTOEV
!WRITE(6,*)ikq, nbnd, nbnd_occ
!WRITE(6,*) et(ikq,:)*RYTOEV

GOTO 124 
DO iw = 1, nwgreen
    green =(0.0d0, 0.0d0)
!HL here we parallelize over Gs  do ig = igstart, igstop
   do ig = 1, ngmsig !zero after this point

!  Solve linear system  to obtain: G^{A}_[G,w](G') 
      rhs(:)  = (0.0d0, 0.0d0)
      rhs(ig) = -(1.0d0, 0.0d0)
      gr_A    = (0.0d0, 0.0d0) 
      conv_root = .true.
      lter = 0

! Going to use same trick as in solve linter where we combine the frequency and the eigenvalue 
! and in this case a small imaginary component

      etc(:,ikq) = CMPLX( et(:,ikq) + w_ryd(iw), eta, kind=DP)

!    npwq should be the number of plane waves that describes this wave function.
!    in SGW nbnd_occ is set to 1
!     call cbcg_solve(cch_psi_all, cg_psi, etc(:,ikq), rhs, gr_A, h_diag, &
!             npwx, npwq, tr_cgsolve, ikq, lter, conv_root, anorm, nbnd_occ(ikq), npol, .true.)
!     IF (conv_root) WRITE(6,'(" Root Converged ")')
!     WRITE(6,*) iw, ig, gr_A(1)
!     Evaluate the G^{N}_{k,G,w}(G') part of the Green's function.

      gr_N(:) = (0.0d0, 0.0d0)
      do ibnd = 1, 4
!       x = w_ryd(iw) - eval(ibnd)
        x = w_ryd(iw) - et(ibnd, ikq)
        dirac = eta / pi / (x**2.d0 + eta**2.d0)

!      no spin factor (see above) 
!HL    I might need to divide this by two since q-espresso has the spin factor in the kweights.
!      gr_N = gr_N + 2.d0 * twopi * ci * conjg ( psi(ig,ibnd) ) * psi(1:ngm,ibnd) * dirac
!      gr_N = gr_N + tpi * ci * conjg( evq(ig,ibnd) ) * evq(1:ngm,ibnd) * dirac
!      do igp = 1, npwq

       do igp = 1, ngmsig
        gr_N(ig) = gr_N(ig) + tpi * ci * conjg( evq(ig,ibnd) ) * evq(igp,ibnd) * dirac
       enddo
      enddo

      gr = gr_A + gr_N

!     write(6,*)gr
!     write(stdout,*) gr

      do igp = 1, ngmsig
        green (ig,igp) = gr (igp)
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

ENDDO !enddo on iw
124 CONTINUE

!HL DEBUG explicitly summing over states to construct G_{0,G,k,w}(G') 
!To be compared with output from sternheimer calculation at the same 
!Wave vectors. 
!Again need to zero eigenvalues here. 

  et(:,ikq) = et(:,ikq) - (6.58D0/13.605)

  !ig1 = 1
  !ig2 = 1
  !gr_exp = (0.0d0, 0.0d0)
  !WRITE(6, '("Summing over nbands", i4)')nbnd
  green(:,:) = (0.0d0, 0.0d0)
  gr_exp(:) = (0.0d0, 0.0d0)
  
 !write(6,*) evq(30,:)
  do iw = 1, nwgreen
    do ig = 7, 7
      do igp = 2, 2
        do ibnd = 1, nbnd
           if (ibnd.le.4) then
             !delta = 1.d0
             delta = -1.d0 ! Only want to look at the Analytic part in the upper half plane.
           else
            delta = -1.d0
           endif

          gr_exp(iw) = gr_exp(iw) + evq(ig,ibnd)*conjg(evq(igp,ibnd)) &
                         / ( w_ryd(iw) - et(ibnd,ikq) - ci * delta * eta)

          ! green(ig,igp) = green(ig,igp) + evq(ig,ibnd)*conjg(evq(igp,ibnd)) &
          !    / ( w_ryd(iw) - et(ibnd,ikq) - ci * delta * eta)

        enddo
        !green_linsys does not include this factor of two.
        !green(ig,igp) = 2.0d0 * gr_exp(iw) 
        green(ig,igp) = gr_exp(iw) 

      enddo
    enddo
     write(6,'(3f15.10)') wgreen(iw), 2.0d0*gr_exp(iw)
    !write(6,'(3f15.10)') wgreen(iw), green(1,1)
    rec0 = (iw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
    write ( iungreen, rec = rec0, iostat = ios) green
  enddo
  CALL stop_clock('greenlinsys')
  STOP
RETURN
END SUBROUTINE green_linsys_test

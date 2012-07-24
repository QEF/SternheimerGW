SUBROUTINE green_linsys (ik0)
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin, pi, tpi, RYTOEV, eps8
  USE cell_base,            ONLY : tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : xk, wk, nkstot
  USE gvect,                ONLY : nrxx, g, nl, ngm, ecutwfc
  USE gsmooth,              ONLY : doublegrid, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, ngms
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et
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
                                   ext_recover, eta, tr2_green

  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iuwfc, lrwfc, iuwfcna, iungreen, lrgrn
  USE eqv,                  ONLY : evq, eprec
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE recover_mod,          ONLY : read_rec, write_rec

  USE mp,                   ONLY : mp_sum

  USE disp,        ONLY : nqs
  USE freq_gw,     ONLY : fpol, fiu, nfs, nfsmax, nwgreen, wgreen
  USE gwsigma,     ONLY : ngmsig, ecutsco, ngmsco
 !PARALLEL
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime, &
                         nproc_pool, nproc, me_pool, my_pool_id, npool
  USE mp,         ONLY: mp_barrier, mp_bcast, mp_sum

  IMPLICIT NONE 

  real(DP) :: thresh, anorm, averlt, dr2
  logical :: conv_root
  COMPLEX(DP) :: gr_A(npwx, 1), rhs(npwx, 1)
  COMPLEX(DP) :: aux1(npwx)
 !HLS should be green(ngmsigw, ngmsigw)
 !i.e. the lower cut off.  
 !COMPLEX(DP) :: ci, cw, green(ngmsig,ngmsig)
  COMPLEX(DP) :: ci, cw, green(ngmsig,ngmsig)
  COMPLEX(DP), ALLOCATABLE :: etc(:,:)
  REAL(DP)    :: eprecloc 
  INTEGER :: iw, igp, iwi
  INTEGER :: iq, ik0
  INTEGER :: rec0, n1
  !REAL(DP) :: dirac, x, delta
  REAL(DP) :: dirac, delta
  !COMPLEX(DP) :: x
  REAL(DP) :: x
  real(DP) :: k0mq(3) 
  real(DP) :: w_ryd(nwgreen)
 !external ccg_psi, cch_psi_all_fix, cch_psi_all_green
  external cg_psi, ccg_psi, cch_psi_all_fix, cch_psi_all_green
  real(DP) , allocatable :: h_diag (:,:)
 !COMPLEX(DP) , allocatable :: h_diag (:,:)

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
!HL need a threshold here for the linear system solver. This could also go in the punch card
!with some default at a later date. 
    REAL(DP) :: tr_cgsolve = 1.0d-4
!Arrays to handle case where nlsco does not contain all G vectors required for |k+G| < ecut
    INTEGER     :: igkq_ig(npwx) 
    INTEGER     :: igkq_tmp(npwx) 
    INTEGER     :: counter
!PARALLEL
    INTEGER :: igstart, igstop, ngpool, ngr, igs
!LINALG
    COMPLEX(DP), EXTERNAL :: zdotc
!HL NOTA BENE:
!Green's function has dimensions npwx, 1 in current notation...
!Since we store in G space there is no truncation error at this stage...
!When we start going to real space might want to transform on to the full density grid to 
!do our GW products.

    allocate (h_diag (npwx, 1))
    allocate (etc(nbnd, nkstot))

    ci = (0.0d0, 1.0d0)
!Convert freq array generated in freqbins into rydbergs.
    w_ryd(:) = wgreen(:)/RYTOEV

    CALL start_clock('greenlinsys')
    where_rec='no_recover'


   if (nksq.gt.1) rewind (unit = iunigk)

!Loop over q in the IBZ_{k}
   do iq = 1, nksq 
      if (lgamma) then
           ikq = iq
          else
           ikq = 2*iq
      endif
      if (nksq.gt.1) then
          read (iunigk, err = 100, iostat = ios) npw, igk
 100      call errore ('green_linsys', 'reading igk', abs (ios) )
      endif
      if(lgamma) npwq=npw 
      if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
 200       call errore ('green_linsys', 'reading igkq', abs (ios) )
      endif
!Need a loop to find all plane waves below ecutsco when igkq takes us outside of this sphere.
!igkq_tmp is gamma centered index up to ngmsco,
!igkq_ig  is the linear index for looping up to npwq.
!need to loop over...
      counter = 0
      igkq_tmp(:) = 0
      igkq_ig(:)  = 0 
      do ig = 1, npwx
         if((igkq(ig).le.ngmsco).and.((igkq(ig)).gt.0)) then
             counter = counter + 1
            !index in total G grid.
             igkq_tmp (counter) = igkq(ig)
            !index for loops 
             igkq_ig  (counter) = ig
         endif
      enddo
!Difference in parallelization routine. Instead of parallelizing over the usual list of G-vectors as in the straight
!forward pilot implementation I need to first generate the list of igkq's within my correlation cutoff
!this gives the number of vectors that requires parallelizing over. Then I split the work (up to counter) between the
!nodes as with the coulomb i.e. igstart and igstop.

#ifdef __PARA
      npool = nproc / nproc_pool
      write(stdout,'("npool", i4, i5)') npool, counter
      if (npool.gt.1) then
      ! number of g-vec per pool and reminder
        ngpool = counter / npool
        ngr = counter - ngpool * npool
      ! the remainder goes to the first ngr pools
        if ( my_pool_id < ngr ) ngpool = ngpool + 1
        igs = ngpool * my_pool_id + 1
        if ( my_pool_id >= ngr ) igs = igs + ngr
      ! the index of the first and the last g vec in this pool
        igstart = igs
        igstop = igs - 1 + ngpool
        write (stdout,'(/4x,"Max n. of G vecs in Green_linsys per pool = ",i5)') igstop-igstart+1
      else
#endif
       igstart = 1
       igstop = counter
#ifdef __PARA
      endif
#endif
! Now the G-vecs up to the correlation cutoff have been divided between pools.
! Calculates beta functions (Kleinman-Bylander projectors), with
! structure factor, for all atoms, in reciprocal space
       call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
! psi_{k+q}(r) is every ikq entry
       call davcio (evq, lrwfc, iuwfc, ikq, - 1)
       do ig = 1, npwq
          g2kin (ig) = ((xk (1,ikq) + g (1, igkq(ig) ) ) **2 + &
                        (xk (2,ikq) + g (2, igkq(ig) ) ) **2 + &
                        (xk (3,ikq) + g (3, igkq(ig) ) ) **2 ) * tpiba2
       enddo
!The 1/|k+G|**2 preconditioner
!       h_diag = 0.d0
!       do ibnd = 1, nbnd_occ (ikq)
!          do ig = 1, npwq
!             h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ikq))
!          enddo
!          IF (noncolin) THEN
!             do ig = 1, npwq
!                h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ikq))
!             enddo
!          END IF
!       enddo
!TPA preconditioner:
!      do ibnd = 1, nbnd_occ
!         do ig = 1, npwq
!            x = g2kin(ig)/(eprec(ikq,ibnd))
!            h_diag(ig,ibnd) = (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0) &
!                            / (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0+16.d0*x**4.d0)
!         enddo
!      enddo
WRITE(6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xk(:,ikq), et(:,ikq)*RYTOEV
WRITE(600+mpime, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xk(:,ikq), et(:,ikq)*RYTOEV

!Set eprecloc to kinetic energy matrix element at top of valence band for the k point.
!NOTA BENE nbnd_occ hardwired to four this should be variable!!
   aux1=(0.d0,0.d0)
   DO ig = 1, npwq
      aux1 (ig) = g2kin (ig) * evq (ig,4)
   END DO
   eprecloc = 1.35d0*zdotc(npwx*npol, evq(1,4), 1, aux1(1),1)
   WRITE(6,'("<|(k+G)^2|>", 1f12.7)') eprecloc*RYTOEV

DO iw = 1, nwgreen
      green =(0.0d0, 0.0d0)
      h_diag(:,:) = 0.d0
      do ibnd = 1, 1
!For all elements up to <Ekin> use standard TPA:
         if (w_ryd(iw).lt.0.0d0) then
            do ig = 1, npwq
!The preconditioner needs to be real and symmetric to decompose as E^T*E
               x = (g2kin(ig)-w_ryd(iw))/eprecloc
               h_diag(ig,ibnd) =  (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0) &
                                 /(27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0+16.d0*x**4.d0)
            enddo
         else 
         !Really choosy preconditioner.
            do ig = 1, npwq
               if(g2kin(ig).gt.w_ryd(iw)) then
                  x = (g2kin(ig) - w_ryd(iw))/eprecloc
                  h_diag(ig,ibnd) =  (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0) &
                                    /(27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0+16.d0*x**4.d0)
               else
                   h_diag(ig,ibnd) = 1.0d0
               endif
            enddo
         endif
      enddo

!for printing out preconditioner.
!      if(mpime.eq.0) then
!         write(800,*)wgreen(iw)
!         do ig = 1, npwq
!            write(800,'(i4,1f12.7)')ig, h_diag(ig,1)
!         enddo
!      endif

!PARA case:
 do ig = igstart, igstop
         rhs(:,:)  = (0.0d0, 0.0d0)
         rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
         gr_A(:,:) = (0.0d0, 0.0d0) 
         lter = 0
         etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
         cw = CMPLX( w_ryd(iw), eta, kind=DP)
! Doing Linear System with Wavefunction cutoff (full density). 
         call  cbcg_solve_fix(cch_psi_all_green, cg_psi, etc(1,ikq), rhs, gr_A, h_diag, &
                              npwx, npwq, tr2_green, ikq, lter, conv_root, anorm, 1, npol, cw, .true.)
       !if(.not.conv_root) write(600+mpime,'("root not converged. ")') 
       !if(.not.conv_root) write(600+mpime,'(f15.10, i4)')anorm, lter 
       !if(.not.conv_root) write(600+mpime,'(f15.10)')    wgreen(iw)
       !Brutal conditions.
       !if(anorm.gt.0.50d0) write(600+mpime,'(f15.10, i4)')anorm, lter 
       !if(anorm.gt.0.50d0) write(600+mpime,'(f15.10)')    wgreen(iw)
       !if(anorm.gt.1.0d0)  write(600+mpime,'("root not converged. ")') 
       !instead of zeroing green's function why not try  G(G,G';\omega) = - Delta(G,G')/((k+g)**2 - \omega)
       if(anorm.gt.1.0d0) gr_A = (0.0d0, 0.0d0) 
       !if(anorm.gt.1.0d0) gr_A(igkq_ig(ig), 1) = (-1.0d0, 0.0d0)/(DCMPLX(g2kin(igkq_ig(ig)), 0.d0)-DCMPLX(w_ryd(iw),eta))
       !alternatively:
       if(anorm.gt.1.0d0) gr_A(igkq_ig(ig), 1) = (1.0d0, 0.0d0)/(DCMPLX(w_ryd(iw),eta))

       do igp = 1, counter
          green (igkq_tmp(ig), igkq_tmp(igp)) = green (igkq_tmp(ig), igkq_tmp(igp)) + gr_A(igkq_ig(igp),1)
          !wonder if this is the problem with LiCl and SiC
          !green (igkq_tmp(ig), igkq_tmp(igp)) = green (igkq_tmp(ig), igkq_tmp(igp)) + conjg(gr_A(igkq_ig(igp),1))
       enddo
 enddo !ig
!Green's Fxn Non-analytic Component:
!do ig = 1, counter
!PARA
    do ig = igstart, igstop
       do igp = 1, counter       
!should be nbnd_occ:
        do ibnd = 1, 4
           x = w_ryd(iw) - et(ibnd, ikq)
           dirac = eta / pi / (x**2.d0 + eta**2.d0)
          !Green should now be indexed (igkq_tmp(ig), igkq_tmp(igp)) according to the
          !large G-grid which extends out to 4*ecutwfc. Keep this in mind when doing
          !ffts and taking matrix elements (especially matrix elements in G space!). 
           green(igkq_tmp(ig), igkq_tmp(igp)) =  green(igkq_tmp(ig), igkq_tmp(igp)) + &
                                                 tpi*ci*(evq(igkq_ig(ig), ibnd))    * &
                                                 conjg(evq(igkq_ig(igp), ibnd))     * dirac
        enddo 
!G=000, G'=000
!if((igkq_tmp(ig).eq.1).and.(igkq_tmp(igp).eq.1)) write(400+mpime, '(3f15.10)') wgreen(iw), green(igkq_tmp(ig), igkq_tmp(igp))
!G=111, G'= 111
!if((igkq_tmp(ig).eq.7).and.(igkq_tmp(igp).eq.7)) write(400+mpime, '(3f15.10)') wgreen(iw), green(igkq_tmp(ig), igkq_tmp(igp))
!G=111, G'= 200
!if((igkq_tmp(ig).eq.7).and.(igkq_tmp(igp).eq.11)) write(400+mpime, '(3f15.10)') wgreen(iw), green(igkq_tmp(ig), igkq_tmp(igp))
!G=200, G'= 020 
!if((igkq_tmp(ig).eq.11).and.(igkq_tmp(igp).eq.13)) write(400+mpime,'(3f15.10)')  wgreen(iw), green(igkq_tmp(ig), igkq_tmp(igp))
       enddo
    enddo 

!Collect G vectors across processors and then write the full green's function to file. 
#ifdef __PARA
    CALL mp_barrier(inter_pool_comm)
!Collect all elements of green's matrix from different
!processors.
    CALL mp_sum (green, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    if(ionode) then
#endif
!HL Original:
!  rec0 = (iw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
!  only want first record. 
!  rec0 = (iw-1) * 1 * nksq + (iq-1) + 1
!  write ( iungreen, rec = rec0, iostat = ios) green
     rec0 = (iw-1) * 1 * nksq + (iq-1) + 1
     CALL davcio(green, lrgrn, iungreen, rec0, +1, ios)
#ifdef __PARA
    endif
    CALL mp_barrier(inter_pool_comm)
!HL 1707
#endif
  ENDDO  ! iw 
ENDDO    ! iq
!Now we have Green's fxn for freq. -wsigma to wsigma, all points in irreducible BZ (depends on sigma)
!and all G vecs up to cutoff ecutsco i.e. a real space description of \Delta r ~ 1/|Gmax|.
!This should contain lots of information...
!is stop_clock causing crashes?
CALL stop_clock('greenlinsys')
RETURN
END SUBROUTINE green_linsys

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
  USE gvect,                ONLY : nrxx, g, nl, ngm, ecutwfc
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
                                   ext_recover, eta
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

  COMPLEX(DP) :: gr_A(ngm, 1), rhs(ngm , 1)
  COMPLEX(DP) :: gr_N(ngm, 1), gr(ngm, 1), ci, cw, green(ngmsig,ngmsig)

  COMPLEX(DP), ALLOCATABLE :: etc(:,:)

  INTEGER :: iw, igp, iwi
  INTEGER :: iq, ik0
  INTEGER :: rec0, n1
  REAL(DP) :: dirac, x, delta
  real(DP) :: k0mq(3) 

  !HL eta should probably be a constant user defined parameter as well...

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

!@10TION
!HL need a threshold here for the linear system solver. This could also go in the punch card
!with some default at a later date. 

        REAL(DP) :: tr_cgsolve = 1.0d-5

        allocate (h_diag ( ngm, nbnd))
        allocate (etc(nbnd, nkstot))
  
        ci = (0.0d0, 1.0d0)

!Convert freq array generated in freqbins into rydbergs.
        w_ryd(:) = wgreen(:)/RYTOEV
        CALL start_clock('greenlinsys')

! WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,ikq)*RYTOEV
! generally the eigenvalues/for each kpoint we want will be stored somewhere in the array xk
! i.e. if we want to look at a specific k-point we need to know where it sits in that array.
! for now I am only interested in k = Gamma, and therefore k-q will always be in the second spot in all these arrays.
! ikks(ik) = 2*ik. ikqs(ik) = 2*ik. ik=1=Gamma.
!        ikk = ikks(ik)
!        ikq = ikqs(ik)

         ikq = 2
         where_rec='no_recover'

!HL for q = Gamma the k+q list coincides with the klist

        if (lgamma) write(6, '("lgamma=.true.")')

        if (nksq.gt.1) rewind (unit = iunigk)

        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
 100        call errore ('green_linsys', 'reading igk', abs (ios) )
        endif


        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
 200        call errore ('green_linsys', 'reading igkq', abs (ios) )
        endif

!      CALL gk_sort( xk(1,ik0), ngm, g, ( ecutwfc / tpiba2 ), ngm - 1, igk, g2kin )
!      CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ), &
!                      ngm , igkq, g2kin )
!      write(6,*) igkq

      if (lgamma) npwq = npw 

      !  write(6,*) npwq
      !  Calculates beta functions (Kleinman-Bylander projectors), with
      !  structure factor, for all atoms, in reciprocal space

       IF (lgamma) THEN
          call init_us_2 (npw, igk, xk (1, ik0), vkb)
       ELSE
          call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
          !write(6,*) vkb
       ENDIF
      
       !write(6,'(3f11.7)')xk(:,:)

       IF (lgamma) THEN
          call davcio (evq, lrwfc, iuwfc, 1, - 1)
       ELSE
          call davcio (evq, lrwfc, iuwfc, ikq, - 1)
       ENDIF

! do ibnd = 1, 4
! write(6, '("Band")')
! write(6,*) evq(:,ibnd)*conjg(evq(:,ibnd))
! enddo

       IF (lgamma) then
            do ig = 1, npw
               g2kin (ig) = ((xk (1,ik0) + g (1, igk(ig) )) **2 + &
                             (xk (2,ik0) + g (2, igk(ig) )) **2 + &
                             (xk (3,ik0) + g (3, igk(ig) )) **2 ) * tpiba2

!               g2kin (ig) = ((xk (1,ik0) + g (1, igk(ig) )) **2 + &
!                             (xk (2,ik0) + g (2, igk(ig) )) **2 + &
!                             (xk (3,ik0) + g (3, igk(ig) )) **2 ) * tpiba2 - (6.5d0/13.605)

               !WRITE (stdout, '("g2kin  ",  3f7.4)') g2kin(ig)
            enddo
       ELSE
            do ig = 1, npwq
               g2kin (ig) = ((xk (1,ikq) + g (1, igkq(ig) ) ) **2 + &
                             (xk (2,ikq) + g (2, igkq(ig) ) ) **2 + &
                             (xk (3,ikq) + g (3, igkq(ig) ) ) **2 ) * tpiba2

!               g2kin (ig) = ((xk (1,ikq) + g (1, igkq(ig) ) ) **2 + &
!                             (xk (2,ikq) + g (2, igkq(ig) ) ) **2 + &
!                             (xk (3,ikq) + g (3, igkq(ig) ) ) **2 ) * tpiba2 - (6.5d0/13.605)
!               WRITE (stdout, '("g2kin  ",  3f7.4)') g2kin(ig)
            enddo
       ENDIF

     ! The G2KIN agrees with EMP as it should since these vectors are of identical magnitude
     ! That suggests there is a problem with the vkb or the local potential. Why is gamma working
     ! but not a generic q-point?     
     !  do ibnd = 1, nbnd_occ (ikq)
     !     ibnd = 4
     !     do ig = 1, npwq
     !        h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ikq))
     !        write(6,'(1f9.7)') h_diag(ig,ibnd)
     !     enddo
     !     IF (noncolin) THEN
     !        do ig = 1, npwq
     !           h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ikq))
     !        enddo
     !     END IF
     !  enddo
     !write(6,'(1f9.7)') eprec (4,2)  
       
    h_diag = 1.d0
    ibnd = 1

    do ig = 1, npwq

       !x = (g2kin(ig))/(eprec(nbnd_occ,ikq))
       !setting reference kinetic energy to top of the valence band.
       !SHIFT
       !x = g2kin(ig)/(eprec(2,1)-(6.2/13.605))

        x = g2kin(ig)/(eprec(2,1))

       !x = g2kin(ig)/(1.28d0)

        h_diag(ig,ibnd) = (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0) &
                   / (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0+16.d0*x**4.d0)

       !write(6,'(1f9.7)') h_diag(ig,ibnd)
    enddo

  !Set zero of energy to top of the valence band.
  !SHIFT

   ! et(:,:) = et(:,:) - (6.2/13.605)
   ! w_ryd(:) = w_ryd(:) + (6.2/13.605)


   !Setting Eigenvalues to zero. 
   !et(:,:) =  et(:,:) - (6.5/13.605)

    IF (lgamma) THEN 
        WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,ik0)*RYTOEV
    ELSE
        WRITE( 6, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xq, et(:,ikq)*RYTOEV
    ENDIF

!DO iw = 1, nwgreen 2nd look for total complots
!WRITE(6,*) 
!Note green's function is zero for almost this whole range! We need a more economical frequency
!description.

DO iw = 1, nwgreen
     green =(0.0d0, 0.0d0)
!Analytic
     do ig = 1, ngmsig !zero after this point
       rhs(:,:)  = (0.0d0, 0.0d0)
       if (lgamma) then
           !rhs(igk(ig), 1) = -(1.0d0, 0.0d0)
            rhs(ig, 1) = -(1.0d0, 0.0d0)
       else 
            rhs(ig, 1) = -(1.0d0, 0.0d0)
       endif

       gr_A(:,:) = (0.0d0, 0.0d0) 
       !conv_root = .true.
       lter = 0
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

       if(lgamma) then 
          !call  cbcg_solve_fix(cch_psi_all_green, cg_psi, etc(1,ik0), rhs, gr_A, h_diag, &
          !      npwx, npwq, tr_cgsolve, ik0, lter, conv_root, anorm, 1, npol, cw, .true.)

           call  cbcg_solve_fix(cch_psi_all_green, cg_psi, etc(1,ik0), rhs, gr_A, h_diag, &
                 npwx, npw, tr_cgsolve, ik0, lter, conv_root, anorm, 1, npol, cw, .true.)
        else
            call  cbcg_solve_fix(cch_psi_all_green, cg_psi, etc(1,ikq), rhs, gr_A, h_diag, &
                  npwx, npwq, tr_cgsolve, ikq, lter, conv_root, anorm, 1, npol, cw, .true.)
       endif
!     IF (.not.conv_root) then
!          WRITE(6,'(" Root Not Converged ")')
!          WRITE( stdout, '(/,5x," iter # ",i3," anorm :",f8.3, &
!         "gr_A(1): ",f8.3)') lter, anorm, real(gr_A(ig))
!     ENDIF


      gr = gr_A 

! do igp = 1, ngmsig

      do igp = 1, npwq
         if(((igkq(igp).lt.ngmsig).and.(igkq(ig)).lt.ngmsig).and.((igkq(igp).gt.0).and.(igkq(ig)).gt.0)) then
         green (igkq(ig), igkq(igp) ) = green (igkq(ig), igkq(igp)) + gr(igp,1)
         endif
      enddo
  enddo !ig

!NON-ANALYTIC
  do ig = 1, ngmsig
       gr_N(:,1) = (0.0d0, 0.0d0)
      do ibnd = 1, 4
          IF (lgamma) then
             x = w_ryd(iw) - et(ibnd, ik0)
          ELSE
             x = w_ryd(iw) - et(ibnd, ikq)
          ENDIF

          dirac = eta / pi / (x**2.d0 + eta**2.d0)

!    no spin factor (see above) 
!    HL HERE WE NEED TO CONSIDER the psi^*(r) psi(r') since the convention in the paper for 
!    the green's function is determined that way. 

       gr_N(1:npwq, 1) = gr_N(1:npwq,1) + tpi * ci * conjg(evq(ig,ibnd)) *  evq(1:npwq,ibnd) * dirac

      enddo !ibnd

    do igp = 1, npwq
     if(((igkq(igp).lt.ngmsig).and.(igkq(ig)).lt.ngmsig).and.((igkq(igp).gt.0) &
        .and.(igkq(ig)).gt.0)) then

     !Acccumulate the full green's function here green_A + green_N.
     !HL: note there are symmetric G-vectors here as well... for instance: q=11 7/11 2/12 are identical. 

            green(igkq(ig), igkq(igp)) =  green(igkq(ig),igkq(igp)) + gr_N(igp, 1)

     endif
    enddo !igp
  enddo !ig

!if (igkq(ig).eq.7.and.igkq(igp).eq.11)write(502,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!if (igkq(ig).eq.7.and.igkq(igp).eq.7)write(505,'(3f15.10)') wgreen(iw), green (igkq(ig),igkq(igp))
!if (igkq(ig).eq.13.and.igkq(igp).eq.13)write(508,'(3f15.10)') wgreen(iw), green (igkq(ig),igkq(igp))
!if (ig.eq.16.and.igp.eq.16)write(507,'(3f15.10)') wgreen(iw), gr_N (igp,1)
!if (ig.eq.11.and.igp.eq.13) write(509,'(3f15.10)') wgreen(iw), green (igkq(ig),igkq(igp))
!if (igkq(ig).eq.1.and.igkq(igp).eq.1)write(500,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!if (igkq(ig).eq.1.and.igkq(igp).eq.7)write(501,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!if (igkq(ig).eq.11.and.igkq(igp).eq.11)write(503,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!if (igkq(ig).eq.11.and.igkq(igp).eq.7)write(504,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!if (igkq(ig).eq.13.and.igkq(igp).eq.13)write(506,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!if (igkq(ig).eq.13.and.igkq(igp).eq.7)write(507,'(3f15.10)') wgreen(iw),green (igkq(ig), igkq(igp) )
!Collect G vectors across processors and then write the full green's function to file. 
! #ifdef __PARA
!  use poolreduce to bring together the results from each pool
!  call poolreduce ( 2 * ngms * ngms, green)
!  if (me.eq.1.and.mypool.eq.1) then
! #endif
!  rec0 = (iw-1) * nk0 * nq + (ik0-1) * nq + (iq-1) + 1

   rec0 = (iw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
   write ( iungreen, rec = rec0, iostat = ios) green

!#ifdef __PARA
!    endif
!#endif
ENDDO  !iw
!STOP

CALL stop_clock('greenlinsys')
RETURN
END SUBROUTINE green_linsys

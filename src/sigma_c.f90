SUBROUTINE sigma_c(ik0) 
! G TIMES W PRODUCT
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE io_files,      ONLY : iunigk, prefix, tmp_dir
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec, trunc_2d
  USE klist,         ONLY : wk, xk
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE eqv,           ONLY : evq, eprec
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, &
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                            wgreen, wsigma, wsigmamin, wsigmamax, &
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gvect,         ONLY : g, ngm, nl
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE modes,         ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : noinv
  USE gwsigma,       ONLY : sigma_c_st
  USE mp_global,     ONLY : mp_global_end
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE mp_pools,      ONLY : nproc_pool, me_pool, my_pool_id, inter_pool_comm, npool
  USE mp_world,      ONLY : nproc, mpime

  IMPLICIT NONE

  COMPLEX(DP)         :: ci, czero
  COMPLEX(DP)         :: phase
  COMPLEX(DP)         :: aux (sigma_c_st%dfftt%nnr)
!For running PWSCF need some variables 
  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limq, inv_q, found
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g_R (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
!v array
!COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:), barcoul_R(:,:)
  REAL(DP) :: qg2, qg, qxy, qz
!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:), greenfr(:,:)
!Integration Variable 
  COMPLEX(DP) :: cprefac
!FREQUENCY GRIDS/COUNTERS
  INTEGER  :: iwim, iw, ikq
  INTEGER  :: iw0, iw0mw, iw0pw
  REAL(DP) :: w_ryd(nwcoul)
!COUNTERS
  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER :: iqstart, iqstop, iqs, nkr
  INTEGER :: iq, ipol, iqrec
  INTEGER :: ikmq, ik0, ik, nkpool
  INTEGER :: rec0, ios
  INTEGER :: counter, ierr
  INTEGER :: inversym, screening
!SYMMETRY
  REAL(DP)              :: xq_ibk(3), xq_ibz(3)
  INTEGER               :: isym, jsym
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
!For G^NA
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ss(3,3)
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
!q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: rcut, spal, zcut
  INTEGER  :: ibnd
!CHECK FOR NAN's
  REAL(DP)     :: ar, ai
!For dirac delta fxn.
  REAL(DP)     :: dirac, x, support

!No global variables
  COMPLEX (DP), ALLOCATABLE :: sigma(:,:,:)
  COMPLEX (DP), ALLOCATABLE :: sigma_g(:,:,:)

!File related:
  character(len=256) :: tempfile, filename
!Complete file name
  integer*8 :: unf_recl

#define DIRECT_IO_FACTOR 8 
! iG(W-v)
   ALLOCATE ( scrcoul_g      (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)     )
   ALLOCATE ( scrcoul_g_R    (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)     )
   ALLOCATE ( scrcoul_pade_g (sigma_c_st%ngmt, sigma_c_st%ngmt)          )
   ALLOCATE ( greenf_g       (sigma_c_st%ngmt, sigma_c_st%ngmt)          )
!These go on the big grid...
   ALLOCATE ( scrcoul        (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr)            )
   ALLOCATE ( greenfr        (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr)            )
!Technically only need gmapsym up to sigma_c_st%ngmt or ngmgrn...
   ALLOCATE ( gmapsym  (ngm, nrot)   )
   ALLOCATE ( eigv     (ngm, nrot)   )
!This is a memory hog...
   ALLOCATE (sigma          (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr, nwsigma) )
   ALLOCATE  (z(nfs), a(nfs), u(nfs))
   w_ryd(:) = wcoul(:)/RYTOEV
   WRITE(6," ")
   WRITE(6,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') ik0, (xk(ipol, ik0), ipol=1,3)
   WRITE(6," ")
   WRITE(6,'(4x, "ngmsco, ", i4, " nwsigma, ", i4)') sigma_c_st%ngmt, nwsigma
   WRITE(6,'(4x, "nrsco, ", i4, " nfs, ", i4)') sigma_c_st%dfftt%nnr, nfs
   ci = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
   sigma(:,:,:) = (0.0d0, 0.0d0)
   CALL start_clock('sigmac')
   CALL gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
   IF(allocated(sigma)) THEN
     WRITE(6,'(4x,"Sigma allocated")')
   ELSE
     WRITE(6,'(4x,"Sigma too large!")')
     CALL mp_global_end()
     STOP
   ENDIF
   WRITE(6,'("nsym, nsymq, nsymbrav ", 3i4)'), nsym, nsymq, nrot 
!Set appropriate weights for points in the brillouin zone.
!Weights of all the k-points are in odd positions in list.
!nksq is number of k points not including k+q.
   wq(:) = 0.0d0
   DO iq = 1, nksq
      IF (lgamma) THEN
         write(6, '(" lgamma ")')
         wq(iq) = 0.5d0*wk(iq) 
      ELSE
         wq(iq) = 0.5d0*wk(2*iq-1) 
      ENDIF
   ENDDO
!Every processor needs access to the files: _gw0si.coul1 and _gw0si.green1
call mp_barrier(inter_pool_comm)

#ifdef __PARA
if(.not.ionode) then
!OPEN coulomb file (only written to by head node).
   filename = trim(prefix)//"."//"coul1"
   tempfile = trim(tmp_dir) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
   open(iuncoul, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
!OPEN green file (only written to by head node as well).
   filename = trim(prefix)//"."//"green1"
   tempfile = trim(tmp_dir) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrgrn, kind=kind(unf_recl))
   open(iungreen, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
endif
#endif

#ifdef __PARA
      npool = nproc / nproc_pool
      write(stdout,'("npool", i4, i5)') npool, nksq
      if (npool.gt.1) then
      ! number of q-vec per pool and remainder
        nkpool = nksq / npool
        nkr = nksq - nkpool * npool
      ! the remainder goes to the first nkr pools
        if ( my_pool_id < nkr ) nkpool = nkpool + 1
        iqs = nkpool * my_pool_id + 1
        if ( my_pool_id >= nkr ) iqs = iqs + nkr
      ! the index of the first and the last g vec in this pool
        iqstart = iqs
        iqstop  = iqs - 1 + nkpool
        write (stdout,'(/4x,"Max n. of Kpoint in Sigma_C per pool = ",i5)') iqstop-iqstart+1
      else
#endif
       iqstart = 1
       iqstop = nksq
#ifdef __PARA
      endif
#endif
!ONLY PROCESSORS WITH K points to process: 
IF(iqstop-iqstart+1.ne.0) THEN
  WRITE(1000+mpime, '("mpime ", i4, "  iqstart, iqstop: ", 2i5)')mpime, iqstart, iqstop
  IF (nksq.gt.1) rewind (unit = iunigk)
  DO iq = iqstart, iqstop
      IF (lgamma) THEN
          ikq = iq
      ELSE
!k+q is in even positions of list (k,k+q)
          ikq = 2*iq
      ENDIF
!  q point for convolution \sum_{q \in IBZ_{k}} G_{k+q} W_{-q}
!  q = (k0 + q) - k0
!  xq_ibk(:) = xk(:,ikq) - xk_kpoints(:, ik0)
!  by actually searching for -q and calculating w of -q  
!  we ensure there is no confusion in which symmetry routine is required.
!  SYMMFIX
!  -q =  k0 - (k0 + q)
     xq_ibk(:) = xk_kpoints(:, ik0) - xk(:, ikq)

!Find which symmetry operation rotates xq_ibk back to
!The irreducible brillouin zone and which q \in IBZ it corresponds to.
!q is stored in the list x_q as positive q but all the calculations have
!been done at -q therefore we are just going to calculate \Sum G_{k+q}W_{-q}
   inv_q=.false.
   call find_q_ibz(xq_ibk, s, iqrec, isym, found_q, inv_q)

   if(lgamma) npwq=npw 

   write(6, *)  
   write(6, '("xq_IBK point")')
   write(6, '(3f11.7)') xq_ibk
   write(6, '("equivalent xq_IBZ point, symop, iqrec")')
   write(6, '(3f11.7, 2i4)') x_q(:,iqrec), isym, iqrec
   write(6,*)

   write(1000+mpime, *)  
   write(1000+mpime, '("xq_IBK point")')
   write(1000+mpime, '(3f11.7)') xq_ibk
   write(1000+mpime, '("equivalent xq_IBZ point, symop, iqrec")')
   write(1000+mpime, '(3f11.7, 2i4)') x_q(:, iqrec), isym, iqrec


!Inverse Dielectric Function is Written to file at this point
!So we read that in, rotate it, and then apply the Coulomb operator.
   scrcoul_g(:,:,:)   = dcmplx(0.0d0, 0.0d0)
   scrcoul_g_R(:,:,:) = dcmplx(0.0d0, 0.0d0)

   if(modielec.and.padecont) PRINT*, "WARNING: PADECONT AND MODIELEC?"

   if(.not.modielec) CALL davcio(scrcoul_g, lrcoul, iuncoul, iqrec, -1)

!Rotate G_vectors for FFT.
!In EPW FG checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig
!I have checked that here as well and it works.
!write(6,'(11i4)')gmapsym(gmapsym(:ngmsco,isym),invs(isym))
   DO igp = 1, sigma_c_st%ngmt
      DO ig = 1, sigma_c_st%ngmt
         IF((gmapsym(ig,isym).le.sigma_c_st%ngmt).and.(gmapsym(igp,isym).le.sigma_c_st%ngmt) &
             .and.(gmapsym(ig,isym).gt.0).and.(gmapsym(ig,isym).gt.0)) then
             DO iwim = 1, nfs
!Rotating dielectric matrix with phase factor 
!for nonsymmorphic space groups normal:
               phase = eigv(ig,isym)*conjg(eigv(igp,isym))
               scrcoul_g_R(ig, igp, iwim) = scrcoul_g(gmapsym(ig,isym), gmapsym(igp,isym),iwim)*phase
             ENDDO
         ENDIF
      ENDDO
   ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Generate bare coulomb:   !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!HL using  sax cutoff
!this should be L_{z}/2
!rcut = 0.50d0*minval(sqrt(sum(at**2,1)))*alat*tpi
!rcut = rcut-rcut/50.0d0
  rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))

!Sax for silicon
!    rcut = 21.329d0
  IF(.not.modielec) THEN
    DO iw = 1, nfs
     DO ig = 1, sigma_c_st%ngmt
!SPHERICAL SCREENING
         qg2 = (g(1,ig) + xq_ibk(1))**2 + (g(2,ig) + xq_ibk(2))**2 + (g(3,ig)+xq_ibk(3))**2
         limq = (qg2.lt.eps8) 
         IF(.not.limq) THEN
             DO igp = 1, sigma_c_st%ngmt
                                          !!!!eps^{-1}(\G,\G';\omegaa  )!!!!!   v(q+G)!!!!
                scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
!!!!!!
!                scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(sqrt(e2*fpi)/(tpiba*qg), 0.0d0)
             ENDDO
         ENDIF
         qg = sqrt(qg2)
         spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
!Normal case using truncated coulomb potential.
         IF(.not.limq) THEN
            DO igp = 1, sigma_c_st%ngmt
                scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(spal, 0.0d0)
            ENDDO
         ELSE
            IF(iw.eq.1) THEN
               scrcoul_g_R(ig, igp, iw) = real(scrcoul_g_R(ig,igp,iw))
            ENDIF
            DO igp = 1, sigma_c_st%ngmt
               scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
            ENDDO
         ENDIF
     ENDDO!ig
    ENDDO!nfs
  ENDIF
!zeroing wings of W again!
    IF(iq.eq.1) THEN
        DO iw = 1, nfs
            DO ig = 2, sigma_c_st%ngmt
               scrcoul_g_R(ig,1,iw)  = dcmplx(0.0d0, 0.d0)
            ENDDO
            DO igp = 2, sigma_c_st%ngmt
               scrcoul_g_R(1,igp,iw) = dcmplx(0.0d0, 0.d0)
            ENDDO
        ENDDO
    ENDIF
    IF(.NOT.modielec) THEN
        IF(godbyneeds) THEN
        DO ig = 1, sigma_c_st%ngmt
          DO igp = 1, sigma_c_st%ngmt 
!For godby-needs plasmon pole the algebra is done assuming real frequency*i.
!that is: the calculation is done at i*wp but we pass a real number as the freq.
               DO iw = 1, nfs
                  z(iw) = dcmplx(aimag(fiu(iw)), 0.0d0)
                  u(iw) = scrcoul_g_R(ig, igp, iw)
               ENDDO
               CALL godby_needs_coeffs(nfs, z, u, a)
               DO iw = 1, nfs 
!Just overwrite scrcoul_g_R with godby-needs coefficients.
                  scrcoul_g_R (ig, igp, iw) = a(iw)
               ENDDO
          ENDDO
        ENDDO
        ELSE IF (padecont) THEN
          DO igp = 1, sigma_c_st%ngmt
           DO ig = 1, sigma_c_st%ngmt
!Pade input points on the imaginary axis
              DO iw = 1, nfs
                 z(iw) = fiu(iw)
                 u(iw) = scrcoul_g_R (ig, igp, iw)
              ENDDO
              call pade_coeff ( nfs, z, u, a)
!Overwrite scrcoul with Pade coefficients to be passed to pade_eval.
             DO iw = 1, nfs 
                scrcoul_g_R (ig, igp, iw) = a(iw)
             ENDDO
           ENDDO !enddo on ig
        ENDDO  !enddo on igp
        ELSE IF(.not.padecont.and..not.godbyneeds) THEN
                 WRITE(6,'("No screening model chosen!")')
        ENDIF
    ENDIF
!Start integration over iw +/- wcoul. 
    WRITE(6,'("Starting Frequency Integration")')
    DO iw = 1, nwcoul
       scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)

      IF(.NOT.modielec) THEN
        DO ig = 1, sigma_c_st%ngmt
           DO igp = 1, sigma_c_st%ngmt
              DO iwim = 1, nfs
                  z(iwim) = fiu(iwim)
                  a(iwim) = scrcoul_g_R (ig,igp,iwim)
              ENDDO
              pade_catch=.false.

              IF (padecont) THEN
                  call pade_eval ( nfs, z, a, dcmplx(w_ryd(iw), eta), scrcoul_pade_g (ig,igp))
              ELSE IF (godbyneeds) THEN
                  scrcoul_pade_g(ig,igp) = a(2)/(dcmplx(w_ryd(iw)**2,0.0d0)-(a(1)-(0.0d0,1.0d0)*eta)**2)
              ELSE
                   WRITE(6,'("No screening model chosen!")')
                   STOP
                   CALL mp_global_end()
              ENDIF
           ENDDO
        ENDDO
      ELSE IF (modielec) THEN
          DO ig = 1, sigma_c_st%ngmt
             CALL mod_diel(ig, xq_ibk, w_ryd(iw), scrcoul_pade_g(ig,ig), 1)
!scrcoul_pade_g(ig,ig) = mod_dielec_(xq,w_ryd,ig)*v_(q+G,truncation)
             qg2 = (g(1,ig) + xq_ibk(1))**2 + (g(2,ig) + xq_ibk(2))**2 + (g(3,ig)+xq_ibk(3))**2
             limq = (qg2.lt.eps8) 
             IF(.not.limq) THEN
                scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
             ENDIF
             qg = sqrt(qg2)
             spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
!Normal case using truncated coulomb potential.
             IF(.not.limq) THEN
                   scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx(spal, 0.0d0)
             ELSE
                   scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
             ENDIF
          ENDDO
      ENDIF
!W(G,G';w)
         scrcoul(:,:) = czero
         CALL fft6(scrcoul_pade_g(1,1), scrcoul(1,1), sigma_c_st,1)
! Now have W(r,r';omega')
! simpson quadrature: int_w1^wN f(w)dw = deltaw * 
! [ 1/3 ( f1 + fN ) + 4/3 sum_even f_even + 2/3 sum_odd f_odd ]
         cprefac = (deltaw/RYTOEV)*wq(iq)*(0.0d0, 1.0d0)/tpi
         IF ( iw/2*2.eq.iw ) THEN
            cprefac = cprefac * 4.d0/3.d0
         ELSE
            cprefac = cprefac * 2.d0/3.d0
         ENDIF

! Presumably the fxn should be zero at the cutoff but might as well be consistent...
         IF (iw.eq.1) THEN
             cprefac = (1.0d0/3.0d0)*(deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
         ELSE IF (iw.eq.nwcoul) THEN
             cprefac = (1.0d0/3.0d0)*(deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
         ENDIF
!Calculate \Sigma(w_0)  = \int \G(w_0 \pm w') W(\w') dw'
         DO iw0 = 1, nwsigma
            iw0mw = ind_w0mw (iw0,iw)
            iw0pw = ind_w0pw (iw0,iw)

            rec0 = (iw0mw-1) * 1 * nksq + (iq-1) + 1
            CALL davcio( greenf_g, lrgrn, iungreen, rec0, -1 )
            greenfr(:,:) = czero
            CALL fft6(greenf_g(1,1), greenfr(1,1), sigma_c_st, +1)
            sigma (:,:,iw0) = sigma (:,:,iw0) + cprefac * greenfr(:,:)*scrcoul(:,:)

            rec0 = (iw0pw-1) * 1 * nksq + (iq-1) + 1
            CALL davcio(greenf_g, lrgrn, iungreen, rec0, -1)
            greenfr(:,:) = czero
            CALL fft6(greenf_g(1,1), greenfr(1,1),sigma_c_st,+1)
            sigma (:,:,iw0) = sigma (:,:,iw0) + cprefac * greenfr(:,:)*scrcoul(:,:)
        ENDDO !on iw0  
    ENDDO ! on frequency convolution over w'
  ENDDO ! end loop iqstart, iqstop 
ENDIF
    DEALLOCATE ( gmapsym          )
    DEALLOCATE ( greenfr          )
    DEALLOCATE ( greenf_g         )
    DEALLOCATE ( scrcoul          )
    DEALLOCATE ( scrcoul_pade_g   )
    DEALLOCATE ( scrcoul_g, scrcoul_g_R )
    DEALLOCATE ( z,a,u )
!#ifdef __PARA
!    CALL mp_barrier(inter_pool_comm)
!    CALL mp_sum(sigma, inter_pool_comm)
!    CALL mp_barrier(inter_pool_comm)
!#endif __PARA
IF (ionode) THEN
  ALLOCATE ( sigma_g (sigma_c_st%ngmt, sigma_c_st%ngmt, nwsigma))
  IF(allocated(sigma_g)) THEN
     WRITE(6,'(4x,"Sigma_g allocated")')
  ELSE
     WRITE(6,'(4x,"Sigma_g too large!")')
     CALL mp_global_end()
     STOP
  ENDIF
  WRITE(6,'(4x,"Sigma in G-Space")')
  sigma_g = (0.0d0,0.0d0)
  DO iw = 1, nwsigma
     CALL fft6(sigma_g(1,1,iw), sigma(1,1,iw), sigma_c_st, -1)
  ENDDO
!Now write Sigma in G space to file. 
  CALL davcio (sigma_g, lrsigma, iunsigma, ik0, 1)
  WRITE(6,'(4x,"Sigma Written to File")')
  CALL stop_clock('sigmac')
  DEALLOCATE ( sigma_g  )
ENDIF !ionode
call mp_barrier(inter_pool_comm)
DEALLOCATE ( sigma    )
RETURN
END SUBROUTINE sigma_c


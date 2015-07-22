SUBROUTINE sym_sigma_c_im(ik0) 
!G TIMES W PRODUCT
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE io_files,      ONLY : iunigk, prefix, tmp_dir
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec, trunc_2d, tmp_dir_coul
  USE klist,         ONLY : wk, xk, nkstot
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE eqv,           ONLY : evq, eprec
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, &
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                            wgreen, wsigma, wsigmamin, wsigmamax, &
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw, &
                            w0pmw, wgtcoul
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gvect,         ONLY : g, ngm, nl
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at,bg
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE modes,         ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : noinv
  USE gwsigma,       ONLY : sigma_c_st
  USE mp_global,     ONLY : mp_global_end
  USE mp_world,      ONLY : nproc, mpime
  USE mp_images,     ONLY : nimage, my_image_id, intra_image_comm,   &
                            me_image, nproc_image, inter_image_comm
  USE mp,            ONLY : mp_sum, mp_barrier
  USE mp_pools,      ONLY : inter_pool_comm

  IMPLICIT NONE

  COMPLEX(DP)         :: ci, czero
  COMPLEX(DP)         :: phase
  COMPLEX(DP)         :: aux (sigma_c_st%dfftt%nnr)
!Sigma arrays
  COMPLEX (DP), ALLOCATABLE :: sigma(:,:,:)
  COMPLEX (DP), ALLOCATABLE :: sigma_g(:,:,:)
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g_R (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:,:), greenfr(:,:)
  COMPLEX(DP) :: cprefac
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
!v array
!COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:), barcoul_R(:,:)
  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: w_ryd(nwcoul)
  REAL(DP) :: xq_ibk(3), xq_ibz(3)
!q-vector of coulomb potential:
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: rcut, spal
!CHECK FOR NAN's
  REAL(DP)     :: ar, ai
!For dirac delta fxn.
  REAL(DP)     :: dirac, x, support, zcut
!FREQUENCY GRIDS/COUNTERS
  INTEGER  :: iwim, iw, ikq
  INTEGER  :: iw0, iw0mw, iw0pw
!COUNTERS
  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER :: iqs, nkr
  INTEGER :: iq, ipol, iqrec
  INTEGER :: ikmq, ik0, ik, nkpool
  INTEGER :: rec0, ios
  INTEGER :: counter, ierr
  INTEGER :: inversym, screening
!SYMMETRY
  INTEGER               :: isym, jsym, isymop, nig0
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
!For G^NA
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ss(3,3)
  INTEGER     :: ibnd
  INTEGER     :: iw0start, iw0stop
!For running PWSCF need some variables 
  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limq, inv_q, found
!File related:
  character(len=256) :: tempfile, filename
!Complete file name
  integer*8 :: unf_recl
  REAL (DP)   :: xk1(3), aq(3)
  INTEGER     :: iq1
  REAL(DP)    :: sxq(3,48), xqs(3,48)
  INTEGER     :: imq, isq(48), nqstar, nkpts
  INTEGER     :: nsq(48), i, ikstar
  REAL(DP)    :: wgt(nsym), xk_un(3,nsym)
  INTEGER     :: ixk1, numxk1

  logical, external :: eqvect

#define DIRECT_IO_FACTOR 8 
! iG(W-v)
   ALLOCATE ( scrcoul_g       (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)    )
   ALLOCATE ( scrcoul_g_R     (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)    )
   ALLOCATE ( scrcoul_pade_g  (sigma_c_st%ngmt, sigma_c_st%ngmt)         )
   ALLOCATE ( greenf_g        (sigma_c_st%ngmt, sigma_c_st%ngmt, 2*nwcoul) )
!These go on the big grid...
   ALLOCATE ( scrcoul        (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr))
   ALLOCATE ( greenfr        (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr))
!Technically only need gmapsym up to sigma_c_st%ngmt or ngmgrn...
   ALLOCATE ( gmapsym  (ngm, nrot)   )
   ALLOCATE ( eigv     (ngm, nrot)   )
!This is a memory hog...
   ALLOCATE (sigma  (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr, nwsigma))
   ALLOCATE  (z(nfs), a(nfs), u(nfs))
   w_ryd(:) = wcoul(:)/RYTOEV

   WRITE(6," ")
   WRITE(6,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
   WRITE(6," ")
   WRITE(6,'(4x, "ngmsco, ", i4, " nwsigma, ", i4)') sigma_c_st%ngmt, nwsigma
   WRITE(6,'(4x, "nrsco, ", i4, " nfs, ", i4)') sigma_c_st%dfftt%nnr, nfs
   zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
   !WRITE(6,'("zcut ", f12.7)'), zcut
   ci = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
   sigma(:,:,:) = (0.0d0, 0.0d0)
   CALL start_clock('sigmac')
   CALL gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
!   if(mpime.eq.1) then
!    do isym=1,nsym
!      do ig =1,sigma_c_st%ngmt
!         write(3000+mpime,'("ig, gmapsym", i5, i5)') ig, gmapsym(ig,isym)
!      enddo
!    enddo 
!   endif
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
!Every processor needs access to the files: _gw0si.coul1 and _gw0si.green1
   call mp_barrier(inter_image_comm)
#ifdef __PARA
!OPEN coulomb file (only written to by head node).
   filename = trim(prefix)//"."//"coul1"
   tempfile = trim(tmp_dir_coul) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
   open(iuncoul, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
#endif
   CALL para_img(nwsigma, iw0start, iw0stop)
   WRITE(6, '(5x, "nwsigma ",i4, " iw0start ", i4, " iw0stop ", i4)') nwsigma, iw0start, iw0stop
   write(6,*) x_q(:,:)
!ONLY PROCESSORS WITH K points to process: 
  IF (nksq.gt.1) rewind (unit = iunigk)
WRITE(6,'("Starting Frequency Integration")')
! kpoints split between pools
!CALL star_q(xk_kpoints(:,ik0), at, bg, nsym, s, invs, nkpts, sxq, isq, nsq, imq, .false. )
!WRITE( 1000+mpime, * )
!WRITE( 1000+mpime, '(5x,a,i4)') 'List of ks in the star:', nkpts
!WRITE( 1000+mpime, '(7x,i4,i4,i4,3f14.9)') (iq1, nsq(iq1), isq(iq1), (sxq(i,iq1), i=1,3), iq1=1,nkpts)
!WRITE( 1000+mpime, * )
!WRITE( 1000+mpime, '(5x,a,i4)') 'List of ks on this processor:', nkstot
!WRITE( 1000+mpime, '(7x,3f14.9)')((xk(i,iq1), i=1,3), iq1=1,nkstot)
!WRITE( 1000+mpime, '(7x,3f14.9)')((x_q(i,iq1), i=1,3), iq1=1,nkstot)
DO iq = 1, nksq
   scrcoul_g(:,:,:)   = dcmplx(0.0d0, 0.0d0)
   if(.not.modielec) CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, -1)
   xq(:) = x_q(:,iq)
   cprefac = 0.50d0*wk(iq)*dcmplx(-1.0d0, 0.0d0)/tpi
   CALL coulpade(scrcoul_g(1,1,1), xq(1))
   nsq(:) = 0
!   CALL find_xk_unique(xq, xk_kpoints, xk_un(1,1), nsq(1), numxk1, wgt)
!   WRITE( 1000+mpime, '(5x,a,i4)') 'List of unique k1s:', numxk1
!   WRITE( 1000+mpime, '(7x,i4, 3f14.9)') (iq1, (xk_un(i,iq1), i=1,3), iq1=1,numxk1)
!   WRITE( 1000+mpime, '(7x,i4, f14.9)') (iq1, wgt(iq1), iq1=1,numxk1)
!   WRITE( 1000+mpime, '(7x,i4, i4)') (iq1, nsq(iq1), iq1=1,numxk1)
!   need the star of q based on the symmetry group of k!
   DO iw0 = iw0start, iw0stop
      write(1000+mpime, '("iw0", i4)') iw0
      DO isymop = 1, nsym
        CALL rotate(xq, aq, s, nsym, invs(isymop))
        xk1 = xk_kpoints(:,ik0) - aq(:)
        nig0 = 1
        inv_q=.false.
        call find_qG_ibz(xk1, s, iqrec, isym, nig0, found_q, inv_q)
        if(inv_q) write(1000+mpime, '("Need to use time reversal")')

        write(1000+mpime, '("xq point, iq, nksq")')
        write(1000+mpime, '(3f11.7, 2i4)') xq(:), iq, isymop
        !write(1000+mpime, '("xk1 point, isym, nig0")')
        !write(1000+mpime, '(3f11.7, 2i4)') xk_un(:, ixk1), numxk1, nig0
        write(1000+mpime, '("xk point, isym, iqrec")')
        write(1000+mpime, '(3f11.7, 2i4)') x_q(:, iqrec), isym, iqrec
        write(1000+mpime, *)
        
        CALL green_linsys_shift_im(greenf_g(1,1,1), iw0, iqrec, 2*nwcoul)

        if(modielec.and.padecont) PRINT*, "WARNING: PADECONT AND MODIELEC?"
!Start integration over iw +/- wcoul.
!Rotate W and initialize necessary quantities for pade_continuation or godby needs.
!Calculate seed system: G(G,G';w=0).
        DO iw = 1, nwcoul
           CALL construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), w_ryd(iw))

           scrcoul = czero
           CALL fft6_c(scrcoul_pade_g(1,1), scrcoul(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isymop, +1)

           greenfr(:,:) = czero
           CALL fft6_g(greenf_g(1,1,iw), greenfr(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isym, nig0, +1)
           sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*(wgtcoul(iw)/RYTOEV)*cprefac*greenfr(:,:)*scrcoul(:,:)
           !if(.not.inv_q) sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*(wgtcoul(iw)/RYTOEV)*cprefac*greenfr(:,:)*scrcoul(:,:)
           !if(inv_q)      sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*(wgtcoul(iw)/RYTOEV)*cprefac*conjg(greenfr(:,:))*scrcoul(:,:)

           greenfr(:,:) = czero
           CALL fft6_g(greenf_g(1,1,iw+nwcoul), greenfr(1,1), sigma_c_st, gmapsym(1,1),eigv(1,1), isym, nig0, +1)
           sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*(wgtcoul(iw)/RYTOEV)*cprefac*greenfr(:,:)*scrcoul(:,:)
           !if(.not.inv_q) sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*(wgtcoul(iw)/RYTOEV)*cprefac*greenfr(:,:)*scrcoul(:,:)
           !if(inv_q)      sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*(wgtcoul(iw)/RYTOEV)*cprefac*conjg(greenfr(:,:))*scrcoul(:,:)
        ENDDO !on iw0  
    ENDDO ! on frequency convolution over w'
  ENDDO
ENDDO

 DEALLOCATE ( gmapsym          )
 DEALLOCATE ( greenfr          )
 DEALLOCATE ( greenf_g         )
 DEALLOCATE ( scrcoul          )
 DEALLOCATE ( scrcoul_pade_g   )
 DEALLOCATE ( scrcoul_g, scrcoul_g_R )
 DEALLOCATE ( z,a,u )
#ifdef __PARA
  CALL mp_barrier(inter_image_comm)
  !CALL mp_sum(sigma, inter_pool_comm)
  CALL mp_sum(sigma, inter_image_comm)
  CALL mp_barrier(inter_image_comm)
#endif __PARA
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
  CALL mp_barrier(inter_image_comm)
  DEALLOCATE ( sigma  )
RETURN
END SUBROUTINE sym_sigma_c_im

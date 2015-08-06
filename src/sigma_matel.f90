SUBROUTINE sigma_matel (ik0)
  USE io_global,            ONLY : stdout, ionode_id, ionode, meta_ionode
  USE io_files,             ONLY : prefix, iunigk, wfc_dir
  USE buffers,              ONLY : get_buffer, close_buffer
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecs,                ONLY : nls
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma, wsigmin, wsigmax, deltaws
  USE klist,                ONLY : xk, wk, nkstot, nks
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma,lrsex, iunsex, iunsigext, lrsigext
  USE control_gw,           ONLY : nbnd_occ, lgamma, do_imag, do_serial, do_sigma_exxG, tmp_dir_coul
  USE wavefunctions_module, ONLY : evc
  USE gwsigma,              ONLY : sigma_x_st, sigma_c_st, nbnd_sig, corr_conv
  USE disp,                 ONLY : xk_kpoints, x_q, nqs
  USE noncollin_module,     ONLY : nspin_mag
  USE eqv,                  ONLY : dmuxc, eprec
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE cell_base,            ONLY : omega, tpiba2, at, bg, alat
  USE buiol,                ONLY : buiol_check_unit
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE fft_custom,           ONLY : fft_cus, set_custom_grid, ggent, gvec_init
  USE mp_pools,             ONLY : inter_pool_comm, npool, kunit, my_pool_id
  USE mp_world,      ONLY : nproc, mpime
  USE save_gw,       ONLY : tmp_dir_save
  USE mp_images,     ONLY : nimage, my_image_id, intra_image_comm
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum

IMPLICIT NONE
  INTEGER                   ::   ig, igp, nw, iw, ibnd, jbnd, ios, ipol, ik0, ir,irp, counter
  REAL(DP)                  ::   w_ryd(nwsigma)
  REAL(DP)                  ::   resig_diag(nwsigma,nbnd_sig), imsig_diag(nwsigma,nbnd_sig),&
                                 et_qp(nbnd_sig), a_diag(nwsigma,nbnd_sig)
  REAL(DP)                  ::   dresig_diag(nwsigma,nbnd_sig), vxc_tr, vxc_diag(nbnd_sig),&
                                 sigma_ex_tr, sigma_ex_diag(nbnd_sig)
  REAL(DP)                  ::   resig_diag_tr(nwsigma), imsig_diag_tr(nwsigma), a_diag_tr(nwsigma),&
                                 et_qp_tr, z_tr, z(nbnd_sig)
  REAL(DP)                  ::   one, zcut
  COMPLEX(DP)               ::   czero, temp
  COMPLEX(DP)               ::   aux(sigma_x_st%ngmt), psic(dfftp%nnr), vpsi(ngm), auxsco(sigma_c_st%ngmt)
  COMPLEX(DP)               ::   ZDOTC, sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                                 sigma_band_ex(nbnd_sig, nbnd_sig), vxc(nbnd_sig,nbnd_sig)
  LOGICAL                   ::   do_band, do_iq, setup_pw, exst, single_line
  INTEGER                   ::   iq, ikq, ikq_head
  COMPLEX(DP), ALLOCATABLE  ::   sigma(:,:,:)
  COMPLEX(DP), ALLOCATABLE  ::   evc_tmp_j(:), evc_tmp_i(:)
  INTEGER, ALLOCATABLE      ::   igkq_ig(:) 
  INTEGER, ALLOCATABLE      ::   igkq_tmp(:) 
!for analytic continuation of selfenergy:
  INTEGER                   ::   nwsigwin, ierr, ng
  REAL(DP), ALLOCATABLE     ::   wsigwin(:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_g_ex(:,:)
  REAL(DP) :: xk_collect(nkstot), wk_collect(nkstot)
!For VXC matrix elements:
  REAL(DP) :: vtxc, etxc, ehart, eth, charge
!arbitrary cutoff
  INTEGER :: sigma_c_ngm
  real(DP), parameter :: eps=1.e-5_dp
  REAL(DP) :: zero(3)
  logical, external :: eqvect
  logical :: found_k
  character (len=256) :: poolnum
  integer*8 :: unf_recl
  INTEGER   :: iunwfc1
  INTEGER   :: kpoolid(nkstot), iqrec1(nkstot)
  INTEGER :: nbase, nksloc, rest, mypoolid
  character (len=256) :: form_str
  character(len=256) :: tempfile, filename

#define DIRECT_IO_FACTOR 8


  ALLOCATE (igkq_tmp(npwx))
  ALLOCATE (igkq_ig(npwx))

  one   = 1.0d0 
  czero = (0.0d0, 0.0d0)
  w_ryd = wsigma/RYTOEV
  nbnd = nbnd_sig 
  zero(:) = 0.d0
  lgamma=.true.
  ikq = 1
  found_k = .false.
  do iq = 1, nqs
     found_k  = (abs(xk_kpoints(1,ik0) - x_q(1,iq)).le.eps).and. &
                (abs(xk_kpoints(2,ik0) - x_q(2,iq)).le.eps).and. &
                (abs(xk_kpoints(3,ik0) - x_q(3,iq)).le.eps)
     if (found_k) then
        ikq_head = iq
        exit
     endif
  enddo

  do iq = 1, nks
     found_k  = (abs(xk_kpoints(1,ik0) - xk(1,iq)).le.eps).and. &
                (abs(xk_kpoints(2,ik0) - xk(2,iq)).le.eps).and. &
                (abs(xk_kpoints(3,ik0) - xk(3,iq)).le.eps)
     if (found_k) then
        ikq = iq
        exit
     endif
  enddo

  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ik0, (xk_kpoints(ipol,ik0) , ipol = 1, 3)
  kpoolid = 0
  iqrec1  = 0

!All pools need access to sigma file now:
  filename = trim(prefix)//"."//"sigma1"
  tempfile = trim(tmp_dir_coul) // trim(filename)
  unf_recl = DIRECT_IO_FACTOR * int(lrsigma, kind=kind(unf_recl))
  open(iunsigma, file = trim(adjustl(tempfile)), iostat = ios, &
  form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)

  filename = trim(prefix)//"."//"sigma_ex1"
  tempfile = trim(tmp_dir_coul) // trim(filename)
  unf_recl = DIRECT_IO_FACTOR * int(lrsex, kind=kind(unf_recl))
  open(iunsex, file = trim(adjustl(tempfile)), iostat = ios, &
  form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)

!ONLY THE POOL WITH THIS KPOINT CALCULATES THE CORRECT MATRIX ELEMENT.
  vxc(:,:) = czero
  sigma_band_ex (:, :) = czero
  sigma_band_c (:,:,:) = czero

  IF (found_k) THEN
      WRITE(1000+mpime,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)
      IF (nksq.gt.1) then
          CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                        npw, igk, g2kin )
      ENDIF
      IF(lgamma) npwq = npw
      call get_buffer (evc, lrwfc, iuwfc, ikq)
      zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
! generate v_xc(r) in real space:
      v%of_r(:,:) = (0.0d0)
      CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
      WRITE(1000+mpime, '("Taking Matels.")')
      WRITE(1000+mpime, '("Taking NPWQ.", i4)') npwq
      WRITE(1000+mpime, '("my_image_id", i4)') my_image_id
      DO jbnd = 1, nbnd_sig
         psic = czero
         DO ig = 1, npwq
            psic ( nls (igkq(ig)) ) = evc(ig, jbnd)
         ENDDO
!Need to do this fft according to igkq arrays and switching between serial/parallel routines. 
         CALL invfft ('Wave', psic(:), dffts)
         DO ir = 1, dfftp%nnr
            psic (ir) = psic(ir) * v%of_r (ir,1)
         ENDDO
         CALL fwfft ('Wave', psic(:), dffts)
         DO ig = 1, npwq
            vpsi(ig) = psic(nls(igkq(ig)))
         ENDDO
         DO ibnd = 1, nbnd_sig
            vxc(ibnd,jbnd) = ZDOTC (npwq, evc (1, ibnd), 1, vpsi, 1)
         ENDDO
      ENDDO
      WRITE(1000+mpime, '(4x,"VXC (eV)")')
      WRITE(1000+mpime, '(8(1x,f7.3))') real(vxc(:,:))*RYTOEV
      WRITE(1000+mpime, '("Max number Plane Waves WFC ", i4)') npwx
      WRITE(1000+mpime, '("Sigma_Ex Matrix Element")') 
      counter  = 0
      igkq_tmp(:) = 0
      igkq_ig(:)  = 0
      DO ig = 1, npwq
         IF((igkq(ig).le.sigma_x_st%ngmt).and.((igkq(ig)).gt.0)) then
          counter = counter + 1
          igkq_tmp (counter) = igkq(ig)
          igkq_ig  (counter) = ig
         ENDIF
      ENDDO
      ALLOCATE (sigma_g_ex (sigma_x_st%ngmt, sigma_x_st%ngmt))
      ALLOCATE (evc_tmp_i  (sigma_x_st%ngmt))
      ALLOCATE (evc_tmp_j  (sigma_x_st%ngmt))
      sigma_g_ex(:,:) = (0.0d0, 0.0d0)
      ios = 0 
      READ( UNIT = iunsex, REC = ik0, IOSTAT = ios ) sigma_g_ex
      !CALL davcio(sigma_g_ex, lrsex, iunsex, 1, -1)
      IF(ios /= 0) THEN
        WRITE(1000+mpime, '(5x, "Could not read Sigma_X file. Have you calculated it?")') 
      ELSE
        do ibnd = 1, nbnd_sig
           evc_tmp_i(:) = czero
          do jbnd = 1, nbnd_sig
             evc_tmp_j(:) = czero
             do ig = 1, counter
                evc_tmp_i(igkq_tmp(ig)) = evc(igkq_ig(ig), ibnd) 
             enddo
             do ig = 1, sigma_x_st%ngmt
                do igp = 1, counter
                   evc_tmp_j(igkq_tmp(igp)) = evc(igkq_ig(igp), jbnd)
                enddo
                do igp = 1, sigma_x_st%ngmt
                   sigma_band_ex (ibnd, jbnd) = sigma_band_ex (ibnd, jbnd) + evc_tmp_j (igp)*sigma_g_ex(igp,ig)*conjg(evc_tmp_i(ig))
                enddo
             enddo
          enddo
        enddo
      ENDIF
      WRITE(1000+mpime,*) 
      write(1000+mpime,'(4x,"Sigma_ex (eV)")')
      write(1000+mpime,'(8(1x,f7.3))') real(sigma_band_ex(:,:))*RYTOEV
      write(1000+mpime,*) 
      write(1000+mpime,'(8(1x,f7.3))') aimag(sigma_band_ex(:,:))*RYTOEV
      DEALLOCATE(sigma_g_ex)
      DEALLOCATE(evc_tmp_i)
      DEALLOCATE(evc_tmp_j)
!MATRIX ELEMENTS OF SIGMA_C:
      WRITE(1000+mpime,*) 
      WRITE(1000+mpime, '("Sigma_C Matrix Element")') 
      ALLOCATE (sigma(sigma_c_st%ngmt,sigma_c_st%ngmt,nwsigma)) 
      ALLOCATE (evc_tmp_i(sigma_c_st%ngmt))
      ALLOCATE (evc_tmp_j(sigma_c_st%ngmt))
      counter     = 0
      igkq_tmp(:) = 0
      igkq_ig(:)  = 0
!For convergence tests corr_conv can be set at input lower than ecutsco.
!This allows you to calculate the correlation energy at lower energy cutoffs
      IF (corr_conv.eq.sigma_c_st%ecutt) THEN
          sigma_c_ngm = sigma_c_st%ngmt
      ELSEIF(corr_conv.lt.sigma_c_st%ecutt) THEN
        DO ng = 1, ngm
           IF ( gl( igtongl (ng) ) .le. corr_conv ) sigma_c_ngm = ng
        ENDDO
      ELSE
        WRITE(6, '("Corr Conv cannot be greater than ecut_sco")')
        STOP
      ENDIF
      WRITE(1000+mpime, *)
      WRITE(1000+mpime, '(5x, "G-Vects CORR_CONV:")')
      WRITE(1000+mpime, '(5x, f6.2, i5)') corr_conv, sigma_c_ngm
      WRITE(1000+mpime, *)
      do ig = 1, npwq
         if((igkq(ig).le.sigma_c_ngm).and.((igkq(ig)).gt.0)) then
             counter = counter + 1
             igkq_tmp (counter) = igkq(ig)
             igkq_ig  (counter) = ig
         endif
      enddo
      sigma = dcmplx(0.0d0, 0.0d0)
      if(do_serial) then
         do iw = 1, nwsigma
            CALL davcio (sigma(:,:, iw), lrsigma, iunsigma, iw, -1)
         enddo
      else
!let's us avoid crash if we haven't calculated one of these things yet:
         READ( UNIT = iunsigma, REC = ik0, IOSTAT = ios ) sigma
!         CALL davcio(sigma, lrsigma, iunsigma, 1, -1)
         if(ios /= 0) then
            WRITE(1000+mpime, '("Could not read Sigma_C file. Have you calculated it?")')
            sigma_band_c (:,:,:) = czero
         else
            sigma_band_c (:,:,:) = czero
            do ibnd = 1, nbnd_sig
               evc_tmp_i(:) = czero
             do jbnd = 1, nbnd_sig
               evc_tmp_j(:) = czero
              do iw = 1, nwsigma
                 do ig = 1, counter
                    evc_tmp_i(igkq_tmp(ig)) = evc(igkq_ig(ig), ibnd)
                 enddo
                 do ig = 1, sigma_c_ngm
                    do igp = 1, counter
                       evc_tmp_j(igkq_tmp(igp)) = evc(igkq_ig(igp), jbnd)
                    enddo
                    do igp = 1, sigma_c_ngm
                       sigma_band_c (ibnd, jbnd, iw) = sigma_band_c (ibnd, jbnd, iw) + evc_tmp_j(ig)*sigma(ig,igp,iw)*conjg(evc_tmp_i(igp))
                    enddo
                 enddo
             enddo
            enddo
           enddo
           DEALLOCATE (sigma)
           DEALLOCATE (evc_tmp_i)
           DEALLOCATE (evc_tmp_j)
           WRITE (1000+mpime,'("Finished Sigma_c")')
         endif
      endif
!Need to broadcast from the current pool to all the nodes
  ENDIF!on pool with K-point
  call mp_barrier(inter_pool_comm)
  call mp_sum(vxc, inter_pool_comm)
  call mp_sum(sigma_band_c, inter_pool_comm)
  call mp_sum(sigma_band_ex, inter_pool_comm)
  call mp_barrier(inter_pool_comm)
  IF(meta_ionode) THEN
     IF(do_imag) then 
!We can set arbitrary \Sigma(\omega) energy windows with analytic continuation:
        nwsigwin  = 1 + ceiling((wsigmax - wsigmin)/deltaws)
        allocate (wsigwin(nwsigwin))
        do iw = 1, nwsigwin
            wsigwin(iw) = wsigmin + (wsigmax-wsigmin)/float(nwsigwin-1)*float(iw-1)
        enddo
        allocate (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))
!print selfenergy on the imaginary axis.
        call print_matel_im(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
!do analytic continuation and print selfenergy on the real axis.
        sigma_band_con(:,:,:) = dcmplx(0.0d0, 0.d0)
        call sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
        call print_matel(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
     ELSE
        call print_matel(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
     ENDIF
  ENDIF
  IF(allocated(sigma_band_con)) deallocate(sigma_band_con)
  IF(allocated(igkq_tmp)) deallocate(igkq_tmp)
  IF(allocated(igkq_ig)) deallocate(igkq_ig)
RETURN
END SUBROUTINE sigma_matel

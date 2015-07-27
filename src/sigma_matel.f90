SUBROUTINE sigma_matel (ik0)
  USE io_global,            ONLY : stdout, ionode_id, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE buffers,              ONLY : get_buffer
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecs,                ONLY : nls
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma, wsigmin, wsigmax, deltaws
  USE klist,                ONLY : xk, wk, nkstot
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma,lrsex, iunsex, iunsigext, lrsigext
  USE control_gw,           ONLY : nbnd_occ, lgamma, do_imag, do_serial, do_sigma_exxG
  USE wavefunctions_module, ONLY : evc
  USE gwsigma,              ONLY : sigma_x_st, sigma_c_st, nbnd_sig, corr_conv
  USE disp,                 ONLY : xk_kpoints
  USE noncollin_module,     ONLY : nspin_mag
  USE eqv,                  ONLY : dmuxc, evq, eprec
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE cell_base,            ONLY : omega, tpiba2, at, bg, alat
  USE buiol,                ONLY : buiol_check_unit
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE fft_custom,           ONLY : fft_cus, set_custom_grid, ggent, gvec_init
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
  INTEGER                   ::   iq, ikq
  COMPLEX(DP), ALLOCATABLE  ::   sigma(:,:,:)
  COMPLEX(DP), ALLOCATABLE  ::   evc_tmp_j(:), evc_tmp_i(:)
  INTEGER, ALLOCATABLE      ::   igkq_ig(:) 
  INTEGER, ALLOCATABLE      ::   igkq_tmp(:) 
!for analytic continuation of selfenergy:
  INTEGER                   ::   nwsigwin, ierr, ng
  REAL(DP), ALLOCATABLE     ::   wsigwin(:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_g_ex(:,:)
!For VXC matrix elements:
  REAL(DP) :: vtxc, etxc, ehart, eth, charge
!arbitrary cutoff
  INTEGER :: sigma_c_ngm
  real(DP), parameter :: eps=1.e-5_dp
  REAL(DP) :: zero(3)
  logical, external :: eqvect
  logical :: found_k

  ALLOCATE (igkq_tmp(npwx))
  ALLOCATE (igkq_ig(npwx))

  one   = 1.0d0 
  czero = (0.0d0, 0.0d0)
  w_ryd = wsigma/RYTOEV
  nbnd = nbnd_sig 
  zero(:) = 0.d0
  lgamma=.true.
  iq = 1
  found_k = .false.
  do while(.not.found_k)
     found_k  = (abs(xk_kpoints(1,ik0) - xk(1,iq)).le.eps).and. &
                (abs(xk_kpoints(2,ik0) - xk(2,iq)).le.eps).and. & 
                (abs(xk_kpoints(3,ik0) - xk(3,iq)).le.eps) 
     if (found_k) then
        ikq = iq
        write(6,'("K POINT FOUND ", 2i4)'),ikq,iq 
     else
        iq = iq + 1
     endif
  enddo
  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ik0, (xk_kpoints(ipol,ik0) , ipol = 1, 3)
  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)
!write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk (ipol,ikq) , ipol = 1, 3)
  IF (ionode) THEN
      IF (nksq.gt.1) then
          CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                        npw, igk, g2kin )
      ENDIF
      if(lgamma) npwq = npw
      IF (.not.lgamma.and.nksq.gt.1) then
        CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ), &
                      npwq, igkq, g2kin )
      ENDIF
  call get_buffer (evc, lrwfc, iuwfc, ikq)
  zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
  WRITE(6,'("zcut ", f12.7)'), zcut
  WRITE(6,'("NBND ", i5)'), nbnd_sig
! generate v_xc(r) in real space:
  v%of_r(:,:) = (0.0d0)
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  vxc(:,:) = (0.0d0, 0.0d0)
  WRITE(6,'("Taking Matels.")')
  WRITE(6,'("Taking NPWQ.", i4)')npwq
  do jbnd = 1, nbnd_sig
     psic = czero
     do ig = 1, npwq
        psic ( nls (igkq(ig)) ) = evc(ig, jbnd)
     enddo
!Need to do this fft according to igkq arrays and switching between serial/parallel routines. 
     CALL invfft ('Wave', psic(:), dffts)
     do ir = 1, dfftp%nnr
        psic (ir) = psic(ir) * v%of_r (ir,1)
     enddo
     CALL fwfft ('Wave', psic(:), dffts)
     do ig = 1, npwq
        vpsi(ig) = psic(nls(igkq(ig)))
     enddo
     do ibnd = 1, nbnd_sig
        vxc(ibnd,jbnd) = ZDOTC (npwq, evc (1, ibnd), 1, vpsi, 1)
     enddo
  enddo
  write(stdout,'(4x,"VXC (eV)")')
  write(stdout,'(8(1x,f7.3))') real(vxc(:,:))*RYTOEV
  WRITE(6,'("Max number Plane Waves WFC ", i4)') npwx
  WRITE(6,'("Sigma_Ex Matrix Element")') 
  counter  = 0
  igkq_tmp(:) = 0
  igkq_ig(:)  = 0
  do ig = 1, npwq
     if((igkq(ig).le.sigma_x_st%ngmt).and.((igkq(ig)).gt.0)) then
      counter = counter + 1
      igkq_tmp (counter) = igkq(ig)
      igkq_ig  (counter) = ig
    endif
  enddo
  if(do_sigma_exxG) GOTO 143
  ALLOCATE (sigma_g_ex (sigma_x_st%ngmt, sigma_x_st%ngmt))
  ALLOCATE (evc_tmp_i  (sigma_x_st%ngmt))
  ALLOCATE (evc_tmp_j  (sigma_x_st%ngmt))
  sigma_g_ex(:,:) = (0.0d0, 0.0d0)

!CALL davcio(sigma_g_ex, lrsex, iunsex, 1, -1)
  ios = 0 
  READ( UNIT = iunsex, REC = 1, IOSTAT = ios ) sigma_g_ex
  if(ios /= 0) then
    WRITE(6, '(5x, "Could not read Sigma_X file. Have you calculated it?")') 
  else
    sigma_band_ex (:, :) = czero
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
  endif

  WRITE(6,*) 
  write(stdout,'(4x,"Sigma_ex (eV)")')
  write(stdout,'(8(1x,f7.3))') real(sigma_band_ex(:,:))*RYTOEV
  write(stdout,*) 
  write(stdout,'(8(1x,f7.3))') aimag(sigma_band_ex(:,:))*RYTOEV

  DEALLOCATE(sigma_g_ex)
  DEALLOCATE(evc_tmp_i)
  DEALLOCATE(evc_tmp_j)
143 CONTINUE
!MATRIX ELEMENTS OF SIGMA_C:
  WRITE(6,*) 
  WRITE(6,'("Sigma_C Matrix Element")') 

  ALLOCATE (sigma(sigma_c_st%ngmt,sigma_c_st%ngmt,nwsigma)) 
  ALLOCATE (evc_tmp_i(sigma_c_st%ngmt))
  ALLOCATE (evc_tmp_j(sigma_c_st%ngmt))

  counter     = 0
  igkq_tmp(:) = 0
  igkq_ig(:)  = 0


  !actually want to 
  if (corr_conv.eq.sigma_c_st%ecutt) then
      sigma_c_ngm = sigma_c_st%ngmt
  elseif(corr_conv.lt.sigma_c_st%ecutt) then
    do ng = 1, ngm
       if ( gl( igtongl (ng) ) .le. corr_conv ) sigma_c_ngm = ng
    enddo
  else
    write(6, '("Corr Conv cannot be greater than ecut_sco")')
    stop
  endif

  WRITE(6,*)
  WRITE(stdout, '(5x, "G-Vects CORR_CONV:")')
  WRITE(stdout, '(5x, f6.2, i5)') corr_conv, sigma_c_ngm
  WRITE(6,*)

  do ig = 1, npwq
     !if((igkq(ig).le.sigma_c_st%ngmt).and.((igkq(ig)).gt.0)) then
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
     READ( UNIT = iunsigma, REC = 1, IOSTAT = ios ) sigma
     if(ios /= 0) then
        print*, "Couldn't read Sigma_C file. Have you calculated it?"
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
        deallocate (sigma) 
        Write(6,'("Finished Sigma_c")')
     endif
  endif
  if (do_imag) then 
!We can set arbitrary \Sigma(\omega) energy windows with analytic continuation:
    nwsigwin  = 1 + ceiling((wsigmax - wsigmin)/deltaws)
    allocate (wsigwin(nwsigwin))
    do iw = 1, nwsigwin
        wsigwin(iw) = wsigmin + (wsigmax-wsigmin)/float(nwsigwin-1)*float(iw-1)
    enddo
    allocate (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))
!print selfenergy on the imaginary axis.
    call print_matel_im(ikq, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
!do analytic continuation and print selfenergy on the real axis.
    sigma_band_con(:,:,:) = dcmplx(0.0d0, 0.d0)
    call sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
    call print_matel(ikq, vxc(1,1), sigma_band_ex(1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
  else
    call print_matel(ikq, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
  endif
ENDIF
if(allocated(sigma_band_con)) deallocate(sigma_band_con)
RETURN
END SUBROUTINE sigma_matel

SUBROUTINE sigma_matel (ik0)
  USE io_global,            ONLY : stdout, ionode_id, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE buffers,              ONLY : get_buffer
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : nls
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma
  USE klist,                ONLY : xk, wk, nkstot
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma,lrsex, iunsex, iunsigext, lrsigext
  USE control_gw,           ONLY : nbnd_occ, lgamma, do_imag, do_serial, do_sigma_exxG
  USE wavefunctions_module, ONLY : evc
  USE gwsigma,              ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE disp,                 ONLY : xk_kpoints
  USE noncollin_module,     ONLY : nspin_mag
  USE eqv,                  ONLY : dmuxc, evq, eprec
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE cell_base,            ONLY : omega, tpiba2, at, bg
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
  REAL(DP)                  ::   one
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
  REAL(DP)                  ::   wsigmax, wsigmin, deltaws
  INTEGER                   ::   nwsigwin, ierr
  REAL(DP), ALLOCATABLE     ::   wsigwin(:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_g_ex(:,:)
!For VXC matrix elements:
  REAL(DP) :: vtxc, etxc, ehart, eth, charge

  ALLOCATE (igkq_tmp(npwx))
  ALLOCATE (igkq_ig(npwx))

  one   = 1.0d0 
  czero = (0.0d0, 0.0d0)
  w_ryd = wsigma/RYTOEV
  nbnd = nbnd_sig 

  print*, lgamma
  if (lgamma) then
      ikq = ik0
  else
      ikq = 2*ik0
  endif

  write(stdout,'(/4x,"k0(",i3," ) = (",3f7.3," )")') ikq, (xk (ipol,ikq) , ipol = 1, 3)

  IF (ionode) THEN
     if (nksq.gt.1) rewind (unit = iunigk)
     if (nksq.gt.1) then
        read (iunigk, err = 100, iostat = ios) npw, igk
 100        call errore ('green_linsys', 'reading igk', abs (ios) )
     endif
  
     if(lgamma) npwq = npw

     if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npwq, igkq
 200             call errore ('green_linsys', 'reading igkq', abs (ios) )
     endif

!if just gamma then psi_{\Gamma} should be first entry in list.
   if (lgamma) then
     call get_buffer (evc, lrwfc, iuwfc, ikq)
   else
!else then psi_{\k+\gamma = \psi_{k}} should be second entry in list.
     call get_buffer (evq, lrwfc, iuwfc, ikq)
   endif

  WRITE(6,'("NBND")')
  WRITE(6,*) nbnd_sig
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

  ierr = buiol_check_unit(iunsex)
  print*, "iunsigma_x", ierr
  CALL davcio(sigma_g_ex, lrsex, iunsex, 1, -1)
  print*, sigma_g_ex(1:5,1:5)

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
            aux(igp) = sigma_g_ex (igp, ig)
         enddo
           sigma_band_ex (ibnd, jbnd) = sigma_band_ex (ibnd, jbnd) + &
           ZDOTC(sigma_x_st%ngmt, evc_tmp_j (1:sigma_x_st%ngmt), 1, aux,1)*evc_tmp_i(ig)
      enddo
    enddo
  enddo

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

  do ig = 1, npwq
     if((igkq(ig).le.sigma_c_st%ngmt).and.((igkq(ig)).gt.0)) then
       counter = counter + 1
       igkq_tmp (counter) = igkq(ig)
       igkq_ig  (counter) = ig
     endif
  enddo

  sigma = dcmplx(0.0d0, 0.0d0)

  ierr = buiol_check_unit(iunsigma)
  print*, "iunsigma", ierr

  if(do_serial) then
     do iw = 1, nwsigma
        CALL davcio (sigma(:,:, iw), lrsigma, iunsigma, iw, -1)
     enddo
  else
     READ( UNIT = iunsigma, REC = lrsigma, IOSTAT = ios ) sigma
     if(ios /= 0) then
        print*, "Couldn't read Sigma_C file. Have you calculated it?"
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
            do ig = 1, sigma_c_st%ngmt
              do igp = 1, counter
                 evc_tmp_j(igkq_tmp(igp)) = evc(igkq_ig(igp), jbnd)
              enddo
              do igp = 1, sigma_c_st%ngmt
               auxsco(igp) = sigma (ig, igp, iw)
              enddo
              sigma_band_c (ibnd, jbnd, iw) = sigma_band_c (ibnd, jbnd, iw) + &
              conjg(evc_tmp_i(ig))*ZDOTC(sigma_c_st%ngmt, conjg(evc_tmp_j (1:sigma_c_st%ngmt)), 1, auxsco, 1)
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
    wsigmin   = -30.0
    wsigmax   =  30.0
    deltaws   =   0.2
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

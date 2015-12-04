  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE sigma_matel (ik0)
  USE io_global,            ONLY : stdout, ionode_id, ionode, meta_ionode
  USE io_files,             ONLY : prefix, iunigk, wfc_dir
  USE buffers,              ONLY : get_buffer, close_buffer
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecs,                ONLY : nls
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma, wsig_wind_min, wsig_wind_max, deltaws
  USE klist,                ONLY : xk, wk, nkstot, nks
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma,lrsex, iunsex, iunsigext, lrsigext
  USE control_gw,           ONLY : nbnd_occ, lgamma, do_imag, do_serial, do_sigma_exxG, tmp_dir_coul
  USE wavefunctions_module, ONLY : evc
  USE gwsigma,              ONLY : sigma_x_st, sigma_c_st, nbnd_sig, corr_conv, exch_conv, &
                                   sigma_band_exg, gcutcorr
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
  USE mp_images,     ONLY : nimage, my_image_id, intra_image_comm,&
                            inter_image_comm
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum

IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)
  COMPLEX(DP), ALLOCATABLE  :: sigma_g_ex(:,:)
  COMPLEX(DP)               ::   czero, temp
  COMPLEX(DP)               ::   aux(sigma_x_st%ngmt), psic(dffts%nnr), vpsi(ngm), auxsco(gcutcorr)
  COMPLEX(DP)               ::   ZdoTC, sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                                 sigma_band_ex(nbnd_sig, nbnd_sig), vxc(nbnd_sig,nbnd_sig)
  COMPLEX(DP), ALLOCATABLE  ::   sigma(:,:,:)
  COMPLEX(DP), ALLOCATABLE  ::   evc_tmp_j(:), evc_tmp_i(:)
  REAL(DP), ALLOCATABLE     ::   wsigwin(:)
  REAL(DP)                  ::   w_ryd(nwsigma)
  REAL(DP)                  ::   resig_diag(nwsigma,nbnd_sig), imsig_diag(nwsigma,nbnd_sig),&
                                 et_qp(nbnd_sig), a_diag(nwsigma,nbnd_sig)
  REAL(DP)                  ::   dresig_diag(nwsigma,nbnd_sig), vxc_tr, vxc_diag(nbnd_sig),&
                                 sigma_ex_tr, sigma_ex_diag(nbnd_sig)
  REAL(DP)                  ::   resig_diag_tr(nwsigma), imsig_diag_tr(nwsigma), a_diag_tr(nwsigma),&
                                 et_qp_tr, z_tr, z(nbnd_sig)
  REAL(DP)                  ::   one, zcut
  REAL(DP)    :: xk_collect(nkstot), wk_collect(nkstot)
  REAL(DP)    :: vtxc, etxc, ehart, eth, charge
  REAL(DP)    :: zero(3)
  INTEGER, ALLOCATABLE      ::   igkq_ig(:) 
  INTEGER, ALLOCATABLE      ::   igkq_tmp(:) 
  INTEGER                   ::   iq, ikq, ikq_head
  INTEGER                   ::   ig, igp, nw, iw, ibnd, jbnd, ios, ipol, ik0, ir,irp, counter
  INTEGER                   ::   nwsigwin, ierr, ng
  INTEGER     :: sigma_c_ngm, sigma_x_ngm
  INTEGER     :: iunwfc1
  INTEGER     :: kpoolid(nkstot), iqrec1(nkstot)
  INTEGER     :: nbase, nksloc, rest, mypoolid
  LOGICAL     ::   do_band, do_iq, setup_pw, exst, single_line
  integer*8 :: unf_recl
  logical, external :: eqvect
  logical :: found_k
  character (len=256) :: poolnum
  character (len=256) :: form_str
  character(len=256) :: tempfile, filename
  real(DP), parameter :: eps=1.e-5_dp

#define DIRECT_IO_FACTOR 8


  allocate (igkq_tmp(npwx))
  allocate (igkq_ig(npwx))

  one   = 1.0d0 
  czero = (0.0d0, 0.0d0)
  w_ryd = wsigma/RYTOEV
  nbnd = nbnd_sig 
  zero(:) = 0.d0
  lgamma=.true.
  ikq = 1
  found_k = .false.

  if((xk_kpoints(1,ik0).eq.0.0).and.(xk_kpoints(2,ik0).eq.0.0).and.(xk_kpoints(3,ik0).eq.0.0))then 
     ikq_head = 1
  else
     ikq_head = 2
  endif

  ikq = ikq_head
  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ik0, (xk_kpoints(ipol,ik0) , ipol = 1, 3)
  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)
  kpoolid = 0
  iqrec1  = 0

!All pools need access to sigma file now:
  filename = trim(prefix)//"."//"sigma1"
  tempfile = trim(tmp_dir_coul) // trim(filename)
  unf_recl = DIRECT_IO_FACTOR * int(lrsigma, kind=kind(unf_recl))
  open(iunsigma, file = trim(adjustl(tempfile)), iostat = ios, &
  form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
  write(1000+mpime,*) tempfile, ios

  if(.not. do_sigma_exxG) then
     filename = trim(prefix)//"."//"sigma_ex1"
     tempfile = trim(tmp_dir_coul) // trim(filename)
     unf_recl = DIRECT_IO_FACTOR * int(lrsex, kind=kind(unf_recl))
     open(iunsex, file = trim(adjustl(tempfile)), iostat = ios, &
     form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
     write(1000+mpime,*) tempfile, ios
  endif

!ONLY THE POOL WITH THIS KPOINT CALCULATES THE CORRECT MATRIX ELEMENT.
  vxc(:,:) = czero
  sigma_band_ex (:, :) = czero
  sigma_band_c (:,:,:) = czero

  if (meta_ionode) THEN
      write(1000+mpime,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)
      CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                    npw, igk, g2kin )
      npwq = npw
      call get_buffer (evc, lrwfc, iuwfc, ikq)
      zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
! generate v_xc(r) in real space:
      v%of_r(:,:) = (0.0d0)
      CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
      write(1000+mpime, '("Taking Matels.")')
      write(1000+mpime, '("Taking NPWQ.", i4)') npwq
      write(1000+mpime, '("my_image_id", i4)') my_image_id
      do jbnd = 1, nbnd_sig
         psic = czero
         do ig = 1, npwq
            psic ( nls (igk(ig)) ) = evc(ig, jbnd)
         enddo
!Need to do this fft according to igkq arrays and switching between serial/parallel routines. 
         CALL invfft ('Wave', psic(:), dffts)
         do ir = 1, dfftp%nnr
            psic (ir) = psic(ir) * v%of_r (ir,1)
         enddo
         CALL fwfft ('Wave', psic(:), dffts)
         do ig = 1, npwq
            vpsi(ig) = psic(nls(igk(ig)))
         enddo
         do ibnd = 1, nbnd_sig
            vxc(ibnd,jbnd) = ZdoTC (npwq, evc (1, ibnd), 1, vpsi, 1)
         enddo
      enddo
      write(1000+mpime, '(4x,"VXC (eV)")')
      write(1000+mpime, '(8(1x,f7.3))') real(vxc(:,:))*RYTOEV
      write(1000+mpime, '("Max number Plane Waves WFC ", i4)') npwx
      write(1000+mpime, '("Sigma_Ex Matrix Element")') 

    if(.not.do_sigma_exxG) then
      allocate (sigma_g_ex (sigma_x_st%ngmt, sigma_x_st%ngmt))
      allocate (evc_tmp_i  (sigma_x_st%ngmt))
      allocate (evc_tmp_j  (sigma_x_st%ngmt))

      if ((exch_conv.eq.sigma_x_st%ecutt) .or. (exch_conv.eq.0.0)) THEN
          sigma_x_ngm = sigma_x_st%ngmt
      else if((exch_conv .lt. sigma_x_st%ecutt) .and. (exch_conv.gt.0.0)) THEN
        do ng = 1, ngm
           if ( gl( igtongl (ng) ) .le. (exch_conv/tpiba2)) sigma_x_ngm = ng
        enddo
      else
        write(6, '("Exch Conv must be greater than zero and less than ecut_sco")')
        STOP
      endif

      counter  = 0
      igkq_tmp(:) = 0
      igkq_ig(:)  = 0

      do ig = 1, npwq
         if((igk(ig).le.sigma_x_ngm).and.((igk(ig)).gt.0)) then
          counter = counter + 1
          igkq_tmp (counter) = igk(ig)
          igkq_ig  (counter) = ig
         endif
      enddo
      sigma_g_ex(:,:) = (0.0d0, 0.0d0)
      ios = 0 
      READ( UNIT = iunsex, REC = ik0, IOSTAT = ios ) sigma_g_ex
      if(ios /= 0) then
        write(1000+mpime, '(5x, "Could not read Sigma_X file. Have you calculated it?")') 
      else
        do ibnd = 1, nbnd_sig
           evc_tmp_i(:) = czero
          do jbnd = 1, nbnd_sig
             evc_tmp_j(:) = czero
             do ig = 1, counter
                evc_tmp_i(igkq_tmp(ig)) = evc(igkq_ig(ig), ibnd) 
             enddo
             do ig = 1, sigma_x_ngm
                do igp = 1, counter
                   evc_tmp_j(igkq_tmp(igp)) = evc(igkq_ig(igp), jbnd)
                enddo
                do igp = 1, sigma_x_ngm
                   sigma_band_ex (ibnd, jbnd) = sigma_band_ex (ibnd, jbnd) + &
&                  evc_tmp_j (igp)*sigma_g_ex(igp,ig)*conjg(evc_tmp_i(ig))
                enddo
             enddo
          enddo
        enddo
      endif
      write(1000+mpime,*) 
      write(1000+mpime,'(4x,"Sigma_ex (eV)")')
      write(1000+mpime,'(8(1x,f7.3))') real(sigma_band_ex(:,:))*RYTOEV
      write(1000+mpime,*) 
      write(1000+mpime,'(8(1x,f7.3))') aimag(sigma_band_ex(:,:))*RYTOEV
      deallocate(sigma_g_ex)
      deallocate(evc_tmp_i)
      deallocate(evc_tmp_j)
else
  do ibnd = 1, nbnd_sig
     sigma_band_ex(ibnd,ibnd) = sigma_band_exg(ibnd)
  enddo
endif
!MATRIX ELEMENTS OF SIGMA_C:
      write(1000+mpime,*) 
      write(1000+mpime, '("Sigma_C Matrix Element")') 
      allocate (sigma(gcutcorr, gcutcorr,nwsigma)) 
      allocate (evc_tmp_i(gcutcorr))
      allocate (evc_tmp_j(gcutcorr))
      counter     = 0
      igkq_tmp(:) = 0
      igkq_ig(:)  = 0
!For convergence tests corr_conv can be set at input lower than ecutsco.
!This allows you to calculate the correlation energy at lower energy cutoffs
      if (corr_conv.eq.sigma_c_st%ecutt) THEN
          sigma_c_ngm = gcutcorr
      else if(corr_conv .lt. sigma_c_st%ecutt .and. corr_conv.gt.0.0) THEN
        do ng = 1, ngm
           if ( gl( igtongl (ng) ) .le. (corr_conv/tpiba2)) sigma_c_ngm = ng
        enddo
      else
        write(6, '("Corr Conv must be greater than zero and less than ecut_sco")')
        stop
      endif

      write(1000+mpime, *)
      write(1000+mpime, '(5x, "G-Vects CORR_CONV:")')
      write(1000+mpime, '(5x, f6.2, i5)') corr_conv, sigma_c_ngm
      write(1000+mpime, *)
      do ig = 1, npwq
         if((igk(ig).le.sigma_c_ngm).and.((igk(ig)).gt.0)) then
             counter = counter + 1
             igkq_tmp (counter) = igk(ig)
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
         write(1000+mpime, *) ios
         if(ios /= 0) then
            write(1000+mpime, '("Could not read Sigma_C file. Have you calculated it?")')
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
                       sigma_band_c (ibnd, jbnd, iw) = sigma_band_c (ibnd, jbnd, iw) +  &
                                 evc_tmp_j(ig)*sigma(ig,igp,iw)*conjg(evc_tmp_i(igp))
                    enddo
                 enddo
             enddo
            enddo
           enddo
           deallocate (sigma)
           deallocate (evc_tmp_i)
           deallocate (evc_tmp_j)
           write (1000+mpime,'("Finished Sigma_c")')
         endif
      endif
!Need to broadcast from the current pool to all the nodes
  endif!on pool with K-point
call mp_barrier(inter_pool_comm)

!Now first pool should always have
!the kpoint we are looking for.
  if(meta_ionode) THEN
     if(do_imag) then 
!We can set arbitrary \Sigma(\omega) energy windows with analytic continuation:
        nwsigwin  = 1 + ceiling((wsig_wind_max - wsig_wind_min)/deltaws)
        allocate (wsigwin(nwsigwin))
        do iw = 1, nwsigwin
            wsigwin(iw) = wsig_wind_min + (wsig_wind_max-wsig_wind_min)/float(nwsigwin-1)*float(iw-1)
        enddo
        allocate (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))
!print selfenergy on the imaginary axis.
        call print_matel_im(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
        !call print_matel_im(2, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
!do analytic continuation and print selfenergy on the real axis.
        sigma_band_con(:,:,:) = dcmplx(0.0d0, 0.d0)
        call sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
        call print_matel(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
     else
        call print_matel(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
     endif
  endif
  if(allocated(sigma_band_con)) deallocate(sigma_band_con)
  if(allocated(igkq_tmp)) deallocate(igkq_tmp)
  if(allocated(igkq_ig))  deallocate(igkq_ig)
  if(allocated(sigma_band_exg)) deallocate(sigma_band_exg)

call mp_barrier(inter_pool_comm)
call mp_barrier(inter_image_comm)
return
end SUBROUTINE sigma_matel

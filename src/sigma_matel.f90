SUBROUTINE sigma_matel (ik0)
  USE io_global,            ONLY : stdout, ionode_id, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,              ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,              ONLY : fpol, fiu, nfs, nwsigma, wsigma
  USE klist,                ONLY : xk, wk, nkstot
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma,lrsex, iunsex, iunsigext, lrsigext
  USE control_gw,           ONLY : nbnd_occ, lgamma, do_imag
  USE wavefunctions_module, ONLY : evc
  USE gwsigma,              ONLY : ngmsig, nbnd_sig, sigma_g_ex, ngmsco, ngmsex
  USE disp,                 ONLY : xk_kpoints
  USE noncollin_module,     ONLY : nspin_mag
  USE eqv,                  ONLY : dmuxc
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE fft_scalar,           ONLY : cfft3ds, cfft3d
  USE fft_base,             ONLY : dffts
  USE fft_parallel,         ONLY : tg_cft3s
  USE cell_base,            ONLY : omega, tpiba2, at, bg
  USE mp,                   ONLY: mp_barrier, mp_bcast, mp_sum, mp_end
  USE mp_global,            ONLY: mp_startup, nimage, npool, intra_image_comm, inter_image_comm, &
                                  nproc_pool, mpime, nproc, my_pool_id, me_pool, &
                                  mp_global_end, inter_pool_comm

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
COMPLEX(DP)               ::   aux(ngmsex), psic(nrxx), vpsi(ngm),auxsco(ngmsco)!, auxr(nrxxs)
COMPLEX(DP)               ::   ZDOTC, sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                               sigma_band_ex(nbnd_sig, nbnd_sig), vxc(nbnd_sig,nbnd_sig)
LOGICAL                   ::   do_band, do_iq, setup_pw, exst, single_line
INTEGER                   ::   iman, nman, ndeg(nbnd_sig), ideg, iq, ikq
COMPLEX(DP), ALLOCATABLE  ::   sigma(:,:,:)
COMPLEX(DP), ALLOCATABLE  ::   evc_tmp_j(:), evc_tmp_i(:)
INTEGER, ALLOCATABLE      ::   igkq_ig(:) 
INTEGER, ALLOCATABLE      ::   igkq_tmp(:) 
COMPLEX(DP)               ::   corr_element
REAL(DP)                  ::   dvoxel
!for analytic continuation of selfenergy:
REAL(DP)                  ::   wsigmax, wsigmin, deltaws
INTEGER                   ::   nwsigwin
REAL(DP), ALLOCATABLE     ::   wsigwin(:)

COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)

!For VXC matrix elements:
REAL(DP) :: vtxc, etxc, ehart, eth, charge

ALLOCATE (igkq_tmp(npwx))
ALLOCATE (igkq_ig(npwx))

     one   = 1.0d0 
     czero = (0.0d0, 0.0d0)
     w_ryd = wsigma/RYTOEV
     nbnd = nbnd_sig 

     iq = 1 
     xq(:) = xk_kpoints(:, ik0)
     lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
     setup_pw = .TRUE.
     do_band = .TRUE.

     if (lgamma) then
           ikq = ik0
        else
           ikq = 2*ik0
     endif

     write(stdout,'(/4x,"k0(",i3," ) = (",3f7.3," )")') ikq, (xk (ipol,ikq) , ipol = 1, 3)
     !WRITE(6,'("Running PWSCF")')
     !CALL run_pwscf(do_band)
     !CALL initialize_gw()

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
     CALL davcio (evc, lrwfc, iuwfc, 1, -1)
   else
!else then psi_{\k+\gamma = \psi_{k}} should be second entry in list.
     CALL davcio (evc, lrwfc, iuwfc, 2, -1)
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
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)

     do ir = 1, nrxx
        psic (ir) = psic(ir) * v%of_r (ir,1)
     enddo

     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)

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

!igkq_tmp is index of G-vector on gamma_centered grid.
!igkq_ig is the index of the g_vector evc(ig) on the k+G grid.
!counter is the total number which fall within some cutoff npwq,
!ngmsex, etc. 

 counter  = 0
 igkq_tmp(:) = 0
 igkq_ig(:)  = 0

 do ig = 1, npwq
      if((igkq(ig).le.ngmsex).and.((igkq(ig)).gt.0)) then
        counter = counter + 1
        igkq_tmp (counter) = igkq(ig)
        igkq_ig  (counter) = ig
    endif
 enddo

!WRITE(6,'("Number of G vectors for Sigma_ex", i4)') counter
!@10TION: only looping up to counter so need to watch that...
!setting these arrays to dim ngmsex lets us calculate all matrix elements with sigma as 
!a vector^{T}*matrix*vector product.

 ALLOCATE (sigma_g_ex (ngmsex, ngmsex))
 ALLOCATE (evc_tmp_i  (ngmsex))
 ALLOCATE (evc_tmp_j  (ngmsex))

 sigma_g_ex(:,:) = (0.0d0, 0.0d0)
 CALL davcio(sigma_g_ex, lrsex, iunsex, 1, -1)
 sigma_band_ex (:, :) = czero

 do ibnd = 1, nbnd_sig
     evc_tmp_i(:) = czero
  do jbnd = 1, nbnd_sig
     evc_tmp_j(:) = czero
     do ig = 1, counter
        evc_tmp_i(igkq_tmp(ig)) = evc(igkq_ig(ig), ibnd) 
     enddo
     do ig = 1, ngmsex
        do igp = 1, counter
           evc_tmp_j(igkq_tmp(igp)) = evc(igkq_ig(igp), jbnd)
        enddo
        do igp = 1, ngmsex
            aux(igp) = sigma_g_ex (igp, ig)
        enddo
           sigma_band_ex (ibnd, jbnd) = sigma_band_ex (ibnd, jbnd) + &
           evc_tmp_i (ig) * ZDOTC(ngmsex, evc_tmp_j (1:ngmsex), 1, aux, 1)
      enddo
  enddo
 enddo

 DEALLOCATE(sigma_g_ex)
 DEALLOCATE(evc_tmp_i)
 DEALLOCATE(evc_tmp_j)

 WRITE(6,*) 
 write(stdout,'(4x,"Sigma_ex (eV)")')
 write(stdout,'(8(1x,f7.3))') real(sigma_band_ex(:,:))*RYTOEV
 write(stdout,*)
 write(stdout,'(8(1x,f7.3))') aimag(sigma_band_ex(:,:))*RYTOEV

!MATRIX ELEMENTS OF SIGMA_C:
 WRITE(6,*) 
 WRITE(6,'("Sigma_C Matrix Element")') 
 ALLOCATE (sigma(ngmsco,ngmsco,nwsigma)) 
 ALLOCATE (evc_tmp_i(ngmsco))
 ALLOCATE (evc_tmp_j(ngmsco))

 counter     = 0
 igkq_tmp(:) = 0
 igkq_ig(:)  = 0

 do ig = 1, npwq
    if((igkq(ig).le.ngmsco).and.((igkq(ig)).gt.0)) then
        counter = counter + 1
        igkq_tmp (counter) = igkq(ig)
        igkq_ig  (counter) = ig
    endif
 enddo

 sigma = dcmplx(0.0d0, 0.0d0)
 do iw = 1, nwsigma
     CALL davcio (sigma(:,:, iw), lrsigma, iunsigma, iw, -1)
 enddo

!!!!!!!!REAL SPACE SIGMA_C
 sigma_band_c (:,:,:) = czero
!!!!!!!!
 WRITE(6,*) 
 WRITE(6,'("Number of G vectors for sigma_corr, npwq", 2i8)') counter, npwq
 WRITE(6,*) 
 sigma_band_c (:,:,:) = czero
  do ibnd = 1, nbnd_sig
     evc_tmp_i(:) = czero
   do jbnd = 1, nbnd_sig
      evc_tmp_j(:) = czero
      do iw = 1, nwsigma
       do ig = 1, counter
             evc_tmp_i(igkq_tmp(ig)) = evc(igkq_ig(ig), ibnd)
       enddo
      do ig = 1, ngmsco
            do igp = 1, counter
               evc_tmp_j(igkq_tmp(igp)) = evc(igkq_ig(igp), jbnd)
            enddo
            do igp = 1, ngmsco
           !With inversion symmetry these are perfectly symmetric ig -> igp
              auxsco(igp) = sigma (ig, igp, iw)
            enddo
           sigma_band_c (ibnd, jbnd, iw) = sigma_band_c (ibnd, jbnd, iw) + &
           (evc_tmp_i(ig))*ZDOTC(ngmsco, evc_tmp_j (1:ngmsco), 1, auxsco, 1)
      enddo
   enddo
  enddo
 enddo
DEALLOCATE (sigma) 

if (do_imag) then 
!Can set arbitrary sigma windows with analytic continuation:
    wsigmax   =  15.0
    wsigmin   = -10.0
    deltaws   =   0.1
    nwsigwin  = 1 + ceiling((wsigmax - wsigmin)/deltaws)
    allocate (wsigwin(nwsigwin))
    do iw = 1, nwsigwin
        wsigwin(iw) = wsigmin + (wsigmax-wsigmin)/float(nwsigwin-1)*float(iw-1)
    enddo
    allocate (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))

    write(6,*) wsigwin(:)
    sigma_band_con(:,:,:) = dcmplx(0.0d0, 0.d0)

!print imaginary representation
!   call print_matel(vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
!do analytic continuation and print real
    call sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
    call print_matel(ikq, vxc(1,1), sigma_band_ex(1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
else
    call print_matel(vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
endif

ENDIF

if(allocated(sigma_band_con)) deallocate(sigma_band_con)
CALL mp_barrier(inter_pool_comm)
RETURN
END SUBROUTINE sigma_matel


  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE sigma_exchG(ik0)
  USE kinds,                ONLY : DP
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,                 ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE control_gw,           ONLY : eta, nbnd_occ, trunc_2d, multishift, lgamma
  USE klist,                ONLY : wk, xk, nkstot, nks
  USE io_files,             ONLY : prefix, iunigk, wfc_dir
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE symm_base,        ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE cell_base,        ONLY : omega, tpiba2, at, bg, tpiba, alat
  USE eqv,              ONLY : evq, eprec
  USE units_gw,         ONLY : iunsex, lrsex, lrwfc, iuwfc
  USE qpoint,           ONLY : xq, npwq, igkq, nksq
  USE gwsigma,          ONLY : sigma_x_st, nbnd_sig, sigma_band_exg
  USE buffers,          ONLY : save_buffer, get_buffer
  USE io_global,        ONLY : stdout, ionode_id, ionode, meta_ionode
  USE gvect,            ONLY : nl, ngm, g, nlm, gstart, gl, igtongl
  USE mp,               ONLY : mp_sum, mp_barrier, mp_bcast
  USE noncollin_module, ONLY : npol, nspin_mag
  USE mp_pools,         ONLY : inter_pool_comm, npool, kunit, my_pool_id
  USE mp_world,         ONLY : nproc, mpime
  USE save_gw,          ONLY : tmp_dir_save
  USE mp_images,        ONLY : nimage, my_image_id, inter_image_comm
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE fft_base,         ONLY: dffts
  USE gvecs,            ONLY : nls
IMPLICIT NONE
!ARRAYS to describe exchange operator.
  COMPLEX(DP), ALLOCATABLE :: eigv(:,:)
  COMPLEX(DP), ALLOCATABLE :: barcoul(:)
  COMPLEX(DP), ALLOCATABLE :: miv(:), mvj(:)
  COMPLEX(DP), ALLOCATABLE :: psi(:), dpsic(:)
  COMPLEX(DP) :: czero, exch_element
  COMPLEX(DP) :: psik(npwx*npol, nbnd_sig)
  COMPLEX(DP) :: psikp(npwx*npol, nbnd_sig)
  COMPLEX(DP) :: dipole(dffts%nnr), matel
  COMPLEX(DP) :: pwg0(sigma_x_st%dfftt%nnr)
  COMPLEX(DP) :: phase
  COMPLEX(DP) :: ZdoTC
  REAL(DP)   :: dvoxel, wgt,nsymm1, sigma_ex_tr
  REAL(DP)   :: qg2, qg, qxy, qz
  REAL(DP)   :: rcut, spal, zcut
  REAL(DP)   :: xq_coul(3)
  REAL(DP)   :: xq_old(3)
  REAL(DP)   :: voxel
  REAL(DP)   :: sigma_ex_diag(nbnd_sig)
  REAL(DP)   :: xk1(3), aq(3)
  REAL(DP)   :: fac
  REAL(DP), parameter :: eps=1.e-5_dp
  INTEGER    :: ikmq, ik0, ik, igkdim
  INTEGER    :: ig, igp, npe, irr, icounter, ir, irp
  INTEGER    :: iq, ipol, ibnd, jbnd, vbnd, counter, counterk
  INTEGER    :: rec0, ios
  INTEGER    :: iman, nman, ndeg(nbnd_sig), ideg, ikq
  INTEGER    :: igkp(npwx) 
  INTEGER    :: igk_ig(npwx) 
  INTEGER    :: igk_tmp(npwx) 
  INTEGER    :: igkq_ig(npwx) 
  INTEGER    :: igkq_tmp(npwx) 
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  INTEGER    :: isym, isymop
  INTEGER    :: iqrec, nig0
  INTEGER    :: ik1old, npwkp
  INTEGER    :: ik1, ikstar, iq1
  LOGICAL    :: do_band, do_iq, setup_pw, exst, limit, single_line
  LOGICAL    :: found_q, inv_q
  LOGICAL    :: found_k

  call start_clock('sigma_exch')

  allocate ( gmapsym  (ngm, nrot)       )
  allocate ( eigv     (ngm, nrot)       )
  if       (.not.allocated(sigma_band_exg)) allocate(sigma_band_exg(nbnd_sig,num_k_pts))
  allocate ( psi  (dffts%nnr)  )
  allocate ( dpsic(dffts%nnr))
  allocate ( barcoul(npwx))

  write(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
  write(6,'(4x,"Occupied bands at Gamma: ",i3)') nbnd_occ(ik0)
  write(6,'(4x,"nksq,nks,kunit ",3i4)') nksq, nks, kunit
  write(6,'(4x,"nsym ",i4)') nsym
  write(stdout,'(4x,"Running Sigma_exchgq")')
  call gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
  czero = (0.0d0, 0.0d0)
  limit =.false.
  wgt = 1.0/omega
  nsymm1  = 1.0/(float(nsym))
  sigma_band_exg(:,ik0) = (0.0d0, 0.0d0)
  ik1old = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <nk | \Sum_{q} nk-q nk-q |nk>                !
! Need to collect k-point on all processors    !
! Then they loop over all local k-q points     !
! and sum matrix element over pools at the end !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  found_k    = .false.
  psikp(:,:) = dcmplx(0.0d0,0.0d0)
  aq(:)      = 0.0d0

! pick up psikp we are interested in !
  if(lgamma) then
    IF(my_pool_id==0) CALL get_buffer(psikp, lrwfc, iuwfc, 1)
  else
    IF(my_pool_id==0) CALL get_buffer(psikp, lrwfc, iuwfc, 2)
  endif
  CALL mp_barrier(inter_pool_comm)
  CALL mp_bcast(psikp, 0, inter_pool_comm)
  CALL gk_sort(xk_kpoints(1,ik0), ngm, g, ( ecutwfc / tpiba2 ), &
               npwkp, igkp, g2kin)
  !write(1000+mpime,*)sum(psikp(:,1)*conjg(psikp(:,1)))

  do ik = 1, nksq
     if(lgamma) ik1 = ik
     if(.not.lgamma) ik1 = 2*ik-1
     psik(:,:)  = dcmplx(0.0d0,0.0d0)
     call get_buffer (psik, lrwfc, iuwfc, ik1)
     do isymop = 1, nsym
        nig0  = 1
        counter  = 0
        igkq_tmp = 0
        igkq_ig  = 0

        call rotate(xk(1,ik1), aq, s, nsym, invs(isymop))
        xq  = xk_kpoints(:, ik0) - aq(:)

        call gk_sort(xk(1,ik1), ngm, g, ( ecutwfc / tpiba2 ), &
                     npw, igk, g2kin)
        npwq = npw 
        igkq = igk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! Construct W(q;G;G'): !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
        if(.not.trunc_2d) THEN
           !do ig = 1, sigma_x_st%ngmt
           do ig = 1, npw
! FFT means we have -G components of vq(G):
              qg = sqrt((g(1,ig)  + xq(1))**2.d0  + (g(2,ig) + xq(2))**2.d0  &
                      + (g(3,ig)  + xq(3))**2.d0)
              qg2 = (g(1,ig)  + xq(1))**2.d0  + (g(2,ig) + xq(2))**2.d0  &
                 + ((g(3,ig)) + xq(3))**2.d0
              limit = (qg.lt.eps8)
              if(.not.limit) then
                 spal = 1.0d0 - cos (rcut * tpiba * qg)
                 barcoul (ig) = &
&                e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
              else
                 barcoul(ig) = (fpi*e2*(rcut**2))/2
              endif
           enddo
        else
            zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat*nq3
            rcut = -2*pi*zcut**2
           !do ig = 1, sigma_x_st%ngmt
            do ig = 1, npw
               qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
               qxy  = sqrt((g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2)
               qz   = sqrt((g(3,ig) + xq(3))**2)
               spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
               if(qxy.gt.eps8) then
                  spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
                  barcoul(ig) = &
&                 dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
               else if(qxy.lt.eps8.and.qz.gt.eps8) then
                  spal = 1.0d0 - cos(tpiba*qz*zcut) - tpiba*qz*zcut*sin(tpiba*qz*zcut)
                  barcoul(ig) = &
&                 dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
               else  
                  barcoul(ig) = dcmplx(rcut, 0.0d0)
               endif
            enddo
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do vbnd = 1, nbnd_occ(ik0)
           psi (:) = (0.d0, 0.d0)
           do ig = 1, npw 
              psi(nls(gmapsym(igkq(ig), invs(isymop)))) = psik(ig, vbnd)
           enddo
           call invfft ('Wave', psi(:), dffts)
!Just diagonal elements of exchange operator
           do ibnd = 1, nbnd_sig
              dpsic(:) = (0.d0, 0.d0)
              do ig = 1, npwkp
                 dpsic(nls(igkp(ig))) = psikp(ig, ibnd)
              enddo
              call invfft ('Wave', dpsic(:), dffts)
              dipole(:) = (0.0d0, 0.0d0)
              do ir = 1, dffts%nnr
                 dipole (ir) = dipole (ir) + dconjg(psi (ir))*dpsic (ir)
              enddo
              call fwfft('Wave', dipole(:), dffts)
              !do ig = 1, sigma_x_st%ngmt
              do ig = 1, npw
                 sigma_band_exg(ibnd,ik0) = sigma_band_exg(ibnd,ik0) - &
&                0.5*wk(ik1)*dconjg(dipole(nls(ig)))*(dipole(nls(ig)))*barcoul(ig)
              enddo 
           enddo !ibnd
        enddo !v\inocc
      enddo !isym
 enddo ! on q

 call mp_barrier(inter_pool_comm)
 call mp_sum (sigma_band_exg(:,ik0), inter_pool_comm)

 sigma_band_exg(:, ik0) = wgt*nsymm1*sigma_band_exg(:, ik0)

 write(stdout,'(4x,"Sigma_ex (eV)")')
 write(stdout,*) real(sigma_band_exg(:, ik0))*RYTOEV
 write(stdout,*) aimag(sigma_band_exg(:, ik0))*RYTOEV
 !write(stdout,'(4x,"Sigma_ex val (eV)", 8(1x,f7.2))')  real(sigma_band_exg(1:8, ik0))*RYTOEV

 deallocate ( barcoul  )
 deallocate ( psi      )
 deallocate ( dpsic    )
 deallocate ( gmapsym  )
 deallocate ( eigv     )

 call stop_clock('sigma_exch')
 9000 format(8(1x,f7.2))
 9005 format(8(1x,f14.7))
END SUBROUTINE sigma_exchG

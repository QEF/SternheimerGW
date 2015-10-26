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
  USE control_gw,           ONLY : eta, nbnd_occ, trunc_2d, multishift
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
  USE mp,               ONLY : mp_sum, mp_barrier
  USE noncollin_module, ONLY : npol, nspin_mag
  USE mp_pools,         ONLY : inter_pool_comm, npool, kunit, my_pool_id
  USE mp_world,         ONLY : nproc, mpime
  USE save_gw,          ONLY : tmp_dir_save
  USE mp_images,        ONLY : nimage, my_image_id, inter_image_comm
  USE fft_interfaces,   ONLY : invfft, fwfft
IMPLICIT NONE
!ARRAYS to describe exchange operator.
  COMPLEX(DP), ALLOCATABLE :: eigv(:,:)
  COMPLEX(DP), ALLOCATABLE :: barcoul(:)
  COMPLEX(DP), ALLOCATABLE :: miv(:), mvj(:)
  COMPLEX(DP), ALLOCATABLE :: psi(:), dpsic(:)
  COMPLEX(DP) :: dipole(sigma_x_st%dfftt%nnr), matel
  COMPLEX(DP) :: czero, exch_element
  COMPLEX(DP) :: psik(npwx*npol, nbnd_sig)
  COMPLEX(DP) :: pwg0(sigma_x_st%dfftt%nnr)
  COMPLEX(DP) :: phase
  COMPLEX(DP) :: ZdoTC
  REAL(DP)   :: dvoxel, wgt,nsymm1, sigma_ex_tr
  REAL(DP)   :: qg2, qg, qxy, qz
  REAL(DP)    :: rcut, spal, zcut
  REAL(DP)   :: xq_coul(3)
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
  INTEGER    :: igk_ig(npwx) 
  INTEGER    :: igk_tmp(npwx) 
  INTEGER    :: igkq_ig(npwx) 
  INTEGER    :: igkq_tmp(npwx) 
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  INTEGER    :: isym, isymop
  INTEGER    :: iqrec, nig0
  INTEGER    :: ik1old
  INTEGER    :: ik1, ikstar, iq1
  LOGICAL    :: do_band, do_iq, setup_pw, exst, limit, single_line
  LOGICAL    :: found_q, inv_q
  LOGICAL    :: found_k

  call start_clock('sigma_exch')

  allocate ( gmapsym  (ngm, nrot)       )
  allocate ( eigv     (ngm, nrot)       )
  allocate ( sigma_band_exg(nbnd_sig)   )
  allocate ( psi(sigma_x_st%dfftt%nnr)  )
  allocate ( dpsic(sigma_x_st%dfftt%nnr))
  allocate ( miv(sigma_x_st%dfftt%nnr), mvj(sigma_x_st%dfftt%nnr))
  allocate ( barcoul(sigma_x_st%ngmt))

  write(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
  write(6,'(4x,"Occupied bands at Gamma: ",i3)') nbnd_occ(ik0)
  write(6,'(4x,"nksq,nks,kunit ",3i4)') nksq, nks, kunit
  write(6,'(4x,"nsym ",i4)') nsym
  call gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
  czero = (0.0d0, 0.0d0)
  limit =.false.
  wgt = 1.0/omega
  nsymm1  = 1.0/(float(nsym))
  sigma_band_exg(:) = (0.0d0, 0.0d0)

  ik1old = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <nk | \Sum_{q} nk-q nk-q |nk>                !
! Need to collect k-point on all processors    !
! Then they loop over all local k-q points     !
! and sum matrix element over pools at the end !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  found_k =.false.
  psik(:,:) = dcmplx(0.0d0,0.0d0)
  do iq = 1, nks
     found_k  = (abs(xk_kpoints(1,ik0) - xk(1,iq)).le.eps).and. &
                (abs(xk_kpoints(2,ik0) - xk(2,iq)).le.eps).and. &
                (abs(xk_kpoints(3,ik0) - xk(3,iq)).le.eps)
     if (found_k) then
        ikq = iq
        exit
     endif
  enddo
  if(found_k) then
     call get_buffer (psik, lrwfc, iuwfc, ikq)
  endif
  call mp_barrier(inter_pool_comm)
  call mp_sum(psik, inter_pool_comm)
  call mp_barrier(inter_pool_comm)
  !call gk_sort(xk_kpoints(1,ik0), ngm, g, ( ecutwfc / tpiba2 ), &
  !             npw, igk, g2kin)
  call gk_sort(xk_kpoints(1,ik0), ngm, g, ( ecutwfc / tpiba2 ), &
               npw, igk, g2kin)
  npwq = npw
  call wavecut(npw, counterk, igk, igk_tmp, igk_ig)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iq = 1, nqs
     do isymop = 1, nsym
        xq(:) = x_q(:,iq)
        call rotate(xq, aq, s, nsym, invs(isymop))
        xk1  = xk_kpoints(:, ik0) - aq(:)
        nig0  = 1
        call find_qG_ibz(xk1, s, iqrec, isym, nig0, found_q, inv_q)
! Checking if kpoint is in this pool:          !
        found_k = .false.
        do ikstar = 1, nks
           found_k  = (abs(xk(1, ikstar) - x_q(1, iqrec)).le.eps).and. &
                      (abs(xk(2, ikstar) - x_q(2, iqrec)).le.eps).and. & 
                      (abs(xk(3, ikstar) - x_q(3, iqrec)).le.eps) 
           if (found_k) then
              ik1 = ikstar 
              exit
           endif
        enddo

        if (found_k) then
            if (ik1.ne.ik1old) then
              call get_buffer (evc, lrwfc, iuwfc, ik1)
            endif
            ik1old = ik1
        else 
            cycle 
        endif
        
        call gk_sort(xk(1,ik1), ngm, g, ( ecutwfc / tpiba2 ), &
                     npw, igk, g2kin)
        npwq = npw
        call wavecut(npwq, counter, igk, igkq_tmp, igkq_ig)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! Construct W(q;G;G'): !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (inv_q) write(6,'("inv_q", i4)') nig0
        barcoul(:) = (0.0d0,0.0d0)
        rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
        if(.not.trunc_2d) THEN
           do ig = 1, sigma_x_st%ngmt
! FFT means we have -G components of vq(G):
              qg = sqrt((-g(1,ig)  + xq(1))**2.d0  + (-g(2,ig) + xq(2))**2.d0  &
                      + (-g(3,ig)  + xq(3))**2.d0)
              qg2 = (-g(1,ig)  + xq(1))**2.d0  + (-g(2,ig) + xq(2))**2.d0  &
                 + ((-g(3,ig)) + xq(3))**2.d0
              limit = (qg.lt.eps8)
              if(.not.limit) then
                 spal = 1.0d0 - cos (rcut * tpiba * qg)
                 barcoul (gmapsym(ig, invs(isymop))) = &
&                e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
              else
                 barcoul(gmapsym(ig,invs(isymop))) = (fpi*e2*(rcut**2))/2
              endif
           enddo
        else
            zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat*nq3
            rcut = -2*pi*zcut**2
            do ig = 1, sigma_x_st%ngmt
               qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
               qxy  = sqrt((g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2)
               qz   = sqrt((g(3,ig) + xq(3))**2)
               spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
               if(qxy.gt.eps8) then
                  spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
                  barcoul(gmapsym(ig,invs(isymop))) = &
&                 dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
               else if(qxy.lt.eps8.and.qz.gt.eps8) then
                  spal = 1.0d0 - cos(tpiba*qz*zcut) - tpiba*qz*zcut*sin(tpiba*qz*zcut)
                  barcoul(gmapsym(ig,invs(isymop))) = &
&                 dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
               else  
                  barcoul(gmapsym(ig,invs(isymop))) = dcmplx(rcut, 0.0d0)
               endif
            enddo
        endif
        if(nig0.gt.1) then
           pwg0(:) = dcmplx(0.0d0, 0.0d0)
           pwg0(sigma_x_st%nlt(nig0)) = dcmplx(1.0d0, 0.0d0)
           CALL invfft('Custom', pwg0(:), sigma_x_st%dfftt)
        else
           pwg0(:) = dcmplx(1.0d0, 0.0d0)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do vbnd = 1, nbnd_occ(ik0)
           psi (:) = (0.d0, 0.d0)
           do ig = 1, counter
              psi(sigma_x_st%nlt(gmapsym(igkq_tmp(ig),invs(isym)))) = evc(igkq_ig(ig), vbnd)
           enddo
           call invfft ('Custom', psi(:), sigma_x_st%dfftt)
!Just diagonal elements of exchange operator
           do ibnd = 1, nbnd_sig
              dpsic(:) = (0.d0, 0.d0)
              do ig = 1, counterk
                 dpsic (sigma_x_st%nlt(igk_tmp(ig))) = psik(igk_ig(ig), ibnd)
              enddo
              call invfft ('Custom', dpsic(:), sigma_x_st%dfftt)
              dipole(:) = (0.0d0, 0.0d0)
              if(.not. inv_q) then
                do ir = 1, sigma_x_st%dfftt%nnr
                   dipole (ir) = dipole (ir) + psi (ir)*pwg0(ir)*conjg(dpsic (ir))
                enddo
              else
!Use T.R. to recover \psi_{+k-q}:
                do ir = 1, sigma_x_st%dfftt%nnr
                   dipole (ir) = dipole (ir) + conjg((psi(ir))*(pwg0(ir)))*conjg(dpsic (ir))
                enddo
              endif
              call fwfft ('Custom', dipole(:), sigma_x_st%dfftt)
              miv(:) = (0.0d0,0.0d0)
              do ig = 1, sigma_x_st%ngmt
                 miv(ig) = dipole(sigma_x_st%nlt(ig))
              enddo
              do ig = 1, sigma_x_st%ngmt
                 sigma_band_exg(ibnd) = sigma_band_exg(ibnd) - &
&                wq(iq)*dconjg(miv(ig))*(miv(ig))*barcoul(ig)
              enddo 
           enddo !ibnd
        enddo !v\inocc
      enddo !isym
 enddo ! on q

 call mp_barrier(inter_pool_comm)
 call mp_sum (sigma_band_exg, inter_pool_comm)  
 sigma_band_exg = wgt*nsymm1*sigma_band_exg

 write(stdout,'(4x,"Sigma_ex (eV)")')
 write(stdout,*) real(sigma_band_exg(:))*RYTOEV
 write(stdout,*) aimag(sigma_band_exg(:))*RYTOEV
 write(stdout,'(4x,"Sigma_ex val (eV)", 8(1x,f7.2))')  real(sigma_band_exg(1:8))*RYTOEV

 call stop_clock('sigma_exch')
 9000 format(8(1x,f7.2))
 9005 format(8(1x,f14.7))
END SUBROUTINE sigma_exchG

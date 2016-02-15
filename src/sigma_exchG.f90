  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
subroutine sigma_exchG(ik0)
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,             ONLY : nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE control_gw,       ONLY : eta, nbnd_occ, trunc_2d, multishift, lgamma
  USE klist,            ONLY : wk, xk, nkstot, nks
  USE io_files,         ONLY : prefix, iunigk, wfc_dir
  USE wvfct,            ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE gvecw,            ONLY : ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE symm_base,        ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot,&
                               copy_sym, inverse_s, s_axis_to_cart
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
  complex(DP), allocatable :: eigv(:,:)
  complex(DP), allocatable :: barcoul(:)
  complex(DP), allocatable :: psi(:), dpsic(:)
  complex(DP) :: czero, exch_element
  complex(DP) :: psik(npwx*npol, nbnd_sig)
  complex(DP) :: psikp(npwx*npol, nbnd_sig)
  complex(DP) :: dipole(dffts%nnr), matel
  complex(DP) :: pwg0(sigma_x_st%dfftt%nnr)
  complex(DP) :: phase
  complex(DP) :: ZdoTC
  real(DP)   :: wgt, nsymm1, sigma_ex_tr
  real(DP)   :: qg2, qg, qxy, qz
  real(DP)   :: rcut, spal, zcut
  real(DP)   :: xq_coul(3)
  real(DP)   :: xq_old(3)
  real(DP)   :: sigma_ex_diag(nbnd_sig)
  real(DP)   :: xk1(3), aq(3)
  real(DP)   :: fac
  real(DP), parameter :: eps=1.e-5_dp
  integer    :: ikmq, ik0, ik, igkdim, nsymq
  integer    :: ig, igp, npe, irr, icounter, ir, irp
  integer    :: iq, ipol, ibnd, jbnd, vbnd
  integer    :: rec0, ios
  integer    :: iman, nman, ndeg(nbnd_sig), ideg, ikq
  integer    :: igkp(npwx) 
  integer    :: igk_ig(npwx) 
  integer    :: igk_tmp(npwx) 
  integer    :: igkq_ig(npwx) 
  integer    :: igkq_tmp(npwx) 
  integer, allocatable  :: gmapsym(:,:)
  integer    :: isym, isymop
  integer    :: iqrec, nig0
  integer    :: npwkp
  integer    :: ik1, ikstar, iq1
  logical    :: limit
  logical    :: found_q, inv_q, minus_q, sym(48)
!For Star of q.
  real(DP) :: sxq(3,48), xqs(3,48)
  integer  :: imq, isq(48), nqstar, nqs
  integer  :: nsq(48), i, nsymrot

  call start_clock('sigma_exch')

  allocate ( gmapsym  (ngm, nrot)       )
  allocate ( eigv     (ngm, nrot)       )
  allocate ( psi  (dffts%nnr))
  allocate ( dpsic(dffts%nnr))
  allocate ( barcoul(npwx))
  if       (.not.allocated(sigma_band_exg)) allocate(sigma_band_exg(nbnd_sig))

  write(stdout,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0,&
      &(xk_kpoints(ipol, ik0), ipol=1,3)
  !write(stdout,'(4x,"Occupied bands at Gamma: ",i3)') nbnd_occ(ik1)
  write(stdout,'(4x,"nksq,nks,kunit ",3i4)') nksq, nks, kunit
  write(stdout,'(4x,"nsym ",i4)') nsym
  write(stdout,'(4x,"Running Sigma_exchgq")')

  call gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

  czero = (0.0d0, 0.0d0)
  limit =.false.
  wgt = 1.0/omega
  nsymm1  = 1.0/(float(nsym))
  !nsymm1  = 1.0/(float(nq1*nq2*nq3))
  sigma_band_exg(:) = (0.0d0, 0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <nk | \Sum_{q} nk-q nk-q |nk>                !
! Need to collect k-point on all processors    !
! Then they loop over all local k-q points     !
! and sum matrix element over pools at the end !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  psikp(:,:) = dcmplx(0.0d0,0.0d0)
  barcoul(:) = dcmplx(0.0d0,0.0d0)
  npwkp   = 0 
  igk(:)  = 0
  igkq(:) = 0
  igkp(:) = 0
  aq(:)   = 0.0d0
! pick up psikp we are interested in !
  if(lgamma) then
    IF(my_pool_id==0) call get_buffer(psikp, lrwfc, iuwfc, 1)
  else
    IF(my_pool_id==0) call get_buffer(psikp, lrwfc, iuwfc, 2)
  endif
  call mp_barrier (inter_pool_comm)
  call mp_bcast   (psikp, 0, inter_pool_comm)
  call mp_barrier (inter_pool_comm)
  call gk_sort    (xk_kpoints(1,ik0), ngm, g, ( ecutwfc / tpiba2 ), &
                   npwkp, igkp, g2kin)

  call gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)
  do ik = 1, nksq
     if(lgamma) ik1 = ik
     if(.not.lgamma) ik1 = 2*ik-1
     psik(:,:)  = dcmplx(0.0d0,0.0d0)
     call get_buffer (psik, lrwfc, iuwfc, ik1)
     call gk_sort(xk(1,ik1), ngm, g, ( ecutwfc / tpiba2 ), &
                  npw, igk, g2kin)
     npwq = npw 
     igkq = igk
!only sum over small group of q
!sym(1:nsym) = .true.
!call smallg_q (xk(1,ik1), 1, at, bg, 1, s, ftau, sym, minus_q)
!nsymq       = copy_sym ( nsym, sym )
!call inverse_s()    
!call s_axis_to_cart()
!call gmap_sym(nsymq, s, ftau, gmapsym, eigv, invs)
     CALL star_q(xk(1,ik1), at, bg, nsym, s, invs, nqs, sxq, isq, nsq, imq, .false. )
!     write( 1000+mpime, * )
!     write( 1000+mpime, '(5x,a,i4)') 'Number of q in the star = ', nqs
!     write( 1000+mpime, '(5x,a)') 'List of q in the star:'
!     write( 1000+mpime, '(7x,i4,i4,i4,3f14.9)') (iq1, nsq(iq1), isq(iq1), (sxq(i,iq1), i=1,3), iq1=1,nqs)
!     write( 1000+mpime, '(7x,i4,i4)') (iq1, isq(iq1), iq1=1,nsym)
     do iq1 = 1, nqs
        nsymrot=0
        do isym=1,nsym
           if (isq(isym) == iq1) then
               nsymrot=nsymrot+1
               if (nsymrot == 1) isymop=isym
           endif
        enddo
        if(nsymrot == 0) then
          call errore('dfile_star','no symmetry relates q at star(q)',iq)
        endif
        nig0     = 1
        igkq_tmp = 0
        igkq_ig  = 0
!       call rotate(xk(1,ik1), aq, s, nsym, invs(isymop))
!       xq  = xk_kpoints(:, ik0) - aq(:)
        xq  = xk_kpoints(:, ik0) - sxq(:,iq1)
!into crystal
        !xk1 = aq
        !call cryst_to_cart(1, xk1, at, -1)
        !write(2000+mpime,'(4f10.5)') wk(ik1), aq(1:3)
        !write(1000+mpime,'("isymop, ", i4)') isymop
        !write(1000+mpime,'(4f10.5)') wk(ik1), aq(1:3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! Construct v(q;G;G'): !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
        barcoul(:) = dcmplx(0.0d0,0.0d0)
        if(.not.trunc_2d) THEN
           do ig = 1, npwx
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
           do ig = 1, npwx
              qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
              qxy  = sqrt((g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2)
              qz   = sqrt((g(3,ig) + xq(3))**2)
              spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
              if(qg2.gt.eps8) then
                 spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*(cos(tpiba*qz*zcut))
                 barcoul(ig) = dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
              else  
                 barcoul(ig) = dcmplx(rcut, 0.0d0)
              endif
           enddo
        endif
        do vbnd = 1, nbnd_occ(ik1)
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
              do ig = 1, npw
                 sigma_band_exg(ibnd) = sigma_band_exg(ibnd) - &
&                0.5d0*wk(ik1)*float(nsq(iq1))*dconjg(dipole(nls(ig)))*(dipole(nls(ig)))*barcoul(ig)
!&                0.5d0*wk(ik1)*dconjg(dipole(nls(ig)))*(dipole(nls(ig)))*barcoul(ig)
              enddo 
           enddo !ibnd
        enddo !v\inocc
      enddo !isym
 enddo ! on q

 call mp_barrier(inter_pool_comm)
 call mp_sum (sigma_band_exg, inter_pool_comm)
 call mp_barrier(inter_pool_comm)
 sigma_band_exg(:) = wgt*nsymm1*sigma_band_exg(:)
 !sigma_band_exg(:) = wgt*sigma_band_exg(:)

 write(stdout,'(4x,"Sigma_ex (eV)")')
 write(stdout,*) real(sigma_band_exg(:))*RYTOEV
 write(stdout,*) aimag(sigma_band_exg(:))*RYTOEV
!write(stdout,'(4x,"Sigma_ex val (eV)", 8(1x,f7.2))')  real(sigma_band_exg(1:8, ik0))*RYTOEV

 deallocate ( gmapsym  )
 deallocate ( eigv     )
 deallocate ( psi      )
 deallocate ( dpsic    )
 deallocate ( barcoul  )

 call stop_clock('sigma_exch')
 9000 format(8(1x,f7.2))
 9005 format(8(1x,f14.7))
end subroutine sigma_exchG

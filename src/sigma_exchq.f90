  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE sigma_exch(ik0)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE control_gw,    ONLY : eta, nbnd_occ, trunc_2d, multishift, lgamma
  USE klist,         ONLY : wk, xk, nkstot, nks
  USE io_files,      ONLY : prefix, iunigk, wfc_dir
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE cell_base,     ONLY : omega, tpiba2, at, bg, tpiba, alat
  USE eqv,           ONLY : evq, eprec
  USE units_gw,      ONLY : iunsex, lrsex, lrwfc, iuwfc
  USE qpoint,        ONLY : npwq, igkq, nksq
  USE gwsigma,       ONLY : sigma_x_st, nbnd_sig, gexcut
  USE buffers,       ONLY : save_buffer, get_buffer
  USE io_global,     ONLY : stdout, ionode_id, ionode, meta_ionode
  USE gvect,         ONLY : nl, ngm, g, nlm, gstart, gl, igtongl
  USE mp,            ONLY : mp_sum, mp_barrier
  USE noncollin_module, ONLY : npol, nspin_mag
  USE mp_pools,      ONLY : inter_pool_comm, npool, kunit, my_pool_id
  USE mp_world,      ONLY : nproc, mpime
  USE save_gw,       ONLY : tmp_dir_save
  USE mp_images,     ONLY : nimage, my_image_id, inter_image_comm
  USE mp_global,     ONLY : mp_global_end
IMPLICIT NONE
! ARRAYS to describe exchange operator.
! q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  COMPLEX(DP) :: sigma_band_ex(nbnd_sig, nbnd_sig)
  COMPLEX(DP), ALLOCATABLE :: sigma_ex(:,:)
  COMPLEX(DP), ALLOCATABLE :: greenf_na(:,:), greenf_nar(:,:)
  COMPLEX(DP), ALLOCATABLE :: barcoul(:,:), barcoulr(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_g_ex(:,:)
  COMPLEX(DP), ALLOCATABLE ::  eigv(:,:)
  COMPLEX(DP) :: ZdoTC
  COMPLEX(DP) :: czero, exch_element
  REAL(DP)    :: rcut, spal, zcut
  REAL(DP)    :: xq_coul(3)
  REAL(DP)    :: xq(3)
  REAL(DP)    :: qg2, qg, qxy, qz
  REAL(DP)    :: xk1(3), aq(3)
  REAL(DP)    :: sxq(3,48), xqs(3,48)
  REAL(DP), PARAMETER :: eps=1.e-5_dp
  INTEGER, ALLOCATABLE     :: gmapsym(:,:)
  INTEGER    :: igkq_ig(npwx), igkq_tmp(npwx) 
  INTEGER    :: ikq
  INTEGER    :: iq1,ik1,iqstart,iqstop
  INTEGER    :: iqrec, nig0, isym
  INTEGER    :: imq, isq(48), nqstar, nkpts
  INTEGER    :: nsq(48), i, ikstar
  INTEGER    :: isymop
  INTEGER    :: iunwfc1
  INTEGER    :: kpoolid(nkstot), iqrec1(nkstot)
  INTEGER    :: nbase, nksloc, rest, mypoolid
  INTEGER    :: ikmq, ik0, ik
  INTEGER    :: ig, igp, ir, irp, kcounter
  INTEGER    :: iq, ipol, ibnd, jbnd, counter, ios
  LOGICAL    :: limit
  LOGICAL    :: found_k
  LOGICAL    :: found_q, inv_q
  INTEGER*8  :: unf_recl
  CHARACTER (len=256) :: form_str 
  CHARACTER (len=256) :: tempfile
  CHARACTER (len=256) :: poolnum

#define DIRECT_IO_FACTOR 8
! Self-Energy grid:
! iGv
  call start_clock('sigma_exch')
  allocate ( sigma_ex    (sigma_x_st%dfftt%nnr, sigma_x_st%dfftt%nnr) )
  allocate ( gmapsym  (ngm, nrot)   )
  allocate ( eigv     (ngm, nrot)   )
  call gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)
  rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
  write(6,'("rcut ", f12.7)') rcut
  !We need to access all the wavefunctions, these should be all collected on the
  !ionode
  write(1000+mpime,*) trim(wfc_dir)
  write(1000+mpime,*) trim(tmp_dir_save)
  write(1000+mpime,*) trim(prefix) 
  write(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
  write(6,'(4x,"Occupied bands at k: ",i3)') nbnd_occ(ik0)
  write(6,*)"   lgamma, gexcut ", lgamma, gexcut
  write(6,'(4x,"nksq,nks,nkstot,kunit ",4i4)') nksq, nks,nkstot, kunit
  kcounter = 0
  czero = (0.0d0, 0.0d0)
  sigma_ex(:,:) = (0.0d0, 0.0d0)
!New pool parallel approach, we cycle if the kpoint isn't on the present
!node, all nodes cycles over full Brillouin zone for q.
!write(1000+mpime,*) nks, nksq, xk(:,1:nks)
  do ik = 1, nksq
     if(lgamma) ik1 = ik
     if(.not.lgamma) ik1 = 2*ik-1
     call get_buffer (evc, lrwfc, iuwfc, ik1)
     npwq = npw
     do isymop = 1, nsym
!Need a loop to find all plane waves below ecutsco when igkq 
!takes us outside of this sphere.  
        nig0  = 1
        counter  = 0
        igkq_tmp = 0
        igkq_ig  = 0

        call rotate(xk(1,ik1), aq, s, nsym, invs(isymop))
        xq(:)   = xk_kpoints(:,ik0) - aq(:)

        call gk_sort(xk(1,ik1), ngm, g, ( ecutwfc / tpiba2 ), &
                      npw, igk, g2kin )
        npwq = npw
        igkq = igk
!igkq = igk
        do ig = 1, npw
           if((igk(ig).le.sigma_x_st%ngmt).and.((igk(ig)).gt.0)) then
           !if((igk(ig).le.gexcut).and.((igk(ig)).gt.0)) then
               counter = counter + 1
               igkq_tmp (counter) = igk(ig)
               igkq_ig  (counter) = ig
           endif
        enddo
        WRITE(1000+mpime,*) counter
        allocate ( greenf_na   (sigma_x_st%ngmt, sigma_x_st%ngmt) )
!psi_{k+q}(r)psi^{*}_{k+q}(r')
        greenf_na = (0.0d0, 0.0d0)
        do ig = 1, counter
           do igp = 1, counter
              do ibnd = 1, nbnd_occ(1)
                 greenf_na(igkq_tmp(ig), igkq_tmp(igp)) = greenf_na(igkq_tmp(ig), igkq_tmp(igp)) +         &
                                                          tpi*(0.0d0, 1.0d0)*conjg(evc(igkq_ig(ig),ibnd))* &
                                                         (evc(igkq_ig(igp), ibnd))
              enddo
           enddo
        enddo
!Fourier transform of green's function
        allocate ( greenf_nar  (sigma_x_st%dfftt%nnr, sigma_x_st%dfftt%nnr)  )
        greenf_nar(:,:) = czero
        call fft6_g(greenf_na(1,1), greenf_nar(1,1), sigma_x_st, gmapsym, eigv(1,1), isymop, 1, 1)
        deallocate(greenf_na)
        allocate ( barcoul  (sigma_x_st%ngmt, sigma_x_st%ngmt) )
        rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
        barcoul(:,:) = (0.0d0,0.0d0)
        if(.not.trunc_2d) THEN
           do ig = 1, sigma_x_st%ngmt
           !do ig = 1, gexcut
              qg = sqrt((g(1,ig)  + xq(1))**2.d0  + (g(2,ig) + xq(2))**2.d0  &
                      + (g(3,ig)  + xq(3))**2.d0)
              qg2 = (g(1,ig)  + xq(1))**2.d0  + (g(2,ig) + xq(2))**2.d0  &
                 + ((g(3,ig)) + xq(3))**2.d0
              limit = (qg.lt.eps8)
              if(.not.limit) then
                 spal = 1.0d0 - cos (rcut * tpiba * qg)
                 barcoul (ig, ig) = e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
              else
                !barcoul(gmapsym(ig, invs(isymop)), gmapsym(ig, invs(isymop))) = (fpi*e2*(rcut**2))/2
                barcoul(ig, ig) = (fpi*e2*(rcut**2))/2
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
               barcoul(ig,ig) = &
&              dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
        else if(qxy.lt.eps8.and.qz.gt.eps8) then
               spal = 1.0d0 - cos(tpiba*qz*zcut) - tpiba*qz*zcut*sin(tpiba*qz*zcut)
               barcoul(ig,ig) = &
&              dcmplx(e2*fpi/(tpiba2*qg2)*spal, 0.0d0)
        else  
               barcoul(ig,ig) = dcmplx(rcut, 0.0d0)
        endif
       enddo
     endif
     allocate (barcoulr    (sigma_x_st%dfftt%nnr,  sigma_x_st%dfftt%nnr))
     barcoulr(:,:) = (0.0d0, 0.0d0)
     call fft6(barcoul(1,1), barcoulr(1,1), sigma_x_st, 1)
     deallocate(barcoul)
     sigma_ex = sigma_ex + 0.5*wk(ik1)*(1.0d0/dble(nsym))*&
&               (0.0d0,1.0d0)/tpi*greenf_nar*barcoulr
     deallocate(barcoulr)
     deallocate(greenf_nar)
   enddo!isym
  enddo!iq
  call mp_barrier(inter_pool_comm)!I think I need these after writes!
  call mp_sum (sigma_ex, inter_pool_comm)  
  call mp_sum (kcounter, inter_pool_comm)  
!  write(1000+mpime,*)"kcounter: ", kcounter
  allocate ( sigma_g_ex  (sigma_x_st%ngmt, sigma_x_st%ngmt) )
  sigma_g_ex(:,:) = (0.0d0,0.0d0)
  if (meta_ionode) THEN
      call fft6(sigma_g_ex, sigma_ex, sigma_x_st, -1)
      call davcio(sigma_g_ex, lrsex, iunsex, ik0,  1)
  endif
  deallocate(sigma_g_ex)
  deallocate (sigma_ex)
  call mp_barrier(inter_pool_comm)!I think I need these after writes!
  call mp_barrier(inter_image_comm)!I think I need these after writes!
  call stop_clock('sigma_exch')
end subroutine sigma_exch

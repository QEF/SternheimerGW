!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
SUBROUTINE sigma_exch(ik0)

  USE buffers,              ONLY : save_buffer, get_buffer
  USE cell_base,            ONLY : omega, tpiba2, at, bg, tpiba, alat
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE control_gw,           ONLY : eta, nbnd_occ, truncation, multishift, lgamma, output
  USE disp,                 ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE eqv,                  ONLY : evq
  USE fft6_module,          ONLY : fwfft6, invfft6
  USE gvect,                ONLY : nl, ngm, g, nlm, gstart, gl, igtongl
  USE gvecw,                ONLY : ecutwfc
  USE gwsigma,              ONLY : sigma_x_st, nbnd_sig, gexcut
  USE io_files,             ONLY : prefix, wfc_dir
  USE io_global,            ONLY : stdout, ionode_id, ionode, meta_ionode
  USE kinds,                ONLY : DP
  USE kinds_gw,             ONLY : i8b
  USE klist,                ONLY : wk, xk, nkstot, nks
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE mp_global,            ONLY : mp_global_end
  USE mp_images,            ONLY : nimage, my_image_id, inter_image_comm
  USE mp_pools,             ONLY : inter_pool_comm, npool, kunit, my_pool_id
  USE mp_world,             ONLY : nproc, mpime
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE qpoint,               ONLY : npwq, igkq, nksq
  USE save_gw,              ONLY : tmp_dir_save
  USE sigma_io_module,      ONLY : sigma_io_write_x
  USE symm_base,            ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE timing_module,        ONLY : time_sigma_x
  USE truncation_module,    ONLY : truncate
  USE units_gw,             ONLY : lrsex, lrwfc, iuwfc, iunsex
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin

  IMPLICIT NONE

! ARRAYS to describe exchange operator.
! q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  COMPLEX(DP) :: sigma_band_ex(nbnd_sig, nbnd_sig)
  COMPLEX(DP), ALLOCATABLE :: sigma_ex(:,:)
  COMPLEX(DP), ALLOCATABLE :: greenf_na(:,:), greenf_nar(:,:)
  COMPLEX(DP), ALLOCATABLE :: barcoul(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_g_ex(:,:)
  COMPLEX(DP), ALLOCATABLE ::  eigv(:,:)
  COMPLEX(DP) :: ZdoTC
  COMPLEX(DP) :: czero, exch_element
  REAL(DP)    :: rcut, spal, zcut
  REAL(DP)    :: xq_coul(3)
  REAL(DP)    :: xq(3), q_G(3)
  REAL(DP)    :: qg2, qg, qxy, qz
  REAL(DP)    :: xk1(3), aq(3)
  REAL(DP)    :: sxq(3,48), xqs(3,48)
  REAL(DP), PARAMETER :: eps=1.e-5_dp
  INTEGER, ALLOCATABLE     :: gmapsym(:,:)
  INTEGER    :: igk(npwx), igkq_ig(npwx), igkq_tmp(npwx) 
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
  INTEGER    :: ig, igp, ir, irp
  INTEGER    :: iq, ipol, ibnd, jbnd, counter, ios
  LOGICAL    :: limit
  LOGICAL    :: found_k
  LOGICAL    :: found_q, inv_q
  INTEGER(i8b) :: unf_recl
  CHARACTER (len=256) :: form_str 
  CHARACTER (len=256) :: tempfile
  CHARACTER (len=256) :: poolnum

#define DIRECT_IO_FACTOR 8
! Self-Energy grid:
! iGv
  CALL start_clock(time_sigma_x)

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

        !
        ! evaluate the bare Coulomb interaction in a truncated geometry
        !
        ALLOCATE(barcoul(sigma_x_st%dfftt%nnr, sigma_x_st%dfftt%nnr))
        barcoul = czero
        !
        DO ig = 1, sigma_x_st%ngmt
          !
          q_G = (xq + g(:,ig)) * tpiba
          barcoul(ig, ig) = truncate(truncation, q_G)
          !
        END DO ! ig

        !
        ! now Fourier transform to real space
        !
        CALL invfft6('Custom', barcoul, sigma_x_st%dfftt, sigma_x_st%nlt, omega)

        !
        ! Sigma = G V
        !
        sigma_ex = sigma_ex + 0.5 * wk(ik1) * (1.0d0 / REAL(nsym, KIND = dp)) * &
&                  CMPLX(0.0d0, 1.0d0, KIND = dp) / tpi * greenf_nar * barcoul

        DEALLOCATE(barcoul)
        DEALLOCATE(greenf_nar)

     END DO ! isym
  END DO ! iq

  call mp_sum (sigma_ex, inter_pool_comm)  

  !
  ! evaluate Fourier transform and write to file
  !
  IF (meta_ionode) THEN
    !
    ! Fourier transform Sigma_x(r, r') -> Sigma_x(G, G')
    !
    CALL fwfft6('Custom', sigma_ex, sigma_x_st%dfftt, sigma_x_st%nlt, omega)
    !
    ! create copy of the array to write to file
    ALLOCATE(sigma_g_ex(sigma_x_st%ngmt, sigma_x_st%ngmt))
    sigma_g_ex = sigma_ex(:sigma_x_st%ngmt, :sigma_x_st%ngmt)
    !
    ! write unformatted
    CALL davcio(sigma_g_ex, lrsex, iunsex, ik0,  1)
    ! write formatted
    CALL sigma_io_write_x(output%unit_sigma, ik0, sigma_g_ex)
    !
    DEALLOCATE(sigma_g_ex)
    !
  END IF ! meta_ionode

  DEALLOCATE(sigma_ex)

  CALL stop_clock(time_sigma_x)

end subroutine sigma_exch

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
SUBROUTINE green_linsys_shift_re (green, mu, iq)

  USE buffers,              ONLY : save_buffer, get_buffer
  USE cell_base,            ONLY : tpiba2
  USE check_stop,           ONLY : check_stop_now
  USE constants,            ONLY : degspin, pi, tpi, RYTOEV, eps8
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, current_iq, ext_recover, &
                                   eta, tr2_green, maxter_green, prec_shift
  USE disp,                 ONLY : nqs, x_q
  USE ener,                 ONLY : ef
  USE eqv_gw,               ONLY : evq, eprectot
  USE freq_gw,              ONLY : fiu, nfs, nfsmax, wgreen, deltaw, &
                                   w0pmw, nwgreen, nwcoul, wcoul
  USE gvect,                ONLY : g, ngm
  USE gvecw,                ONLY : ecutwfc
  USE gwsigma,              ONLY : sigma_c_st, ecutsco, ecutprec
  USE io_files,             ONLY : prefix
  USE io_global,            ONLY : stdout, ionode
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, wk, nkstot, nks, igk_k, ngk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_bands,             ONLY : nproc_bgrp, ntask_groups
  USE mp_global,            ONLY : nproc_pool_file, &
                                   nproc_bgrp_file, nproc_image_file
  USE mp_images,            ONLY : nimage, my_image_id, intra_image_comm,   &
                                   me_image, nproc_image, inter_image_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_world,             ONLY : nproc, mpime
  USE nlcc_gw,              ONLY : nlcc_any
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE timing_module,        ONLY : time_green
  USE units_gw,             ONLY : iuwfc, lrwfc, iuwfcna, iungreen, lrgrn
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin, et

  IMPLICIT NONE 

  !should be freq blocks...
  COMPLEX(DP) :: gr_A_shift(npwx, nwgreen)
  COMPLEX(DP) :: gr_A(npwx, 1), rhs(npwx , 1)
  COMPLEX(DP) :: gr(npwx, 1), ci, cw 
  COMPLEX(DP) :: green(sigma_c_st%ngmt, sigma_c_st%ngmt, nwgreen)
  COMPLEX(DP), ALLOCATABLE :: etc(:,:)

  REAL(DP) :: dirac, x, delta, support
  REAL(DP) :: k0mq(3) 
  REAL(DP) :: w_ryd(nwgreen)
  REAL(DP) , allocatable :: h_diag (:,:)
  REAL(DP)               :: eprec_gamma
  REAL(DP) :: thresh, anorm, averlt, dr2, sqrtpi
  REAL(DP) :: tr_cgsolve = 1.0d-4
  REAL(DP) :: ehomo, elumo, mu
  INTEGER :: iw, igp, iw0
  INTEGER :: iq, ik0
  INTEGER :: rec0, n1, gveccount
  INTEGER, ALLOCATABLE      :: niters(:)
  INTEGER :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ig,         & ! counter on G vectors
             ndim,       & ! integer actual row dimension of dpsi
             is,         & ! counter on spin polarizations
             nt,         & ! counter on types
             na,         & ! counter on atoms
             nrec, nrec1,& ! the record number for dvpsi and dpsi
             ios,        & ! integer variable for I/O control
             mode          ! mode index
   INTEGER     :: igk(npwx)
   INTEGER     :: igkq_ig(npwx) 
   INTEGER     :: igkq_tmp(npwx) 
   INTEGER     :: counter
   INTEGER :: igstart, igstop, ngpool, ngr, igs, ngvecs

   REAL(DP) :: gam(3)
   REAL(DP) :: xk1(3)

   LOGICAL :: conv_root
   EXTERNAL cg_psi, ch_psi_all_green

   CALL start_clock(time_green)

   allocate  (h_diag (npwx, 1))
   allocate  (etc(nbnd_occ(1), nkstot))
   ci = (0.0d0, 1.0d0)
!Convert freq array generated in freqbins into rydbergs.
   do iw =1, nwgreen
      w_ryd(iw) = w0pmw(1,iw)/RYTOEV
   enddo
   where_rec='no_recover'
   ikq = iq

   CALL gk_sort_safe(xk(:,ikq), ngm, g, (ecutwfc / tpiba2 ), &
                npw, igk, g2kin)
   igkq = igk
   npwq = npw

   ! store the extra elements in the expanded arrays
   igk_k(:, nks + 1) = igkq
   ngk(nks + 1)      = npw
   xk(:, nks + 1)    = xk(:,ikq)

    counter = 0
    igkq_tmp(:) = 0
    igkq_ig(:)  = 0 
    do ig = 1, npwx
       if((igkq(ig).le.sigma_c_st%ngmt).and.((igkq(ig)).gt.0)) then
           counter = counter + 1
           igkq_tmp (counter) = igkq(ig)
           igkq_ig  (counter) = ig
       endif
    enddo
    igstart = 1
    igstop = counter
    ngvecs = igstop-igstart + 1
    if(.not.allocated(niters)) ALLOCATE(niters(ngvecs))
    niters = 0 
    call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
    call get_buffer (evc, lrwfc, iuwfc, ikq)
    DO ig = 1, npwq
       g2kin (ig) = ((xk (1,ikq) + g (1, igkq(ig) ) ) **2 + &
                     (xk (2,ikq) + g (2, igkq(ig) ) ) **2 + &
                     (xk (3,ikq) + g (3, igkq(ig) ) ) **2 ) * tpiba2
    ENDDO
    green  = (0.0d0, 0.0d0)
    h_diag = 0.d0
    do ig = 1, npwq
       if(g2kin(ig).le.ecutprec) then
          h_diag(ig,1) =  1.0d0
       else
         if(prec_shift) then
            h_diag(ig,1)= 1.d0/max(1.0d0, g2kin(ig)/(eprectot(nbnd_occ(1),ikq)))
         else
            h_diag(ig,1) =  1.0d0
         endif
       endif
    enddo
    gveccount = 1
    gr_A_shift = (0.0d0, 0.d0)
    niters(:) = 0
    do ig = igstart, igstop
          rhs(:,:)  = (0.0d0, 0.0d0)
          rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
          gr_A(:,:) = (0.0d0, 0.0d0)
          lter = 0
          etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
          cw = CMPLX( 0, 0.0d0, kind=DP) 
          conv_root = .true.
          anorm = 0.0d0
          call cbcg_solve_green(ch_psi_all_green, cg_psi, etc(1,ikq), rhs, gr_A, h_diag,   &
                                npwx, npwq, tr2_green, ikq, lter, conv_root, anorm, 1, npol, &
                                cw , niters(gveccount), .true.)
          call green_multishift_im(npwx, npwq, nwgreen, niters(gveccount), 1, w_ryd(1), mu, gr_A_shift)
          if (ig.eq.igstop) write(1000+mpime,*) niters(gveccount)
          if (niters(gveccount).ge.maxter_green) then
                WRITE(1000+mpime, '(5x,"Gvec: ", i4)') ig
                gr_A_shift(:,:) = dcmplx(0.0d0,0.0d0)
          endif
          do iw = 1, nwgreen
             do igp = 1, counter
                green (igkq_tmp(ig), igkq_tmp(igp),iw) = green (igkq_tmp(ig), igkq_tmp(igp),iw) + &
                                                         gr_A_shift(igkq_ig(igp),iw)
                do ibnd = 1, nbnd_occ(ikq)
                   x = et(ibnd, ikq) - w_ryd(iw)
                   dirac = eta / pi / (x**2.d0 + eta**2.d0)
                   green(igkq_tmp(ig), igkq_tmp(igp), iw) = green(igkq_tmp(ig), igkq_tmp(igp), iw) + &
                                                            tpi*ci*conjg(evc(igkq_ig(ig), ibnd))   * &
                                                            evc(igkq_ig(igp), ibnd) * dirac
                enddo 
             enddo !igp
          enddo !iw
          gveccount = gveccount + 1
    enddo!ig
if(allocated(niters)) DEALLOCATE(niters)
if(allocated(h_diag)) DEALLOCATE(h_diag)
if(allocated(etc))    DEALLOCATE(etc)

  CALL stop_clock(time_green)

END SUBROUTINE green_linsys_shift_re

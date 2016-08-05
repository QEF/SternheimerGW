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
SUBROUTINE green_linsys_shift_im (green, xk1, iw0, mu, nwgreen)

  USE cell_base,            ONLY : tpiba2
  USE check_stop,           ONLY : check_stop_now
  USE constants,            ONLY : degspin, pi, tpi, RYTOEV, eps8
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, current_iq, ext_recover, &
                                   eta, tr2_green, maxter_green, prec_shift,&
                                   multishift, high_io
  USE disp,                 ONLY : nqs, x_q
  USE ener,                 ONLY : ef
  USE eqv_gw,               ONLY : evq, eprectot
  USE freq_gw,              ONLY : fiu, nfs, nfsmax, wgreen, deltaw, w0pmw
  USE gvect,                ONLY : g, ngm
  USE gvecw,                ONLY : ecutwfc
  USE gwsigma,              ONLY : sigma_c_st, ecutsco, ecutprec, gcutcorr
  USE io_files,             ONLY : prefix
  USE io_global,            ONLY : stdout, ionode
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, wk, nkstot, igk_k, nks, ngk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_bands,             ONLY : nproc_bgrp, ntask_groups
  USE mp_global,            ONLY : nproc_pool_file, &
                                   nproc_bgrp_file, nproc_image_file
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_images,            ONLY : nimage, my_image_id, intra_image_comm,   &
                                   me_image, nproc_image, inter_image_comm
  USE mp_world,             ONLY : nproc, mpime
  USE nlcc_gw,              ONLY : nlcc_any
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE parallel_module,      ONLY : parallel_task 
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE timing_module,        ONLY : time_green
  USE units_gw,             ONLY : iuwfc, lrwfc, iuwfcna, iungreen, lrgrn
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin, current_k

  IMPLICIT NONE 

 !should be freq blocks...
 !complex(DP) :: gr_A_shift(npwx, nwgreen)
  integer :: nwgreen
  complex(DP), allocatable :: gr_A_shift(:,:)
  complex(DP) :: gr_A(npwx, 1), rhs(npwx, 1)
  complex(DP) :: gr(npwx, 1), ci, cw 
  complex(DP) :: green(gcutcorr, gcutcorr, nwgreen)
  complex(DP), allocatable :: etc(:,:)
  real(DP) :: dirac, x, delta, support
  real(DP) :: k0mq(3) 
  real(DP) :: w_ryd(nwgreen)
  real(DP) , allocatable :: h_diag (:,:)
  real(DP)               :: eprec_gamma
  real(DP) :: thresh, anorm, averlt, dr2, sqrtpi
  real(DP) :: ehomo, elumo, mu
  real(DP) :: gam(3)
  real(DP), intent(in) :: xk1(3)
  integer :: iw, igp, iw0
  integer :: iq, ik0
  integer :: rec0, n1, gveccount
  integer, allocatable      :: niters(:)
  integer, allocatable :: num_task(:)
  integer :: kter,       & ! counter on iterations
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
  integer, target  :: igk(npwx)
  integer  :: igkq_ig(npwx) 
  integer  :: igkq_tmp(npwx) 
  integer  :: counter
  integer  :: igstart, igstop, ngpool, ngr, igs, ngvecs

  logical :: conv_root
  external   cg_psi, ch_psi_all_green

  call start_clock(time_green)

  ALLOCATE  (h_diag (npwx, 1))
  ALLOCATE  (etc(nbnd_occ(1), nkstot))
  if(multishift) ALLOCATE(gr_A_shift(npwx, nwgreen))
  ci = (0.0d0, 1.0d0)
! write(1000+mpime, *) mu
! Convert freq array generated in freqbins into rydbergs.
  do iw =1, nwgreen
     w_ryd(iw) = w0pmw(1,iw)/RYTOEV
  enddo
  where_rec='no_recover'
! hl arb k.
  call gk_sort_safe(xk1(1), ngm, g, (ecutwfc / tpiba2 ), &
                npw, igk, g2kin)
  igkq => igk
  npwq = npw

  ! store the extra elements in the expanded arrays
  igk_k(:, nks + 1) = igkq
  ngk(nks + 1)      = npw
  xk(:, nks + 1)    = xk1

! Need a loop to find all plane waves below ecutsco when igkq takes us outside of this sphere.
! igkq_tmp is gamma centered index up to ngmsco,
! igkq_ig  is the linear index for looping up to npwq.
  counter = 0
  igkq_tmp(:) = 0
  igkq_ig(:)  = 0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! TODO: Should only do Unique G-vectors same as in W       !!!
!!!!!!!       G_{Sk}(G,G')= G_{k}(S^-1G,S^-1G') where Sk = k     !!!
!!!!!!!       This was not implemented on move from 4.2.1 to 5.0 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ig = 1, npw
     if((igk(ig).le.gcutcorr).and.((igk(ig)).gt.0)) then
         counter = counter + 1
         igkq_tmp (counter) = igk(ig)
         igkq_ig  (counter) = ig
     endif
  enddo
!Now the G-vecs up to the correlation cutoff have been divided between pools.
!Calculates beta functions (Kleinman-Bylander projectors), with
!structure factor, for all atoms, in reciprocal space
!call init_us_2 (npwq, igkq, x_q (1, ikq), vkb)
!HL arb k 
    call init_us_2 (npw, igk, xk1(1), vkb)
    do ig = 1, npw
       g2kin (ig) = ((xk1 (1) + g (1, igk(ig) ) ) **2 + &
                     (xk1 (2) + g (2, igk(ig) ) ) **2 + &
                     (xk1 (3) + g (3, igk(ig) ) ) **2 ) * tpiba2
    enddo
    green  = dcmplx(0.0d0, 0.0d0)
    h_diag = 0.d0
    if(multishift) then
      do ig = 1, npw
           h_diag(ig,1) =  1.0d0
      enddo
    else
      do ig = 1, npw
!h_diag(ig,1)= 1.d0/max(1.0d0, g2kin(ig)/(eprectot(nbnd_occ(1),ikq)))
!HL arbk
         h_diag(ig,1)= 1.d0/max(1.0d0, g2kin(ig)/eprectot(nbnd_occ(1), 1))
      enddo
    endif
!On first frequency block we do the seed system with BiCG:
    gveccount = 1
    if(multishift) gr_A_shift = (0.0d0, 0.d0)
    call parallel_task(inter_image_comm, counter, igstart, igstop, num_task)
!allocate list to keep track of the number of residuals for each G-vector:
    ngvecs = igstop-igstart + 1
    if(.not.allocated(niters)) ALLOCATE(niters(ngvecs))
    niters(:) = 0
    do ig = igstart, igstop
!Doing Linear System with Wavefunction cutoff (full density) for each perturbation. 
          if(multishift) then
             rhs(:,:)            = (0.0d0, 0.0d0)
             rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
             gr_A(:,:)           = (0.0d0, 0.0d0)
             lter                = 0
             etc(:, :)           = CMPLX( 0.0d0, 0.0d0, kind=DP)
             cw                  = CMPLX( 0.0d0, 0.0d0, kind=DP) 
             conv_root           = .true.
             anorm               = 0.0d0
             call cbcg_solve_green(ch_psi_all_green, cg_psi, etc(1,1), rhs, gr_A, h_diag,   &
                                   npwx, npw, tr2_green, 1, lter, conv_root, anorm, 1, npol, &
                                   cw , niters(gveccount), .true.)
             call green_multishift_im(npwx, 2*gcutcorr, nwgreen, niters(gveccount), 1, w_ryd(1), mu, gr_A_shift)
             if (niters(gveccount).ge.maxter_green) then
                   WRITE(1000+mpime, '(5x,"Gvec: ", i4)') ig
                   gr_A_shift(:,:) = dcmplx(0.0d0,0.0d0)
             endif
             do iw = 1, nwgreen
               do  igp = 1, counter
                   green (igkq_tmp(ig), igkq_tmp(igp),iw) = green (igkq_tmp(ig), igkq_tmp(igp),iw) + &
                                                            gr_A_shift(igkq_ig(igp),iw)
               enddo
             enddo
             gveccount = gveccount + 1
          else if(.not.multishift) then
             do iw = 1, nwgreen
                rhs(:,:)  = (0.0d0, 0.0d0)
                rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
                gr_A(:,:) = (0.0d0, 0.0d0)
                lter = 0
                etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
                cw = CMPLX( mu, w_ryd(iw), kind=DP) 
                conv_root = .true.
                anorm = 0.0d0
                call cbcg_solve(ch_psi_all_green, cg_psi, etc(1,1), rhs, gr_A, h_diag,   &
                                npwx, npwq, tr2_green, 1, lter, conv_root, anorm, 1, npol, &
                                cw, maxter_green, .true.)
                if (lter.ge.maxter_green) then
                   WRITE(1000+mpime, '(5x,"Gvec: ", i4)') ig
                   gr_A = dcmplx(0.0d0,0.0d0)
                endif
                do igp = 1, counter
                   green (igkq_tmp(ig), igkq_tmp(igp),iw) = green (igkq_tmp(ig), igkq_tmp(igp),iw) + &
                                                            gr_A(igkq_ig(igp),1)
                enddo
             enddo
             gveccount = gveccount + 1
          endif
    enddo !igstart
#ifdef __PARA
    call mp_barrier(inter_image_comm)
    call mp_sum(green, inter_image_comm)
#endif __PARA
if(allocated(niters))     deallocate(niters)
if(allocated(h_diag))     deallocate(h_diag)
if(allocated(etc))        deallocate(etc)
if(allocated(gr_A_shift)) deallocate(gr_A_shift)

  call stop_clock(time_green)

end subroutine green_linsys_shift_im

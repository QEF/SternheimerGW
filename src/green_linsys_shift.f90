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
SUBROUTINE green_linsys_shift (ik0)
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin, pi, tpi, RYTOEV, eps8
  USE cell_base,            ONLY : tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : xk, wk, nkstot
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, current_iq, ext_recover, &
                                   eta, tr2_green, maxter_green
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iuwfc, lrwfc, iuwfcna, iungreen, lrgrn
  USE eqv,                  ONLY : evq
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE disp,                 ONLY : nqs
  USE freq_gw,              ONLY : fpol, fiu, nfs, nfsmax, nwgreen, wgreen, deltaw
  USE gwsigma,              ONLY : sigma_c_st
  USE gvect,         ONLY : g
  USE mp,            ONLY : mp_sum, mp_barrier
  USE mp_images,     ONLY : nimage, my_image_id, intra_image_comm,   &
                            me_image, nproc_image, inter_image_comm
  USE mp_global,     ONLY : nproc_pool_file, &
                            nproc_bgrp_file, nproc_image_file
  USE mp_bands,      ONLY : nproc_bgrp, ntask_groups
  USE mp_world,      ONLY : nproc, mpime

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

  INTEGER :: iw, igp, iwi
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
!HL need a threshold here for the linear system solver. This could also go in the punch card
!with some default at a later date. 
!Arrays to handle case where nlsco does not contain all G vectors required for |k+G| < ecut
    INTEGER     :: igkq_ig(npwx) 
    INTEGER     :: igkq_tmp(npwx) 
    INTEGER     :: counter

!PARALLEL
    INTEGER :: igstart, igstop, ngpool, ngr, igs, ngvecs

!tmp number of blocks
    INTEGER :: nblocks, block
    LOGICAL :: conv_root, testnan
    EXTERNAL cg_psi, ch_psi_all_green, testnan

     allocate  (h_diag (npwx, 1))
     allocate  (etc(nbnd_occ(ik0), nkstot))

    ci = (0.0d0, 1.0d0)
    nblocks = 1
   !We support the numerical delta fxn in a x eV window...
    support = 2.0d0/RYTOEV

!Convert freq array generated in freqbins into rydbergs.
    w_ryd(:) = wgreen(:)/RYTOEV
    CALL start_clock('greenlinsys')
    where_rec='no_recover'
    if (nksq.gt.1) rewind (unit = iunigk)
!Loop over q in the IBZ_{k}
do iq = 1, nksq 
        if (lgamma) then
            ikq = iq
            else
            ikq = 2*iq
        endif

        if (nksq.gt.1) then
            read (iunigk, err = 100, iostat = ios) npw, igk
 100        call errore ('green_linsys', 'reading igk', abs (ios) )
        endif
    
        if(lgamma) npwq=npw 

      if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
 200       call errore ('green_linsys', 'reading igkq', abs (ios) )
      endif
     
!write(1000+mpime,*) igkq 
!Need a loop to find all plane waves below ecutsco when igkq takes us outside of this sphere.
!igkq_tmp is gamma centered index up to ngmsco,
!igkq_ig  is the linear index for looping up to npwq.
!need to loop over...
      counter = 0
      igkq_tmp(:) = 0
      igkq_ig(:)  = 0 
      do ig = 1, npwx
         if((igkq(ig).le.sigma_c_st%ngmt).and.((igkq(ig)).gt.0)) then
             counter = counter + 1
            !index in total G grid.
             igkq_tmp (counter) = igkq(ig)
            !index for loops 
             igkq_ig  (counter) = ig
         endif
      enddo
      CALL para_img(counter, igstart, igstop)
      WRITE(6, '(5x, "iq ",i4, " igstart ", i4, " igstop ", i4)') iq, igstart, igstop
!allocate list to keep track of the number of residuals for each G-vector:
      ngvecs = igstop-igstart + 1
      if(.not.allocated(niters)) ALLOCATE(niters(ngvecs))
      niters = 0 

! Now the G-vecs up to the correlation cutoff have been divided between pools.
! Calculates beta functions (Kleinman-Bylander projectors), with
! structure factor, for all atoms, in reciprocal space
       call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
       call davcio (evq, lrwfc, iuwfc, ikq, - 1)
       do ig = 1, npwq
          g2kin (ig) = ((xk (1,ikq) + g (1, igkq(ig) ) ) **2 + &
                        (xk (2,ikq) + g (2, igkq(ig) ) ) **2 + &
                        (xk (3,ikq) + g (3, igkq(ig) ) ) **2 ) * tpiba2
       enddo
       WRITE(6, '(4x,"k0+q = (",3f12.7," )",10(3x,f7.3))') xk(:,ikq), et(:,ikq)*RYTOEV
       WRITE(6, '(4x,"tr2_green for green_linsys",e10.3)') tr2_green
       green  = (0.0d0, 0.0d0)
       h_diag = 0.d0
!No preconditioning with multishift
     do ig = 1, npwx
        h_diag(ig,1) =  1.0d0
     enddo
!On first frequency block we do the seed system with BiCG:
     gveccount = 1 
     do ig = igstart, igstop
       do block = 1, nblocks
             rhs(:,:)  = (0.0d0, 0.0d0)
             rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
             gr_A(:,:) = (0.0d0, 0.0d0) 
             lter = 0
             etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
             cw = CMPLX( 0, 0, kind=DP)
             conv_root = .true.
             anorm = 0.0d0
!Doing Linear System with Wavefunction cutoff (full density) for each perturbation. 
             if (block.eq.1) then
                 call cbcg_solve_green(ch_psi_all_green, cg_psi, etc(1,ikq), rhs, gr_A, h_diag,  &
                                       npwx, npwq, tr2_green, ikq, lter, conv_root, anorm, 1, npol, &
                                       cw, niters(gveccount))
             endif
             call green_multishift(npwx, npwq, nwgreen, niters(gveccount), 1, gr_A_shift)
             if (.not.conv_root) WRITE(1000+mpime, '(5x,"Gvec", 10i4, i4)') niters, ig
             if (niters(gveccount).ge.maxter_green) then
                 if (.not.conv_root) WRITE(1000+mpime, '(5x,"Gvec", i4)') ig
                 gr_A_shift(:,:) = dcmplx(0.0d0,0.0d0)
             endif
             do iw = 1, nwgreen
                do igp = 1, counter
                   if( testnan(real(gr_A_shift(igkq_ig(igp),iw))).or.testnan(aimag(gr_A_shift(igkq_ig(igp),iw)))) then
                       gr_A_shift (igkq_ig(igp), iw) = (0.0d0, 0.0d0)
                       write(1000+mpime, '("green nan")')
                   endif
                enddo
             enddo

             do iw = 1, nwgreen
                do igp = 1, counter
                   green (igkq_tmp(ig), igkq_tmp(igp),iw) = green (igkq_tmp(ig), igkq_tmp(igp),iw) + &
                                                            gr_A_shift(igkq_ig(igp),iw)
                enddo
             enddo
         gveccount = gveccount + 1
!Green's Fxn Non-analytic Component:
!HLGREEN TEST
         do iw = 1, nwgreen
           do igp = 1, counter
!should be nbnd_occ:
            do ibnd = 1, nbnd_occ(ikq)
               x = et(ibnd, ikq) - w_ryd(iw)
               dirac = eta / pi / (x**2.d0 + eta**2.d0)
               green(igkq_tmp(ig), igkq_tmp(igp), iw) = green(igkq_tmp(ig), igkq_tmp(igp), iw) + &
                                                        tpi*ci*conjg(evq(igkq_ig(ig), ibnd))   * &
                                                        evq(igkq_ig(igp), ibnd) * dirac
            enddo 
           enddo!igp
         enddo!iw
       enddo !blocks
     enddo !ig
#ifdef __PARA
    CALL mp_barrier(inter_image_comm)
    CALL mp_sum(green, inter_image_comm)
    CALL mp_barrier(inter_image_comm)
    if(ionode) then
#endif
      do iw = 1, nwgreen
         rec0 = (iw-1) * 1 * nksq + (iq-1) + 1
         CALL davcio(green(:,:,iw), lrgrn, iungreen, rec0, +1, ios)
      enddo
#ifdef __PARA
    endif
    CALL mp_barrier(inter_image_comm)
#endif
ENDDO !iq

if(allocated(niters)) DEALLOCATE(niters)
if(allocated(h_diag)) DEALLOCATE(h_diag)
if(allocated(etc))    DEALLOCATE(etc)
CALL stop_clock('greenlinsys')
RETURN
END SUBROUTINE green_linsys_shift

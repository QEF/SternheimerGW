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
SUBROUTINE do_stern()
  USE io_global,  ONLY : stdout, ionode_id, meta_ionode
  USE kinds,      ONLY : DP
  USE disp,       ONLY : nqs, num_k_pts, xk_kpoints, w_of_q_start, x_q
  USE gwsigma,    ONLY : sigma_c_st, gcutcorr
  USE gwsymm,     ONLY : ngmunique, ig_unique, use_symm, sym_friend, sym_ig
  USE control_gw, ONLY : done_bands, reduce_io, recover, tmp_dir_gw,&
                          ext_restart, bands_computed, bands_computed, nbnd_occ, &
                          do_q0_only, solve_direct, tinvert, lrpa, do_epsil
  USE freq_gw,    ONLY : nfs
  USE units_gw,   ONLY : lrcoul, iuncoul
  USE klist,      ONLY : lgauss
  USE mp_global,  ONLY : inter_image_comm, intra_image_comm, &
                         my_image_id, nimage, root_image
  USE mp,         ONLY : mp_sum, mp_barrier
  USE mp_world,   ONLY : mpime
  USE noncollin_module, ONLY : noncolin, nspin_mag

IMPLICIT NONE

  INTEGER :: iq, ik, ig, igstart, igstop, ios, iq1, iq2
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g(:,:,:,:)
  LOGICAL :: do_band, do_iq, setup_pw, exst, do_matel, lgamma
  COMPLEX(DP), ALLOCATABLE :: eps_m(:)

  allocate ( scrcoul_g( gcutcorr, gcutcorr, nfs, nspin_mag))
  allocate ( ig_unique( gcutcorr) )
  allocate ( sym_ig(gcutcorr))
  allocate ( sym_friend(gcutcorr))

  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  do_matel = .TRUE.

  if(lgauss) write(stdout, '(//5x,"SYSTEM IS METALLIC")')
  if(.not.do_epsil) then
      iq1 = w_of_q_start
      iq2 = nqs
  else
  ! In case we want to trace a line through the brillouin zone
  ! or get the screening for a particular grid q points (i.e. coulomb matel).
      iq1 = w_of_q_start
      iq2 = num_k_pts
  endif
    
  do iq = iq1, iq2
!Perform head of dielectric matrix calculation.
     call start_clock ('epsilq')
    do_matel = .FALSE.
    scrcoul_g(:,:,:,:) = (0.0d0, 0.0d0)
     lgamma = ALL(x_q(:,iq) == 0)
     if (lgamma) THEN
        allocate(eps_m(nfs))
        eps_m(:) = dcmplx(0.0d0,0.0d0)
        if(my_image_id.eq.0) THEN
           call prepare_q0(do_band, do_iq, setup_pw, iq)
           call run_nscf(do_band, do_matel, iq)
           call initialize_gw()
           call coulomb_q0G0(iq, eps_m)
           scrcoul_g(1,1,:,1) = eps_m
           write(stdout,'(5x, "epsM(0) = ", f12.7)') eps_m(1)
           write(stdout,'(5x, "epsM(iwp) = ", f12.7)') eps_m(2)
           call clean_pw_gw(iq, .FALSE.)
        endif
        call mp_barrier(inter_image_comm)
    endif
    call prepare_q(do_band, do_iq, setup_pw, iq)
    call run_nscf(do_band, do_matel, iq)
    call initialize_gw()
    if(use_symm) THEN
      write(stdout,'("")')
      write(stdout,'(5x, "SYMMETRIZING COULOMB Perturbations")')
      write(stdout,'("")')
      call stern_symm()
    else
      ngmunique = gcutcorr
      do ig = 1, gcutcorr
         ig_unique(ig) = ig
      enddo
    endif
       if(nimage.gt.1) then
          call para_img(ngmunique, igstart, igstop)
       else
          igstart = 1
          igstop = ngmunique
       endif
       write(stdout, '(5x, "iq ",i4, " igstart ", i4, " igstop ", i4)')iq, igstart, igstop
       call coulomb(iq, igstart, igstop, scrcoul_g)
       if(nimage.gt.1) THEN
          call mp_sum(scrcoul_g, inter_image_comm)
       endif
!Only the meta_image should write to file
       if (meta_ionode) THEN
         call unfold_w(scrcoul_g,iq)
         if(solve_direct.and.tinvert) write(1000+mpime, '("UNFOLDING, INVERTING, WRITING W")')
         if(solve_direct.and.tinvert) call invert_epsilon(scrcoul_g, lgamma, eps_m)
         call davcio(scrcoul_g, lrcoul, iuncoul, iq, +1, ios)
       endif
       if(allocated(eps_m)) deallocate(eps_m)
       call mp_barrier(inter_image_comm)
       call clean_pw_gw(iq, .FALSE.)
       if(do_q0_only) GOTO 126
       call print_clock ('epsilq')
       call stop_clock ('epsilq')
  enddo
126 continue 
   write(stdout, '("Finished Calculating Screened Coulomb")')
   deallocate( scrcoul_g )
   deallocate( ig_unique )
   deallocate( sym_ig )
   deallocate( sym_friend )
end SUBROUTINE do_stern

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
program gw
!-----------------------------------------------------------------------
!... This is the main driver of the Sternheimer-GW code.
!-----------------------------------------------------------------------
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init
  USE control_gw,      ONLY : do_sigma_exx, do_sigma_exxG, do_sigma_matel, do_coulomb,&
                              do_green, multishift, do_sigma_c, do_q0_only,&
                              do_imag, lgamma, output
  USE freq_gw,          ONLY : nwsigma, nwsigwin
  USE gwsigma,          ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE io_files,         ONLY : diropn
  USE units_gw,         ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
  USE wvfct,            ONLY : nbnd
  USE disp,             ONLY : num_k_pts, w_of_k_start, w_of_k_stop
  USE input_parameters, ONLY : max_seconds, force_symmorphic
  USE pp_output_mod,    ONLY : pp_output_open_all

  IMPLICIT NONE

  integer             :: iq, ik, ierr
  character (LEN=9)   :: codepw = 'PW'
  character (LEN=9)   :: code   = 'SGW'
  character (LEN=256) :: auxdyn
  logical             :: do_band, exst, do_matel

! Initialize MPI, clocks, print initial messages
  call mp_startup ( start_images=.true. )
  call environment_start ( code )
  call sgw_opening_message () 
! Initialize GW calculation, Read Ground state information.
  call gwq_readin()
  call check_stop_init()
  call check_initial_status(auxdyn)
! Initialize frequency grids, FFT grids for correlation
! and exchange operators, open relevant GW-files.
  call freqbins()
  call sigma_grids()
  call opengwfil()
! Calculation W
  if(do_coulomb) call do_stern()
  ik = 1
  do_band  = .TRUE.
  do_matel = .TRUE.
  if (do_q0_only) goto 127
! Calculation of CORRELATION energy \Sigma^{c}_{k}=\sum_{q}G_{k-q}{W_{q}-v_{q}}:
  if (do_imag) then
      do ik = w_of_k_start, w_of_k_stop
         call run_nscf(do_band, do_matel, ik)
         call initialize_gw()
         if (do_sigma_c.and.multishift) call diropn(iunresid, 'resid', lrresid, exst)
         if (do_sigma_c.and.multishift) call diropn(iunalphabeta, 'alphbet', lralphabeta, exst)
         if (do_sigma_c) call sigma_c_im(ik)
         if (do_sigma_c.and.multishift) then
            close(unit = iunresid, status = 'DELETE')
            close(unit = iunalphabeta, status = 'DELETE')
         endif
! Calculation of EXCHANGE energy \Sigma^{x}_{k}= \sum_{q}G_{k}{v_{k-S^{-1}q}}:
         if (do_sigma_exx .and. .not.do_sigma_exxG) then   
             call sigma_exch(ik)
         else if(do_sigma_exx .and. do_sigma_exxG) then
             call sigma_exchg(ik)
         endif
! Calculation of Matrix Elements <n\k| V^{xc}, \Sigma^{x}, \Sigma^{c}(iw) |n\k>:
         if (do_sigma_matel) then
           if (ik == w_of_k_start) then         
             call pp_output_open_all(num_k_pts, nbnd_sig, nwsigwin, nwsigma, output)
           end if
           call sigma_matel(ik)
         end if
         call clean_pw_gw(ik, .TRUE.)
      enddo
  else
      do ik = w_of_k_start, w_of_k_stop
         call run_nscf(do_band, do_matel, ik)
         call initialize_gw()
         if(do_sigma_c.and.multishift) call diropn(iunresid, 'resid', lrresid, exst)
         if(do_sigma_c.and.multishift) call diropn(iunalphabeta, 'alphbet',lralphabeta, exst)
         if(do_sigma_c) call sigma_c_re(ik)
         if (do_sigma_c.and.multishift) then
            close(unit = iunresid, status = 'DELETE')
            close(unit = iunalphabeta, status = 'DELETE')
         endif
         if (do_sigma_exx .and. .not.do_sigma_exxG) then   
             call sigma_exch(ik)
         else if(do_sigma_exx .and. do_sigma_exxG) then
             call sigma_exchg(ik)
         endif
         if (do_sigma_matel) call sigma_matel(ik)
         call clean_pw_gw(ik, .TRUE.)
      enddo
  endif
  127 continue
  call close_gwq(.TRUE.)
  call stop_gw( .TRUE. )
end program gw

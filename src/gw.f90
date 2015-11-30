!-----------------------------------------------------------------------
! Sternheimer-GW code: Electronic structure code for 
! performing Quasiparticle calculations.
! Copyright (C) 2015 Henry Lambert, Feliciano Giustino
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
program gw
!-----------------------------------------------------------------------
!... This is the main driver of the Sternheimer-GW code.
!-----------------------------------------------------------------------
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init
  USE control_gw,      ONLY : do_sigma_exx, do_sigma_exxG, do_sigma_matel, do_coulomb,&
                              do_green, multishift, do_sigma_c, do_q0_only,&
                              do_imag, lgamma
  USE gwsigma,          ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE io_files,         ONLY : diropn
  USE units_gw,         ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
  USE wvfct,            ONLY : nbnd
  USE disp,             ONLY : num_k_pts, w_of_k_start, w_of_k_stop
  USE input_parameters, ONLY : max_seconds, force_symmorphic

  IMPLICIT NONE

  INTEGER :: iq, ik, ierr
  LOGICAL :: do_band, do_iq, setup_pw, exst, do_matel
  CHARACTER (LEN=9)   :: code = 'SGW'
  CHARACTER (LEN=256) :: auxdyn

! Initialize MPI, clocks, print initial messages
  call mp_startup ( start_images=.true. )
  call environment_start ( code )

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

  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  do_matel = .TRUE.
  ik = 1

  if (do_q0_only) goto 127
! Calculation of Correlation energy \Sigma^{c}_{k}= \sum_{q}G_{k-q}{W_{q}-v_{q}}
  if (do_imag) then
      do ik = w_of_k_start, num_k_pts
         call run_nscf(do_band, do_matel, ik)
         call initialize_gw()
         if(do_sigma_c.and.multishift) call diropn(iunresid, 'resid', lrresid, exst)
         if(do_sigma_c.and.multishift) call diropn(iunalphabeta, 'alphbet',lralphabeta, exst)
         if(do_sigma_c) call sigma_c_im(ik)
         if(do_sigma_c.and.multishift) then
            close(unit = iunresid, status = 'DELETE')
            close(unit = iunalphabeta, status = 'DELETE')
         endif
         if(do_sigma_exx .and. .not.do_sigma_exxG) then   
            call sigma_exch(ik)
         else if(do_sigma_exx .and. do_sigma_exxG) then
            call sigma_exchg(ik)
         endif
!Calculate <n\k| V^{xc}, \Sigma^{x}, \Sigma^{c}(iw) |n\k>
         if(do_sigma_matel) call sigma_matel(ik)
         CALL clean_pw_gw(ik, .TRUE.)
      enddo
  else
      do ik = w_of_k_start, num_k_pts
         call run_nscf(do_band, do_matel, ik)
         call initialize_gw()
         if(do_sigma_c.and.multishift) call diropn(iunresid, 'resid', lrresid, exst)
         if(do_sigma_c.and.multishift) call diropn(iunalphabeta, 'alphbet',lralphabeta, exst)
         if(do_sigma_c) call sigma_c_re(ik)
      enddo
  endif
  127 continue
  call close_gwq(.TRUE.)
  call stop_gw( .TRUE. )
end program gw

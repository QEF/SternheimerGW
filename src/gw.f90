  !-----------------------------------------------------------------------
  ! Sternheimer-GW code: Electronic structure code for 
  ! performing Quasiparticle calculations.
  ! Copyright (C) 2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
PROGRAM gw
  !-----------------------------------------------------------------------
  !... This is the main driver of the Sternheimer-GW code. It
  !... closely mirrors the construction of the PHonon code.
  !-----------------------------------------------------------------------
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init
  USE control_gw,      ONLY : do_sigma_exx, do_sigma_matel, do_coulomb,&
                              do_green, multishift, do_sigma_c, do_q0_only,&
                              do_imag, lgamma
  USE gwsigma,         ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE io_files,        ONLY : diropn
  USE units_gw,        ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
  USE wvfct,           ONLY : nbnd
  USE disp,            ONLY : num_k_pts

  IMPLICIT NONE

  INTEGER :: iq, ik, ierr
  CHARACTER (LEN=9)   :: code = 'SGW'
  CHARACTER (LEN=256) :: auxdyn
  LOGICAL :: do_band, do_iq, setup_pw, exst, do_matel
  !
  ! Initialize MPI, clocks, print initial messages
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( code )
  ! Initialize GW calculation, Read Ground state information.
  CALL gwq_readin()
  CALL check_stop_init()
  CALL check_initial_status(auxdyn)
  ! Initialize frequency grids, FFT grids for correlation
  ! and exchange operators, open relevant GW-files.
  CALL freqbins()
  CALL sigma_grids()
  CALL opengwfil()
  ! Calculation W
  IF(do_coulomb) CALL do_stern()
  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  do_matel = .TRUE.
  ik = 1
  if(do_q0_only) GOTO 127
  CALL run_nscf(do_band, do_matel, ik)
  CALL initialize_gw()
  IF(do_green.and.multishift) CALL diropn(iunresid, 'resid', lrresid, exst)
  IF(do_green.and.multishift) CALL diropn(iunalphabeta, 'alphbet',lralphabeta, exst)
  ! Calculation \Sigma_{k}= \sum_{q}G_{k-q}W_{q}
  IF (do_imag) THEN
      DO ik = 1, num_k_pts
         IF(do_sigma_c) CALL sigma_c_im(ik)
      ENDDO
  ELSE
      DO ik = 1, num_k_pts
         IF(do_sigma_c) CALL sigma_c_re(ik)
      ENDDO
  ENDIF
  IF(do_green.and.multishift) then
     CLOSE(UNIT = iunresid, STATUS = 'DELETE')
     CLOSE(UNIT = iunalphabeta, STATUS = 'DELETE')
  ENDIF
  do ik = 1, num_k_pts
     IF(do_sigma_exx)   CALL sym_sigma_exch(ik)
  enddo
  do ik = 1, num_k_pts
     IF(do_sigma_matel) CALL sigma_matel(ik)
  enddo
  127 CONTINUE
  CALL close_gwq(.TRUE.)
  CALL stop_gw( .TRUE. )
END PROGRAM gw

PROGRAM gw
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the Sternheimer-GW code. It
  ! ... closely mirrors the construction of the PHonon code.
  !-----------------------------------------------------------------------
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init
  USE control_gw,      ONLY : do_sigma_exx, do_sigma_matel, do_coulomb,&
                              do_green, multishift, do_sigma_c
  USE gwsigma,         ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE io_files,        ONLY : diropn
  USE units_gw,        ONLY : iunresid, lrresid, iunalphabeta, lralphabeta
  USE wvfct,           ONLY : nbnd

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
!Initialize GW calculation, Read Ground state information.
  CALL gwq_readin()
  CALL check_stop_init()
  CALL check_initial_status(auxdyn)
!Initialize frequency grids, FFT grids for correlation
!and exchange operators, open relevant GW-files.
  CALL freqbins()
  CALL sigma_grids()
  CALL opengwfil()
  IF(do_coulomb) CALL do_stern()
  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  do_matel = .TRUE.
  ik=1
  CALL run_nscf(do_band, do_matel, ik)
  CALL initialize_gw()
  IF(do_green.and.multishift) CALL diropn(iunresid, 'resid', lrresid, exst)
  IF(do_green.and.multishift) CALL diropn(iunalphabeta, 'alphbet',lralphabeta, exst)
  IF(do_green)   CALL green_linsys_shift(ik)
  IF(do_sigma_c) CALL sigma_c(ik)
  IF(do_green.and.multishift) then
     CLOSE(UNIT = iunresid, STATUS = 'DELETE')
     CLOSE(UNIT = iunalphabeta, STATUS = 'DELETE')
  ENDIF
  IF(do_sigma_exx)   CALL sigma_exch(ik)
  IF(do_sigma_matel) CALL sigma_matel(ik)
  CALL stop_gw( .TRUE. )
END PROGRAM gw

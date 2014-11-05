PROGRAM gw
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the Sternheimer-GW code. It
  ! ... closely mirrors the construction of the PHonon code.
  !-----------------------------------------------------------------------
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init
  USE control_gw,      ONLY : do_sigma_exx, do_sigma_matel
  USE gwsigma,         ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE wvfct,           ONLY : nbnd

  IMPLICIT NONE
  INTEGER :: iq, ik, ierr
  CHARACTER (LEN=9)   :: code = 'SGW'
  CHARACTER (LEN=256) :: auxdyn
  LOGICAL :: do_band, do_iq, setup_pw, exst
  !
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

  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  nbnd = nbnd_sig
  ik=1
  CALL run_pwscf_green(do_band,ik)
  CALL initialize_gw()
  if(do_sigma_exx)   CALL sigma_exch(ik)
  if(do_sigma_matel) CALL sigma_matel(ik)
  CALL stop_gw( .TRUE. )
END PROGRAM gw

PROGRAM gw
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the Sternheimer-GW code. It
  ! ... closely mirrors the construction of the PHonon code.
  !-----------------------------------------------------------------------
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init
  USE control_gw,      ONLY : do_sigma_exx

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

  CALL gwq_readin()

  CALL check_stop_init()

  CALL check_initial_status(auxdyn)

  CALL freqbins()

  CALL sigma_grids()

  CALL opengwfil()

  ik=1
  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
! Generates small group of k and then forms IBZ_{k}.
  CALL run_pwscf_green(do_band)
  CALL initialize_gw()
  if(do_sigma_exx)   CALL sigma_exch(ik)
! if(do_sigma_matel) CALL sigma_matel(ik)

  CALL stop_gw( .TRUE. )
END PROGRAM gw

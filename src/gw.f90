PROGRAM gw
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the Sternheimer-GW code. It
  ! ... closely mirrors the construction of the PHonon code.
  !-----------------------------------------------------------------------

  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE check_stop,      ONLY : check_stop_init

  IMPLICIT NONE

  INTEGER :: iq, ierr
  LOGICAL :: do_band, do_iq, setup_pw
  CHARACTER (LEN=9)   :: code = 'SGW'
  CHARACTER (LEN=256) :: auxdyn
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

!  if(do_sigma_exx)   CALL sigma_exch(ik)
!  if(do_sigma_matel) CALL sigma_matel(ik)

  CALL stop_gw( .TRUE. )

END PROGRAM gw

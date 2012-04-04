!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE run_pwscf_green(do_band)
!-----------------------------------------------------------------------
!
! This is the driver for when gw calls pwscf.
!
!
  USE control_flags,   ONLY : conv_ions, twfcollect
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE input_parameters,ONLY : pseudo_dir
  USE control_flags,   ONLY : restart
  USE qpoint,          ONLY : xq
  USE control_gw,      ONLY : done_bands, reduce_io, recover, tmp_dir_gw, &
                              ext_restart, bands_computed
  USE save_gw,         ONLY : tmp_dir_save
  USE control_flags,   ONLY: iprint
  USE gvect,           ONlY: ecutwfc
  USE gwsigma,         ONLY: ecutsco, ecutsex
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: dirname, file_base_in, file_base_out
  !
  LOGICAL, INTENT(IN) :: do_band
  !
  LOGICAL :: exst
  !
  CALL start_clock( 'PWSCF' )
  !
  CALL clean_pw( .FALSE. )
  !
  CALL close_files()
  !

  !From now on, work only on the _gw virtual directory
  !Somehow this statement got deleted.

  tmp_dir=tmp_dir_gw

  ! write(6,*)tmp_dir_gw
  ! ... Setting the values for the nscf run

  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  restart = ext_restart
  pseudo_dir= TRIM( tmp_dir_save ) // TRIM( prefix ) // '.save'

! write(6,*) pseudo_dir

  CALL restart_from_file()
  conv_ions=.true.

 !Generate all eigenvectors in IBZ_{k}.
  CALL setup_nscf_green (xq)

  CALL init_run()

  IF (do_band) write(6,'("Calling PW electrons")')
  IF (do_band) CALL electrons()

  IF (.NOT.reduce_io.and.do_band) THEN
     twfcollect=.FALSE. 
     CALL punch( 'all' )
     done_bands=.TRUE.
  ENDIF

  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )


  CLOSE( UNIT = 4, STATUS = 'DELETE' )

  ext_restart=.FALSE.

  CALL close_files()

  bands_computed=.TRUE.


  CALL stop_clock( 'PWSCF' )

  RETURN
END SUBROUTINE run_pwscf_green

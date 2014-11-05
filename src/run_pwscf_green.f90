!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE run_pwscf_green(do_band, ik)
!-----------------------------------------------------------------------
!
! This is the driver for when gw calls pwscf.
!
!
  USE kinds,              ONLY : DP
  USE control_flags,   ONLY : conv_ions, twfcollect
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir, seqopn
  USE lsda_mod,        ONLY : nspin
  USE input_parameters,ONLY : pseudo_dir
  USE control_flags,   ONLY : restart
  USE fft_base,        ONLY : dffts
  USE qpoint,          ONLY : xq
  USE control_gw,      ONLY : done_bands, reduce_io, recover, tmp_dir_gw, &
                              ext_restart, bands_computed, lgamma
  USE save_gw,         ONLY : tmp_dir_save
  USE control_flags,   ONLY : iprint
  USE mp_bands,        ONLY : ntask_groups
  USE disp,            ONLY : xk_kpoints
  USE gwsigma,         ONLY : sigma_x_st, sigma_c_st, nbnd_sig
  USE wvfct,           ONLY : nbnd
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: dirname, file_base_in, file_base_out
  !
  INTEGER   :: ik
  !
  LOGICAL, INTENT(IN) :: do_band
  !
  LOGICAL :: exst
  !
  CALL start_clock( 'PWSCF' )
  !
  CALL clean_pw( .FALSE. )
  !
  CALL close_files( .true. )
  !
  !From now on, work only on the _gw virtual directory
  wfc_dir=tmp_dir_gw
  tmp_dir=tmp_dir_gw
  !
  xq(:) = xk_kpoints(:, ik)
  lgamma = ( (ABS(xq(1))<1.D-12).AND.(ABS(xq(2))<1.D-12).AND.(ABS(xq(3))<1.D-12) )
  !...Setting the values for the nscf run
  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  restart = ext_restart
  conv_ions=.true.
  !Generate all eigenvectors in IBZ_{k}.
  CALL setup_nscf_green (xq)
  CALL init_run()
  IF (do_band) CALL non_scf ( )
  !
  CALL punch( 'all' )
  !
  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  ext_restart=.FALSE.
  !
  CALL close_files(.true.)
  !
  bands_computed=.TRUE.
  !
  !  PWscf has run with task groups if available, but in the phonon 
  !  they are not used, apart in particular points, where they are
  !  activated.
  !
  IF (ntask_groups > 1) dffts%have_task_groups=.FALSE.
  !
  CALL stop_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE run_pwscf_green

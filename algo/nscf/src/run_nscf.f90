!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
MODULE run_nscf_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE run_nscf(do_band, do_matel, ik, config)
!-----------------------------------------------------------------------
!
! This is the driver for when gw calls pwscf.
!
!
  USE basis,             ONLY: starting_wfc, starting_pot, startingconfig
  USE cell_base,         ONLY: at, bg
  USE check_stop,        ONLY: check_stop_now
  USE control_flags,     ONLY: io_level, conv_ions, twfcollect, restart
  USE control_gw,        ONLY: reduce_io, tmp_dir_gw
  USE control_lr,        ONLY: lgamma
  USE disp,              ONLY: xk_kpoints, nqs
  USE fft_base,          ONLY: dffts, dfftp
  USE fft_types,         ONLY: fft_type_allocate
  USE gvect,             ONLY: gcutm
  USE gvecs,             ONLY: gcutms
  USE gwsigma,           ONLY: nbnd_sig
  USE input_parameters,  ONLY: pseudo_dir, force_symmorphic
  USE io_files,          ONLY: tmp_dir, wfc_dir, seqopn, iunwfc
  USE io_global,         ONLY: stdout
  USE klist,             ONLY: nks, nkstot
  USE mp_bands,          ONLY: intra_bgrp_comm
  USE qpoint,            ONLY: xq
  USE setup_nscf_module, ONLY: setup_nscf_green, sigma_config_type
  USE wvfct,             ONLY: nbnd
  !
  IMPLICIT NONE
  !
  !> must be present if do_matel is set to contain the configuration of sigma
  TYPE(sigma_config_type), INTENT(OUT), ALLOCATABLE, OPTIONAL :: config(:)
  !
  INTEGER   :: ik
  !
  LOGICAL, INTENT(IN) :: do_band, do_matel
  !
  LOGICAL :: exst, opend
  !
  CALL start_clock( 'PWSCF' )
  !
  CALL clean_pw( .FALSE. )
  !
  CALL close_files( .true. )
  !
  IF (do_matel) THEN
    xq(:) = xk_kpoints(:, ik)
    IF (.NOT.PRESENT(config)) &
      CALL errore(__FILE__, "config is not optional when do_matel is set", 1)
  END IF
  lgamma = ( (ABS(xq(1))<1.D-8).AND.(ABS(xq(2))<1.D-8).AND.(ABS(xq(3))<1.D-8) )
 !From now on, work only on the _gw virtual directory
  wfc_dir=tmp_dir_gw
  tmp_dir=tmp_dir_gw
 !
 !...Setting the values for the nscf run
  startingconfig = 'input'
  starting_pot   = 'file'
  starting_wfc   = 'atomic'
  restart        = .FALSE.
  conv_ions=.true.
! Generate all eigenvectors in IBZ_{k} for Green's function or IBZ_{q} otherwise.
  if(do_matel) nbnd = nbnd_sig
  !
  CALL fft_type_allocate(dfftp, at, bg, gcutm, intra_bgrp_comm)
  CALL fft_type_allocate(dffts, at, bg, gcutms, intra_bgrp_comm)
  !
  IF (do_matel) THEN
    CALL setup_nscf_green(xq, config)
  ELSE
    CALL setup_nscf(xq)
  END IF
  CALL init_run()
  WRITE( stdout, '(/,5X,"Calculation of q = ",3F12.7)') xq
  IF (do_band) CALL non_scf ( )
  IF (.NOT.reduce_io.and.do_band) THEN
     IF (nks == 1 .and. io_level < 1) THEN
       ! punch opens the wavefunction file, so we need to close them if they are open
       INQUIRE(UNIT=iunwfc,OPENED=opend)
       IF (opend) CLOSE (UNIT=iunwfc, STATUS='keep')
     END IF
     twfcollect=.FALSE.
     CALL punch( 'all' )
  ENDIF


  IF(do_matel.and.nkstot.ne.nqs) THEN
    WRITE(stdout,'("WARNING: You have given a kpoint not in original BZ.'// &
                   'This could mean full symmetry is not exploited.")') 
  ENDIF

  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  CALL close_files(.true.)
  !
  CALL stop_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE run_nscf

END MODULE run_nscf_module

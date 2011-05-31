!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine gwq_recover
  !-----------------------------------------------------------------------
  !
  !    This subroutine tests if a xml restart file exists with the
  !    information of where the code stopped and, if appropriate the
  !    partial dynamical matrix and the partial effective charges.
  !    if (rec_code>2) done_irr, comp_irr
  !    info on calculated irreps - overrides initialization in gwq_setup.
  !    The xml file is in the
  !    directory _gwprefix.gwsave. The xml file contains
  ! where_rec  a string with information of the point where the calculation
  !            stopped
  !   rec_code_read  where_rec     status description
  !
  !  -1000                    Nothing has been read. There is no recover file.
  !  -40         gwq_setup    Only the displacements u have been read from file
  !  -30         gwq_init     u and dyn(0) read from file
  !  -25                      not active yet. Restart in solve_e_fpol
  !  -20         solve_e      all previous. Stopped within solve_e. There 
  !                           should be a recover file.
  !  -10         solve_e2     epsilon and zstareu are available if requested. 
  !                           Stopped within solve_e2. There should be a 
  !                           recover file.
  !   2          gwescf       all previous, raman tenson and elop tensor are
  !                           available if required.
  !   10         solve_linter all previous. Stopped within solve linter. 
  !                           There should be a recover file.
  !   20         gwqscf       all previous dyn_rec(irr) and zstarue0(irr) are
  !                           available.
  !   30         dynmatrix    all previous, dyn and zstarue are available.

  !
  ! The logic of the gwonon code recover is the following:
  ! The recover variable is read from input and never changed. If it is
  ! false it disables completely the recover.
  ! The control of the code is given by the arrays:
  ! comp_iq, done_iq : for each q point if it has to be calculated or
  !                    if it is already available. These are calculated 
  !                    only once by check_initial_status or read from file
  !                    by the same routine.
  ! comp_irr, done_irr : for each irreducible representation if it has
  !                      to be calculated or if it is already calculated.
  !                      The latter variables are valid only for the current
  !                      q and are calculated in gwq_setup and modified here
  !                      if something is on the file.
  ! epsil, done_epsil, zeu, done_zeu, zue, done_zue, lraman, done_lraman,
  ! elop, done_elop ... control the electric field calculations. These are
  ! set by prepare_q, or read from file by gwq_setup.
  !
  ! The position where the code stopped is in the variable rec_code_read
  ! defined above. This variable allows to exit from a routine if the quantity 
  ! calculated by this routine is already saved on file. 
  ! It is the responsibility of the routine (not of the calling code)
  ! to known if it has to make the calculation or just exit because the
  ! value of rec_code_read is too high.
  !
  ! if rec_code_read = (-25), -20, -10, 10  
  !    It is expected that an unformatted recover file exists. 
  !    The following data are in the 
  !    unformatted file and are read by
  !    routines solve_e (-20), solve_e2 (-10), solve_linter (10):
  ! iter, dr2, convt
  !    info on status of linear-response calculation for a given irrep.
  ! dvscfin
  !    self-consistent potential for current iteration and irrep
  ! if (okpaw) dbecsum
  !    the change of the D coefficients calculated so far.
  ! if (okvan) int1, int2, int3
  !    arrays used with US potentials : int1 and int2 calculated in dvanqq, 
  !    int3 calculatec in newdq (depends upon self-consistency)
  !
  !    rec_code_read is valid only for the first q. For the following q
  !    it is reset to -1000 in clean_pw_gw. So the recover file allows to
  !    restart only the current q. However information on other q could
  !    be available in the directory gwsave, so this routine reads the
  !    appropriate files and reset comp_irr and done_irr if appropriate.
  !
  !    NB: The restart of the electron-gwonon part is not available yet.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE gw_restart,    ONLY : gw_readfile
  USE control_gw,    ONLY : epsil, rec_code_read, all_done, where_rec,&
                            zeu, done_epsil, done_zeu, ext_recover, recover
  USE partial,       ONLY : comp_irr, done_irr
  USE modes,         ONLY : nirr

  !
  implicit none
  !
  integer :: irr, ierr
  ! counter on representations
  ! error code
  logical :: exst
  character(len=256) :: filename

  ierr=0
  IF (recover) THEN 
     CALL gw_readfile('data',ierr)
     IF (rec_code_read==-40) THEN
        WRITE( stdout, '(/,4x," Modes are read from file ")') 
     ELSEIF (rec_code_read==-20) THEN
        WRITE( stdout, '(/,4x," Restart in Electric Field calculation")')
     ELSEIF (rec_code_read==-10) then
        WRITE( stdout, '(/,4x," Restart in Raman calculation")') 
     ELSEIF (rec_code_read==2) THEN
        WRITE( stdout, '(/,4x," Restart after Electric Field calculation")')
     ELSEIF (rec_code_read==10.OR.rec_code_read==20) then
        WRITE( stdout, '(/,4x," Restart in Phonon calculation")')
     ELSEIF (rec_code_read==30) then
        WRITE( stdout, '(/,4x," Restart after Phonon calculation")')
     ELSE
        call errore ('gwq_recover', 'wrong restart data file', -1)
        ierr=1
     ENDIF
  ENDIF
!
  ext_recover = ext_recover .AND. ierr==0
!
! The case in which everything has been already calculated and we just
! recollect all the results must be treated in a special way (it does 
! not require any initialization). 
! We check here if everything has been done
!

  RETURN
END SUBROUTINE gwq_recover

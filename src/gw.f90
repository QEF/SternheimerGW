! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
PROGRAM gw
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the GW code.
  ! ... It reads all the quantities calculated by pwscf, it
  ! ... checks if some recover file is present and determines
  ! ... which calculation needs to be done. Finally, it makes
  ! ... a loop over the q points. At a generic q, if necessary, it
  ! ... recalculates the band structure by calling pwscf again.
  ! NC = norm conserving pseudopotentials
  ! US = ultrasoft pseudopotentials
  ! PAW = projector augmented-wave
  ! [1] LDA, [2] [1]+GGA, [3] [2]+LSDA/sGGA, [4] [3]+Spin-orbit/nonmagnetic,
  ! [5] [4]+Spin-orbit/magnetic
  !

  USE io_global,       ONLY : stdout
  USE control_gw,      ONLY : bands_computed, nbnd_occ
  USE wvfct,           ONLY : nbnd
 !HL epsil, trans, elph,
  USE disp,            ONLY : nqs
  USE output,          ONLY : fildrho
  USE check_stop,      ONLY : check_stop_init
  USE gw_restart,      ONLY : gw_writefile, destroy_status_run
  USE save_gw,         ONLY : clean_input_variables
  USE mp_global,       ONLY: mp_startup, nimage
  USE path_io_routines, ONLY : io_path_start
  USE environment,     ONLY: environment_start
  USE freq_gw,     ONLY : nfs, nwsigma
  USE units_gw,    ONLY : iuncoul, iungreen, lrgrn, lrcoul, iunsigma, lrsigma
  USE gwsigma,    ONLY : ngmsig
  
  IMPLICIT NONE
  INTEGER :: iq, ik
  INTEGER :: ios
  LOGICAL :: do_band, do_iq, setup_pw, exst
  CHARACTER (LEN=9)   :: code = 'GW'
  CHARACTER (LEN=256) :: auxdyn

  ! Initialize MPI, clocks, print initial messages
  ! /Modules/mp.f90

#ifdef __PARA
  CALL mp_startup ( )
  IF (nimage>1) CALL io_path_start()
#endif
   
  !/Modules/environment.f90 prints out all parallel information and opening message.
    CALL environment_start ( code )
    CALL gwq_readin()

    WRITE(stdout, '(/5x, "Finished reading variables")')

   !HL
   ! Check stop init Modules/check_stop.f90
   ! This module contains functions to check if the code should
   ! be smoothly stopped.

    CALL check_stop_init()

   ! This routine checks the initial status of the GW run, initializes the qmesh, and prepares
   ! the control of the dispersion calculation. 

    CALL check_initial_status(auxdyn)

   !Generate frequency grid for GW convolution and G-vector refold mapping.

    CALL freqbins()

   !CALL refold()
   !Generate cutoff for Sigma and gvector correspondence.
    CALL ggensig()

!   Coulomb file
    iuncoul = 28
    lrcoul = 2 * ngmsig * ngmsig * nfs
    CALL diropn (iuncoul, 'coul', lrcoul, exst)

!   Green's function file
    iungreen = 31
    lrgrn  = 2 * ngmsig * ngmsig

    CALL diropn (iungreen, 'green', lrgrn, exst)

!   Sigma file
    iunsigma = 32
    lrsigma = 2 * ngmsig * ngmsig * nwsigma
    CALL diropn(iunsigma, 'sigma', lrsigma, exst)

!CALCULATE W(r,r';iw)
GOTO 123
   DO iq = 1, nqs
        !comparing vkbs and g2kins
        !ik = 1
        !CALL prepare_kmq(do_band, do_iq, setup_pw, iq, ik)
        CALL prepare_q(do_band, do_iq, setup_pw, iq)
        CALL run_pwscf(do_band)
        !Initialize the quantities which do not depend on
        !the linear response of the system
        CALL initialize_gw()
       !CALCULATE W(r,r';iw)
        CALL coulomb(iq)
        CALL clean_pw_gw(iq)
   END DO

   WRITE(stdout, '("Finished Calculating Screened Coulomb")') 

!CALCULATE G(r,r'; w) 
    WRITE(stdout, '(/5x, "GREEN LINEAR SYSTEM SOLVER")')
    DO ik = 1, 1
      DO iq = 1, nqs

!For debug we calculate alot more states.
!       nbnd = 40
!       nbnd_occ = 4 

        CALL prepare_kmq(do_band, do_iq, setup_pw, iq, ik)
        CALL run_pwscf(do_band)
        CALL initialize_gw()

!       CALL green_linsys_test(ik, iq)
        CALL green_linsys(ik, iq)

!       WRITE(stdout, '(/5x, "Done Green_linsys")') 
        CALL clean_pw_gw(iq)
      ENDDO
    ENDDO

    WRITE(stdout, '("Finished Calculating Greens Function")') 
!123 CONTINUE

123 CONTINUE

    DO ik = 1, 1
       CALL gw_product(ik)
    ENDDO

    WRITE(6, '("Finished CALCULATING SIGMA")') 

    DO ik = 1, 1
       CALL sigma_matel(ik) 
    ENDDO
STOP

   CALL gw_writefile('init',0)
   CALL clean_input_variables()
   CALL collect_grid_files()
   CALL destroy_status_run()

   IF (bands_computed) CALL print_clock_pw()
   CALL stop_gw( .TRUE. )
END PROGRAM gw

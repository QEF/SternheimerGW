!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE prepare_kmq(do_band, do_iq, setup_pw, iq, ik)
  !-----------------------------------------------------------------------

  !  This routine prepares a few variables that are needed to control
  !  the GW run after the k-q point has been decided, but before
  !  doing the band calculation. 
    
  !  In particular if ldisp=true it sets:
  !  xq : the q point for the GW calculation
  !  current_iq : the current q point
  !  do_iq : if .true. q point has to be calculated
  !  setup_pw : if .true. the pw_setup has to be run
  !  do_band : if .true. the bands need to be calculated before phonon
  
  USE control_flags,   ONLY : modenum
  USE io_global,       ONLY : stdout, ionode
  USE klist,           ONLY : lgauss, xk
  USE qpoint,          ONLY : xq
  USE disp,            ONLY : x_q, done_iq, rep_iq, done_rep_iq, comp_iq
  USE control_gw,      ONLY : ldisp, lgamma, epsil, trans, zue, zeu, &
                              start_irr, last_irr, current_iq, &
                              done_bands
  USE freq_gw,         ONLY : fpol
  USE output,          ONLY : fildyn
  USE gw_restart,      ONLY : gw_writefile

  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq, ik
  LOGICAL, INTENT(OUT) :: do_band, do_iq, setup_pw
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER :: irr
  !
  do_iq=.TRUE.
  !
  ! Case 1) This q point is not calculated because not requested in this run
  !
  IF ( comp_iq(iq)==0 ) THEN
     do_iq=.FALSE.
     RETURN
  ENDIF
  !
  !  Case 2) This q point is not calculated because it has too few 
  !          representation and the starting representation is larger 
  !          than the number of available representations
  !
  current_iq = iq
  
  ! ... set the name for the output file
  ! ... set the "q" point in this case it is actually the "k-q" point but it is easier to
  ! stay in line with all the other routines as they are for now. Ideally could modify run_pwscf
  ! so it does everything in one shot for the green's function. (lots of useless other kpoints at the moment.)
  ! Hack fix for q = 0 again
  !@10TION

IF ( ldisp ) THEN

!HL obsolete name for output file.
!   fildyn = TRIM( auxdyn ) // TRIM( int_to_char( iq ) )

!   HL Don't need to reset gamma!  
     IF ( x_q(1,iq) == 0.D0 .AND. x_q(2,iq) == 0.D0 .AND. x_q(3,iq) == 0.D0 ) x_q(1,iq) = 0.00d0

     xq(:)  = xk (:, ik) - x_q(:,iq)

     lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )

     !epsil = .FALSE.
     !zue   = .FALSE.
     !zeu   = .FALSE.

ENDIF


  ! write(6,*) xq, x_q(:,iq), xk(:,ik)
  ! ... for q /= 0 no calculation of the dielectric tensor,
  ! ...            Born eff. charges, electro-optic, raman or
  ! ...            frequency dependent tensor
  !
  !    epsil = .FALSE.
  !    zue   = .FALSE.
  !    zeu   = .FALSE.
  !    fpol  =  .FALSE.
  !
  !  Save the current status of the run: all the flags, the list of q,
  !  and the current q, the fact that we are before the bands
  !
  CALL gw_writefile('init',0)
  !
  ! ... In the case of q != 0, we make first a non selfconsistent run
  !
  !HL this setup_pw and do bands are obsolete, they were for checking
  !that certain phonon modes had been computed etc.  
  !setup_pw = (.NOT.lgamma.OR.modenum /= 0).AND..NOT. done_bands

  do_band=.FALSE.

  !DO irr=start_irr, MIN(ABS(last_irr),rep_iq(iq))
  !   IF (done_rep_iq(irr,iq) /= 1) THEN
  !      do_band=.TRUE.
  !      EXIT
  !   ENDIF
  !ENDDO

!There are two special cases. When start_irr=0 and last_irr=0 we generate only
!the displacement patterns, and do not calculate the bands. If this q
!has been already calculated we only diagonalize the dynamical matrix
!IF ( start_irr == 0 .AND. last_irr == 0 ) do_band=.FALSE.
!IF ( done_iq(iq) == 1 ) do_band=.FALSE.

!@HL bit of a mess with old PH way of keeping track of which q points have been calculated.
!I'm just commenting all of this and saying if prepare kq is being called its probably because
!I want to calculate the eigenvectors at these k/q points.  

  setup_pw = .TRUE.
  do_band = .TRUE.

  WRITE( stdout, '(/,5X,"Calculation of k-q = ",3F12.7)') xq
  WRITE(6,*) setup_pw, do_band, lgamma

  RETURN
  !
END SUBROUTINE prepare_kmq

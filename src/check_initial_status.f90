!
! Copyright (C) 2012-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE check_initial_status(auxdyn)
  !-----------------------------------------------------------------------
  ! This routine checks the initial status of the phonon run and sets
  ! the variables that control the run, dealing with the image
  ! and GRID parallelization features of the phonon code.
  ! The size of the grid is determined by the following variables:
  ! nqs : the number of q points
  ! x_q : the coordinates of the q points
  ! nfs : the number of imaginary frequencies
  ! fiu : which frequencies 
  ! The flags that control which tensors to calculate
  !
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : modenum
  USE ions_base,       ONLY : nat
  USE io_files,        ONLY : tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE scf,             ONLY : rho
  USE disp,            ONLY : nqs, x_q, nq1, nq2, nq3
  USE qpoint,          ONLY : xq
  USE output,          ONLY : fildyn
  USE control_gw,      ONLY : ldisp, recover, done_bands,  &
                              start_q, last_q, current_iq, tmp_dir_gw, lgamma, &
                              ext_recover, ext_restart
  USE save_gw,         ONLY : tmp_dir_save
  USE units_gw,        ONLY : iudyn
  USE io_rho_xml,      ONLY : write_rho
  USE mp_images,       ONLY : nimage, intra_image_comm
  USE io_global,       ONLY : ionode, ionode_id
  USE io_files,        ONLY : prefix
  USE mp,              ONLY : mp_bcast
  USE xml_io_base,     ONLY : create_directory
  USE mp_global,       ONLY : mp_global_end
  !
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: auxdyn, filename
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER :: iq, iq_start, ierr
  !
  tmp_dir=tmp_dir_gw
  !
  ! If this not a recover run, we generate the q mesh. Otherwise at this
  ! point the code has read the q mesh from the files contained in 
  ! prefix.phsave
  !
  IF (.NOT.recover) THEN
     !
     ! recover file not found or not looked for
     !
     current_iq=1
     IF (ldisp) THEN
        !
        ! ... Calculate the q-points for the dispersion
        !
        CALL q_points()
        !
     ELSE
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xq(:)
        !
     END IF
  ENDIF
  !  Create a new directory where the gw variables are saved and copy
  !  the charge density there.
     CALL write_rho( rho, nspin )
  !  The following commands deal with the flag lqdir=.true. In this case
  !  each q point works on a different directory. We create the directories
  !  if they do not exist and copy the self consistent charge density
  !  there.
  !
  !DO iq = 1,nqs
  !   IF (.NOT.comp_iq(iq)) CYCLE
  !   lgamma = lgamma_iq(iq) 
  !   !
     ! ... each q /= gamma works on a different directory. We create them
     ! here and copy the charge density inside
     !
  !   IF ((.NOT.lgamma.OR. newgrid).AND.lqdir) THEN
  !      tmp_dir_phq= TRIM (tmp_dir_ph) //TRIM(prefix)//&
  !                        & '.q_' // TRIM(int_to_char(iq))//'/'
  !      filename=TRIM(tmp_dir_phq)//TRIM(prefix)//'.save/charge-density.dat'
  !      IF (ionode) inquire (file =TRIM(filename), exist = exst)
        !
  !      CALL mp_bcast( exst, ionode_id, intra_image_comm )
        !
  !      IF (.NOT. exst) THEN
  !         CALL create_directory( tmp_dir_phq )
  !         tmp_dir=tmp_dir_gwq
  !         CALL write_rho( rho, nspin )
  !         tmp_dir=tmp_dir_save
  !      ENDIF
  !   ENDIF
  !ENDDO
  RETURN
  END SUBROUTINE check_initial_status

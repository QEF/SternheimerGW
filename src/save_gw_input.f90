!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2017 Quantum ESPRESSO group,
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
MODULE save_gw
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data saved by the
  !     GW code to restart smoothly.
  !
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: restore_gw_input_variables, clean_input_variables
  !
  INTEGER, PRIVATE :: nat_todo_save, nrapp_save
  INTEGER, ALLOCATABLE, PRIVATE :: list_save(:), atomo_save(:) 
  CHARACTER(LEN=256), PUBLIC :: tmp_dir_save
  !
  !
  CONTAINS
    !
    SUBROUTINE restore_gw_input_variables(  )
      !------------------------------------------------------------------------
      !
      USE io_files,   ONLY : tmp_dir
      USE partial,    ONLY : nat_todo, nrapp
      !
      IMPLICIT NONE
      !
      nat_todo=nat_todo_save
      nrapp=nrapp_save
!      list=list_save
!      atomo=atomo_save
      tmp_dir=tmp_dir_save

      RETURN
    END SUBROUTINE restore_gw_input_variables

    SUBROUTINE clean_input_variables()
    IMPLICIT NONE
    RETURN
    END SUBROUTINE clean_input_variables
    !
END MODULE save_gw

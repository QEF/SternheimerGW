!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
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
SUBROUTINE sgw_opening_message()

    USE kinds, ONLY: DP
    USE io_global,  ONLY :  stdout, ionode, ionode_id
    USE sgw_version, ONLY: sgw_version_number, sgw_svn_revision

    IMPLICIT NONE

    WRITE( stdout, '(/5X, "You are using SGW version: ", A)') sgw_version_number
    WRITE( stdout, '(/5X, "With svn revision number: ", A)') sgw_svn_revision

    WRITE( stdout, '(/5X,"This program is part of the open-source Quantum ", &
         &  "ESPRESSO suite", &
         &/9X," Please also cite H. Lambert and F. Giustino, Phys. Rev. B ",&
         &    " 88. 075117 (2013);", &
         &/9X," URL http://www.sternheimer.org.uk")' )
    RETURN
END SUBROUTINE sgw_opening_message

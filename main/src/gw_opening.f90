!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
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
!> routines to print the SternheimerGW logo and an opening message
!------------------------------------------------------------------------------
MODULE gw_opening

  IMPLICIT NONE

  PRIVATE

  PUBLIC gw_opening_message, gw_opening_logo

CONTAINS

  !> Print the SternheimerGW logo.
  SUBROUTINE gw_opening_logo()

    USE global_version, ONLY: version_number, svn_revision
    USE gw_version,     ONLY: gw_version_number
    USE io_global,      ONLY: stdout, meta_ionode

    ! print the logo
    IF (meta_ionode) THEN
      WRITE(stdout, '(a)')
      WRITE(stdout, '(a)') '            ######## T E R N H E I M E R           http://www.sternheimergw.org/'
      WRITE(stdout, '(a)') '           ###     ###'
      WRITE(stdout, '(a)') '          ###       ##                  #######        #                      #'
      WRITE(stdout, '(a)') '          ###                         ##       ###     #          ##          #'
      WRITE(stdout, '(a)') '           ###                      ##            #     #         ##         #'
      WRITE(stdout, '(a)') '            #####                  #                     #       #  #       #'
      WRITE(stdout, '(a)') '   ///         #####      \\\     #                      #       #  #       #'
      WRITE(stdout, '(a)') '  ///             ####     \\\    #                       #     #    #     #'
      WRITE(stdout, '(a)') ' ///                ###     \\\   #         #########     #     #    #     #'
      WRITE(stdout, '(a)') '(((                  ###     )))  #                 #      #   #      #   #'
      WRITE(stdout, '(a)') ' \\\    ###          ###    ///    #                #      #   #      #   #'
      WRITE(stdout, '(a)') '  \\\    ###        ###    ///      #              #        # #        # #'
      WRITE(stdout, '(a)') '   \\\    ###      ###    ///        #           ##         ###        ###'
      WRITE(stdout, '(a)') '           ##########                 ###########            #          #'
      WRITE(stdout, '(a)')
    END IF

    ! adjust the version number so that the SternheimerGW version is printed
    version_number = gw_version_number
    svn_revision   = 'unknown'
    
  END SUBROUTINE gw_opening_logo

  !> Write the opening message.
  !!
  !! Print the version number used and the git revision if known.
  SUBROUTINE gw_opening_message()

    USE gw_version,  ONLY: gw_git_describe
    USE io_global,   ONLY: stdout

    IMPLICIT NONE

    WRITE(stdout,'(a)')
    WRITE(stdout,'(5x,a)') 'Please cite SternheimerGW as:'
    WRITE(stdout,'(9x,a)') 'M. Schlipf, H. Lambert, N. Zibouche, and F. Giustino, SternheimerGW,'
    WRITE(stdout,'(9x,a)') 'paper in preparation, URL http://www.sternheimergw.org'
    WRITE(stdout,'(a)')
    WRITE(stdout,'(5x,a)') 'To increase the reproducibility of your results you can mention the'
    WRITE(stdout,'(5x,a)') 'git description of this version: ' // gw_git_describe
    WRITE(stdout,'(a)')
    WRITE(stdout,'(5x,a)') 'other relevant papers:'
    WRITE(stdout,'(9x,a)') 'H. Lambert and F. Giustino, Phys. Rev. B 88, 075117 (2013)'
    WRITE(stdout,'(9x,a)') 'F. Giustino, M. L. Cohen, and S. G. Louie, Phys. Rev. B 81, 115105 (2010)'

  END SUBROUTINE gw_opening_message

END MODULE gw_opening

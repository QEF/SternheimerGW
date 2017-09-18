!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2017
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
!> Provides a wrapper around the various routines to perform an analytic continuation.
!!
!! Depending on the setting either a Godby-Needs PP model, a Pade approximant,
!! or the AAA algorithm is used to expand the given quantity from a few points
!! to the whole complex plane.
MODULE analytic_module

  IMPLICIT NONE

  !> Use the Godby-Needs plasmon pole model.
  INTEGER, PARAMETER :: godby_needs = 1

  !> Use the conventional Pade approximation.
  INTEGER, PARAMETER :: pade_approx = 2

  !> Use the robust Pade approximation.
  INTEGER, PARAMETER :: pade_robust = 3

  !> Use the AAA rational approximation.
  INTEGER, PARAMETER :: aaa_approx = 4

CONTAINS

END MODULE analytic_module

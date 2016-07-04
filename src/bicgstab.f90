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
!> Shifted BiCGStab(l) solver for linear equation \f$(A + \lambda I) x = b\f$.
!!
!! Implements the shifted BiCGStab(l) algorithm according to the paper of
!! Fromme, Computing **70**, 87 (2003). The general idea is that the matrix
!! \f$A\f$ and \f$A + \lambda I\f$ span the same Krylov subspace. Hence, we can
!! solve the linear equation of all linear problems at the cost of a single one.
!!
MODULE bicgstab_module

  IMPLICIT NONE

  PRIVATE

  PUBLIC bicgstab

CONTAINS

  !> Main driver routine of the algorithm.
  !!
  !! This subroutine implements the "Algorithm 2" of Fromme's paper.
  !!
  SUBROUTINE bicgstab

    !> Counter for number of iterations.
    INTEGER iter

    !> Maximum number of iterations.
    INTEGER, PARAMETER :: max_iter = 10000

    !
    ! initialization seed system
    !

    !
    ! initialization shifted systems
    !

    ! loop until solution is found
    DO iter = 1, max_iter

      !
      ! perform BiCG part (Algorithm 3)
      !

      !
      ! perform MR part (Algorithm 4)
      !

    END DO ! iter

    IF (iter > max_iter) THEN
      CALL errore(__FILE__, "BiCGstab algorithm did not converge in given&
                           & number of iterations", max_iter)
    END IF

  END SUBROUTINE bicgstab

END MODULE bicgstab_module

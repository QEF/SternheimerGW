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
!> This module provides the functionality to expand the igk_k array by one
!! element.
!!
!! This is necessary when the current k-point is not in the mesh of the
!! non-selfconsistent calculation.
!!
MODULE expand_igk_module

CONTAINS

  !> Extends the global igk_k array by one element.
  SUBROUTINE expand_igk

    USE klist, ONLY : igk_k, ngk, nks
    USE wvfct, ONLY : npwx

    !> temporary array storing the old igk_k
    INTEGER, ALLOCATABLE :: igk_k_tmp(:,:)

    !> temporary array storing the old ngk
    INTEGER, ALLOCATABLE :: ngk_tmp(:)

    ! trivial case - already appropriate size
    IF (SIZE(igk_k, 2) > nks) RETURN

    ! copy old igk and ngk to temporary array
    ALLOCATE(igk_k_tmp(npwx, nks))
    ALLOCATE(ngk_tmp(nks))
    igk_k_tmp = igk_k
    ngk_tmp = ngk

    ! reallocate igk_k and ngk increasing the size by one
    DEALLOCATE(igk_k)
    DEALLOCATE(ngk)
    ALLOCATE(igk_k(npwx, nks + 1))
    ALLOCATE(ngk(nks + 1))

    ! copy the temporary array in the new igk_k and ngk
    igk_k(:,:nks) = igk_k_tmp
    ngk(:nks) = ngk_tmp

    ! free the memory
    DEALLOCATE(igk_k_tmp)
    DEALLOCATE(ngk_tmp)

  END SUBROUTINE expand_igk

END MODULE expand_igk_module

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
!> Wrapper routine for PW’s gk_sort routine
!!
!! gk_sort of PW works well as long as the k-point mesh is consistent
!! throughout the calculation. However, for the k points not in the original
!! mesh of the ground-state calculation \f$|k + G| < E_{\text{cut}}\f$
!! might be fulfilled for more G vectors than for any k-points in the original
!! mesh. In this case, the original gk_sort would crash, because the allocated
!! array is to small. Solution: temporarily increase the array boundary for
!! gk_sort and decrease it afterwards.
!!
SUBROUTINE gk_sort_safe( k, ngm, g, ecut, ngk, igk, gk )

  USE kinds,     ONLY : dp
  USE wvfct,     ONLY : npwx

  IMPLICIT NONE

  REAL(dp), INTENT(IN)  :: k(3)      !< k point for which the G vectors are determined
  INTEGER,  INTENT(IN)  :: ngm       !< total number of G vectors
  REAL(dp), INTENT(IN)  :: g(3,ngm)  !< the 3 coordinates of the G vectors
  REAL(dp), INTENT(IN)  :: ecut      !< the energy cutoff \f$E_{\text{cut}}\f$
  INTEGER,  INTENT(OUT) :: ngk       !< number of G vectors for which \f$|k + G| < E_{\text{cut}}\f$
  INTEGER,  INTENT(OUT) :: igk(npwx) !< map between index in G vector array and sorted lengths of \f$|k + G\f$
  REAL(dp), INTENT(OUT) :: gk(npwx)  !< sorted lengths \f$|k + G|\f$ 

  INTEGER :: npwx_def

  INTEGER, ALLOCATABLE  :: igk_ext(:)
  REAL(dp), ALLOCATABLE :: gk_ext(:)

  ! temporarily increase array boundary
  npwx_def = npwx
  npwx = nint( npwx * 1.1 )

  ! allocate arrays with extended array size
  ALLOCATE( igk_ext(npwx), gk_ext(npwx) )

  ! sort the g-vectors with the increased array size
  CALL gk_sort( k, ngm, g, ecut, ngk, igk_ext, gk_ext )
  
  ! reset to the default
  npwx = npwx_def

  ! cut arrays to smaller size
  ngk = min( ngk, npwx )
  igk = igk_ext(:npwx)
  gk  = gk_ext(:npwx)

  ! deallocate arrays
  DEALLOCATE( igk_ext, gk_ext )

END SUBROUTINE gk_sort_safe

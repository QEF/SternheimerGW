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
SUBROUTINE rotate(xq, raq, s, nsym, isym)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE cell_base, ONLY : bg, at

  IMPLICIT NONE

  REAL(DP), PARAMETER  :: accep=1.e-5_dp
  REAL(DP), INTENT(IN) :: xq (3)
  INTEGER, INTENT(IN)  :: s (3, 3, 48), nsym
  REAL(DP)             :: wrk (3), aq (3), raq (3), zero (3)
  integer :: isym, ipol, jpol

  aq = xq
  call cryst_to_cart (1, aq, at, - 1)
  raq = 0.d0
  do ipol = 1, 3
     do jpol = 1, 3
        raq (ipol) = raq (ipol) + DBLE (s (ipol, jpol, isym) )*aq (jpol)
     enddo
  enddo
  call cryst_to_cart (1, raq, bg, +1)
END SUBROUTINE

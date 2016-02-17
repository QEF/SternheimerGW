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
LOGICAL ELEMENTAL FUNCTION testnan(x) RESULT (res)

  USE kinds, ONLY : DP

#if (defined(__INTEL) || __GNUC__ > 4)
  USE, INTRINSIC :: ieee_arithmetic
#endif

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: x

#if (defined(__INTEL) || __GNUC__ > 4)
  res = ieee_is_nan(x)
#else
#ifdef __GFORTRAN
  res = isnan(x)
#else
  res = (x /= x)
#endif
#endif

END FUNCTION testnan

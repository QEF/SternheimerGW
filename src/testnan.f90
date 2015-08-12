  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
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

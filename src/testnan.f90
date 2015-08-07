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

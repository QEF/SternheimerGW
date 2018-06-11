! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

MODULE sort_module

  USE constant_module, ONLY: dp

  IMPLICIT NONE

  PUBLIC insertion_sort
  PRIVATE

  INTERFACE swap
    MODULE PROCEDURE swap_real, swap_integer
  END INTERFACE swap
    
CONTAINS

  SUBROUTINE insertion_sort(array, map)
    !
    REAL(dp), INTENT(INOUT) :: array(:)
    INTEGER, INTENT(OUT), ALLOCATABLE, OPTIONAL :: map(:)
    INTEGER ii
    !
    CALL init_map_if_present(SIZE(array), map)
    DO ii = 2, SIZE(array)
      CALL swap_from_back_until_sorted(array(1:ii), map)
    END DO
    !
  END SUBROUTINE insertion_sort

  SUBROUTINE init_map_if_present(size, map)
    !
    INTEGER, INTENT(IN) :: size
    INTEGER, INTENT(OUT), ALLOCATABLE, OPTIONAL :: map(:)
    INTEGER ii
    !
    IF (PRESENT(map)) THEN
      ALLOCATE(map(size))
      map = [(ii, ii = 1, size)]
    END IF
    !
  END SUBROUTINE init_map_if_present

  SUBROUTINE swap_from_back_until_sorted(array, map)
    !
    REAL(dp), INTENT(INOUT) :: array(:)
    INTEGER, INTENT(INOUT), OPTIONAL :: map(:)
    INTEGER ii
    !
    DO ii = SIZE(array), 2, -1
      IF (sorted(array(ii-1:ii))) EXIT
      CALL swap(array(ii-1:ii))
      IF (PRESENT(map)) CALL swap(map(ii-1:ii))
    END DO
    !
  END SUBROUTINE swap_from_back_until_sorted 

  LOGICAL FUNCTION sorted(two_value)
    !
    REAL(dp), INTENT(IN) :: two_value(2)
    sorted = two_value(1) < two_value(2)
    !
  END FUNCTION sorted

  SUBROUTINE swap_real(two_value)
    !
    REAL(dp), INTENT(INOUT) :: two_value(2)
    REAL(dp) tmp
    !
    tmp = two_value(1)
    two_value(1) = two_value(2)
    two_value(2) = tmp
    !
  END SUBROUTINE swap_real

  SUBROUTINE swap_integer(two_value)
    !
    INTEGER, INTENT(INOUT) :: two_value(2)
    INTEGER tmp
    !
    tmp = two_value(1)
    two_value(1) = two_value(2)
    two_value(2) = tmp
    !
  END SUBROUTINE swap_integer

END MODULE sort_module

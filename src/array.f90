! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

MODULE array_module

  IMPLICIT NONE

  INTERFACE allocate_copy_from_to
    MODULE PROCEDURE allocate_copy_real, allocate_copy_real_2d, allocate_copy_real_3d, &
      allocate_copy_complex, allocate_copy_complex_2d, allocate_copy_complex_3d, &
      allocate_copy_integer, allocate_copy_integer_2d, allocate_copy_integer_3d
  END INTERFACE

  INTERFACE allocate_init_to
    MODULE PROCEDURE allocate_init_real, allocate_init_real_2d, allocate_init_real_3d, &
      allocate_init_complex, allocate_init_complex_2d, allocate_init_complex_3d, &
      allocate_init_logical
  END INTERFACE allocate_init_to

  INTEGER, PARAMETER :: no_error = 0
  INTEGER, PARAMETER :: array_size_error = 1

  PRIVATE
  PUBLIC allocate_copy_from_to, allocate_init_to, &
    no_error, array_size_error

CONTAINS

  SUBROUTINE allocate_copy_real(source, destination)

    USE constant_module, ONLY: dp

    REAL(dp), INTENT(IN) :: source(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: destination(:)

    ALLOCATE(destination(SIZE(source)))
    destination = source

  END SUBROUTINE allocate_copy_real

  SUBROUTINE allocate_copy_real_2d(source, destination)

    USE constant_module, ONLY: dp

    REAL(dp), INTENT(IN) :: source(:,:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: destination(:,:)

    ALLOCATE(destination(SIZE(source, 1), SIZE(source, 2)))
    destination = source

  END SUBROUTINE allocate_copy_real_2d

  SUBROUTINE allocate_copy_real_3d(source, destination)

    USE constant_module, ONLY: dp

    REAL(dp), INTENT(IN) :: source(:,:,:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: destination(:,:,:)

    ALLOCATE(destination(SIZE(source, 1), SIZE(source, 2), SIZE(source, 3)))
    destination = source

  END SUBROUTINE allocate_copy_real_3d

  SUBROUTINE allocate_copy_complex(source, destination)

    USE constant_module, ONLY: dp

    COMPLEX(dp), INTENT(IN) :: source(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: destination(:)

    ALLOCATE(destination(SIZE(source)))
    destination = source

  END SUBROUTINE allocate_copy_complex

  SUBROUTINE allocate_copy_complex_2d(source, destination)

    USE constant_module, ONLY: dp

    COMPLEX(dp), INTENT(IN) :: source(:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: destination(:,:)

    ALLOCATE(destination(SIZE(source, 1), SIZE(source, 2)))
    destination = source

  END SUBROUTINE allocate_copy_complex_2d

  SUBROUTINE allocate_copy_complex_3d(source, destination)

    USE constant_module, ONLY: dp

    COMPLEX(dp), INTENT(IN) :: source(:,:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: destination(:,:,:)

    ALLOCATE(destination(SIZE(source, 1), SIZE(source, 2), SIZE(source, 3)))
    destination = source

  END SUBROUTINE allocate_copy_complex_3d

  SUBROUTINE allocate_copy_integer(source, destination)

    INTEGER, INTENT(IN) :: source(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: destination(:)

    ALLOCATE(destination(SIZE(source)))
    destination = source

  END SUBROUTINE allocate_copy_integer

  SUBROUTINE allocate_copy_integer_2d(source, destination)

    INTEGER, INTENT(IN) :: source(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: destination(:,:)

    ALLOCATE(destination(SIZE(source, 1), SIZE(source, 2)))
    destination = source

  END SUBROUTINE allocate_copy_integer_2d

  SUBROUTINE allocate_copy_integer_3d(source, destination)

    INTEGER, INTENT(IN) :: source(:,:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: destination(:,:,:)

    ALLOCATE(destination(SIZE(source, 1), SIZE(source, 2), SIZE(source, 3)))
    destination = source

  END SUBROUTINE allocate_copy_integer_3d

  SUBROUTINE allocate_init_real(array_size, init_value, array)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size
    REAL(dp), INTENT(IN) :: init_value
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: array(:)

    ALLOCATE(array(array_size))
    array = init_value

  END SUBROUTINE allocate_init_real

  SUBROUTINE allocate_init_real_2d(array_size, init_value, array, info)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size(:)
    REAL(dp), INTENT(IN) :: init_value
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: array(:,:)
    INTEGER, INTENT(OUT), OPTIONAL :: info

    IF (array_size_incorrect(2, SIZE(array_size), info)) RETURN

    ALLOCATE(array(array_size(1), array_size(2)))
    array = init_value

  END SUBROUTINE allocate_init_real_2d

  SUBROUTINE allocate_init_real_3d(array_size, init_value, array, info)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size(:)
    REAL(dp), INTENT(IN) :: init_value
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: array(:,:,:)
    INTEGER, INTENT(OUT), OPTIONAL :: info

    IF (array_size_incorrect(3, SIZE(array_size), info)) RETURN

    ALLOCATE(array(array_size(1), array_size(2), array_size(3)))
    array = init_value

  END SUBROUTINE allocate_init_real_3d

  SUBROUTINE allocate_init_complex(array_size, init_value, array)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size
    COMPLEX(dp), INTENT(IN) :: init_value
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: array(:)

    ALLOCATE(array(array_size))
    array = init_value

  END SUBROUTINE allocate_init_complex

  SUBROUTINE allocate_init_complex_2d(array_size, init_value, array, info)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size(:)
    COMPLEX(dp), INTENT(IN) :: init_value
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: array(:,:)
    INTEGER, INTENT(OUT), OPTIONAL :: info

    IF (array_size_incorrect(2, SIZE(array_size), info)) RETURN

    ALLOCATE(array(array_size(1), array_size(2)))
    array = init_value

  END SUBROUTINE allocate_init_complex_2d

  SUBROUTINE allocate_init_complex_3d(array_size, init_value, array, info)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size(:)
    COMPLEX(dp), INTENT(IN) :: init_value
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: array(:,:,:)
    INTEGER, INTENT(OUT), OPTIONAL :: info

    IF (array_size_incorrect(3, SIZE(array_size), info)) RETURN

    ALLOCATE(array(array_size(1), array_size(2), array_size(3)))
    array = init_value

  END SUBROUTINE allocate_init_complex_3d

  SUBROUTINE allocate_init_logical(array_size, init_value, array)

    USE constant_module, ONLY: dp

    INTEGER, INTENT(IN) :: array_size
    LOGICAL, INTENT(IN) :: init_value
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: array(:)

    ALLOCATE(array(array_size))
    array = init_value

  END SUBROUTINE allocate_init_logical

  LOGICAL FUNCTION array_size_incorrect(desired, actual, info)

    USE assert_module, ONLY: assert

    INTEGER, INTENT(IN) :: desired, actual
    INTEGER, INTENT(OUT), OPTIONAL :: info

    array_size_incorrect = desired /= actual

    IF (PRESENT(INFO)) THEN
      IF (array_size_incorrect) THEN
        info = array_size_error
      ELSE
        info = no_error
      END IF

    ELSE
      CALL assert(.NOT.array_size_incorrect, &
                  "allocate_init: dimension of array_size inconsistent with array shape")
    END IF

  END FUNCTION array_size_incorrect

END MODULE array_module

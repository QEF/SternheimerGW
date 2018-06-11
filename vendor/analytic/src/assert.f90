! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

MODULE assert_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC assert

CONTAINS

  SUBROUTINE assert(test, message)

    USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit

    LOGICAL, INTENT(IN) :: test
    CHARACTER(LEN=*), INTENT(IN) :: message

#if defined(__DEBUG)

    IF (test) RETURN 

    WRITE(error_unit,*) message
    CALL abort()

#endif

  END SUBROUTINE assert

END MODULE assert_module

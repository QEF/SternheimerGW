! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

PROGRAM test_assert

  USE assert_module, ONLY: assert

  IMPLICIT NONE

  CALL assert(.TRUE., "This test should pass.")
  CALL assert(.FALSE., "Intentionally broken.")

END PROGRAM test_assert

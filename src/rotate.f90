  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
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
  call cryst_to_cart (1, raq, bg,  1)
END SUBROUTINE

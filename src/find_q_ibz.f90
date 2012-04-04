SUBROUTINE find_q_ibz(xq_ibk, s, iq, isym)
!Routine finds which symmetry operation folds xq to the required qpoint, xq_req.
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, nqs

  IMPLICIT NONE

  REAL(DP)                :: xq_ibz(3)
  REAL(DP)                :: x_q_loc(3)
  REAL(DP)                :: xq_ibk_locr(3)
  REAL(DP), INTENT(IN)    :: xq_ibk(3)
  REAL(DP)                :: xq_ibk_loc(3)

  !_loc so i don't operate on arrays passed to subroutine
  !that would require rotating the vectors twice each time. 

  LOGICAL                 :: found_q
  INTEGER,  INTENT(OUT) :: iq, isym
  INTEGER                 :: s(3,3,48), invs(48), ism1
  INTEGER                 :: i, j, k
  REAL(DP), PARAMETER     :: eps=1.0d-5


  !The logic of this routine is: 
  !xq_ibk_locr is the symmetry rotated xq_ibk 
  !x_q_loc is the q point in the IBZ in crystal co-ordinates.
  !if x_q_loc = xq_ibk_locr then we have found the symmetry operation
  !which rotates x_q_ibk to i_q_ibz. 
  !R(q_{IBK}) = xq_{IBZ}

  !Transform xq into crystal co-ordinates. 
   xq_ibk_loc(:) = xq_ibk 
   CALL cryst_to_cart(1, xq_ibk_loc(:), at, -1)
   found_q=.false.
DO iq = 1, nqs
  !Transform xq into crystal co-ordinates. 
    x_q_loc(:) = x_q(:,iq)
   CALL cryst_to_cart(1, x_q_loc(:), at, -1)

   DO isym = 1, 48
      ism1 = invs(isym)
      xq_ibk_locr(1) = s (1, 1, isym) * xq_ibk_loc(1)  + s (1, 2, isym) * xq_ibk_loc(2) + s (1, 3, isym) * xq_ibk_loc(3)
      xq_ibk_locr(2) = s (2, 1, isym) * xq_ibk_loc(1)  + s (2, 2, isym) * xq_ibk_loc(2) + s (2, 3, isym) * xq_ibk_loc(3)
      xq_ibk_locr(3) = s (3, 1, isym) * xq_ibk_loc(1)  + s (3, 2, isym) * xq_ibk_loc(2) + s (3, 3, isym) * xq_ibk_loc(3)
 
      found_q  = (abs(x_q_loc(1) - xq_ibk_locr(1)).lt.eps).and. &
                 (abs(x_q_loc(2) - xq_ibk_locr(2)).lt.eps).and. & 
                 (abs(x_q_loc(3) - xq_ibk_locr(3)).lt.eps) 
      if (found_q) return
   END DO
END DO

if (.not.found_q) then 
    write(6,'("q_point not found in IBZ.")')
    STOP
endif
END SUBROUTINE

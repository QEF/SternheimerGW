SUBROUTINE find_q_ibz(xq_ibk, s, iq, isym, found_q)
! Routine finds which symmetry operation folds xq to the required qpoint, xq_req.
! While it is true we can use the full symmetry group of the crystal for W
! regardless of still need to use only valid symops 

!For crystals without inversion symmetry we pray for a symmetry operations  Sq -> -q +G

  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, nqs
  USE modes,         ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
 !  7) computes the variables needed to pass to the pattern representation
 !     tmq    the matrix of the symmetry which sends q -> -q + G
 !     gi     the G associated to each symmetry operation
 !     gimq   the G of the q -> -q+G symmetry
 !     irgq   the small group indices
 !     nsymq  the order of the small group of q
 !     irotmq the index of the q->-q+G symmetry
 !     nirr   the number of irreducible representation
 !     npert  the dimension of each irreducible representation
 !     nmodes the number of modes
 !     minus_q true if there is a symmetry sending q -> -q+G
  USE symm_base,     ONLY : t_rev, irt, ftau, nrot, nsym, &
                            time_reversal, sname, d1, d2, d3, &
                            copy_sym, s_axis_to_cart
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime, npool, &
                            nproc_pool, me_pool, my_pool_id, nproc

  IMPLICIT NONE

  REAL(DP)                :: xq_ibz(3)
  REAL(DP)                :: x_q_loc(3)
  REAL(DP)                :: xq_ibk_locr(3)
  REAL(DP), INTENT(IN)    :: xq_ibk(3)
  REAL(DP)                :: xq_ibk_loc(3)

  !_loc so i don't operate on arrays passed to subroutine
  !that would require rotating the vectors twice each time. 

  LOGICAL                 :: found_q, s_minus_q
 !variable to acknowledge we've rotated q back the IBZ
 !dummy variable for routine to find Sq = -q + G

  INTEGER,  INTENT(OUT)   :: iq, isym
  INTEGER                 :: s(3,3,48), invs(48), ism1
  INTEGER                 :: i, j, k
  REAL(DP), PARAMETER     :: eps=1.0d-5


  !The logic of this routine is: 
  !xq_ibk_locr is the symmetry rotated xq_ibk 
  !x_q_loc is the q point in the IBZ in crystal co-ordinates.
  !if x_q_loc = xq_ibk_locr then we have found the symmetry operation
  !which rotates x_q_ibk to i_q_ibz. 
  !R(q_{IBK}) = xq_{IBZ}
  !q_{IBK} = R^{-1} xq_{IBZ}
  !Hit the vector we need with every symm op R until it equals a q point in the IBZ.

  !Possible solutions to not using symmetry
  !1) Give up, join a monastery and commit myself to a life of good deeds (simplest). 
  !Transform xq into cartesian co-ordinates. 

   xq_ibk_loc(:) = xq_ibk 
   CALL cryst_to_cart(1, xq_ibk_loc(:), at, -1)
   found_q=.false.

DO iq = 1, nqs
  !Transform xq into cartesian co-ordinates. 
   x_q_loc(:) = x_q(:,iq)
   CALL cryst_to_cart(1, x_q_loc(:), at, -1)
   DO isym = 1, nsym
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

if ((.not.found_q).and.(.not.invsymq)) then 
   !Search for a symmop that gives Sq = -q+G.
   ! xq_ibk_loc(:) = xq_ibk 
   ! s_minus_q=.true.
   !CALL sgq(xq_ibk_loc, at, bg, s, nsym, irgq, nsymq, irotmq, s_minus_q, gi, gimq)
    write(1000+mpime,'("q_point not found in IBZ.")')
   !CALL mp_global_end()
   !STOP
endif

RETURN
END SUBROUTINE

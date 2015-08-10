SUBROUTINE find_qg_ibz(xq_ibk, s, iq, isym, nig0, found_q, inv_q)
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE klist,         ONLY : nks, nkstot, ngauss, degauss, xk, wk
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
  USE gvect,          ONLY : g

  IMPLICIT NONE

  REAL(DP)                :: xq_ibz(3)
  REAL(DP)                :: x_q_loc(3)
  REAL(DP)                :: xq_ibk_locr(3)
  REAL(DP), INTENT(IN)    :: xq_ibk(3)
  REAL(DP)                :: xq_ibk_loc(3)
  !_loc so i don't operate on arrays passed to subroutine
  !that would require rotating the vectors twice each time. 
  LOGICAL                 :: found_q, s_minus_q, inv_q
 !variable to acknowledge we've rotated q back the IBZ
 !dummy variable for routine to find Sq = -q + G
  INTEGER,  INTENT(OUT)   :: iq, isym
  INTEGER,  INTENT(INOUT) :: nig0
 !integer of plane wave that lets us map \psi_{k-\q} to \psi_{\k-\q+G_0\in IBZ} 
  INTEGER                 :: s(3,3,48), invs(48), ism1
  INTEGER                 :: i, j, k, invsym, ig
  REAL(DP)                :: g_loc(3)
  REAL(DP), PARAMETER     :: eps=1.0d-5
  !For COULMAT:
  !K_{IBZ}+G = SK'
   xq_ibk_loc(:) = xq_ibk 
  !Transform xq into crystal co-ordinates. 
   CALL cryst_to_cart(1, xq_ibk_loc(:), at, -1)
   found_q=.false.
   nig0 = 0
   inv_q =.false.
   !DO iq = 1, nks
   DO iq = 1, nkstot
!no pooling!
      x_q_loc(:) = xk(:,iq)
      CALL cryst_to_cart(1, x_q_loc(:), at, -1)
    !!DO ig = 1 , 27 WORKS FOR CUBIC CRYSTAL
      DO ig = 1 , 125 !CAC6 IS very freaky c.f. FG's cuprates stuff
          g_loc = g(:,ig)
          call cryst_to_cart(1, g_loc, at, -1)
        DO isym = 1, nsym
           ism1 = invs(isym)
           xq_ibk_locr(1) = s (1, 1, isym) * (xq_ibk_loc(1) + g_loc(1)) + s (1, 2, isym) * (xq_ibk_loc(2) + g_loc(2))  + s (1, 3, isym) * (xq_ibk_loc(3) + g_loc(3))
           xq_ibk_locr(2) = s (2, 1, isym) * (xq_ibk_loc(1) + g_loc(1)) + s (2, 2, isym) * (xq_ibk_loc(2) + g_loc(2))  + s (2, 3, isym) * (xq_ibk_loc(3) + g_loc(3))
           xq_ibk_locr(3) = s (3, 1, isym) * (xq_ibk_loc(1) + g_loc(1)) + s (3, 2, isym) * (xq_ibk_loc(2) + g_loc(2))  + s (3, 3, isym) * (xq_ibk_loc(3) + g_loc(3))
           found_q  = (abs(x_q_loc(1) - xq_ibk_locr(1)).le.eps).and. &
                      (abs(x_q_loc(2) - xq_ibk_locr(2)).le.eps).and. & 
                      (abs(x_q_loc(3) - xq_ibk_locr(3)).le.eps) 
           if (found_q) then
               nig0  = ig
               inv_q = .false.
               return
           endif
        END DO
      END DO
   END DO

!If crystal lacks inversion symmetry we find G_{-k} using time reversal
!symmetry.
!IF (.not.found_q) then 
!   xq_ibk_loc(:) = -xq_ibk 
!  !Transform xq into crystal co-ordinates. 
!   CALL cryst_to_cart(1, xq_ibk_loc(:), at, -1)
!   DO iq = 1, nkstot
!      !x_q_loc(:) = xk(:,iq)
!!Seeing if we can use time reversal
!      x_q_loc(:) = x_q(:,iq)
!      CALL cryst_to_cart(1, x_q_loc(:), at, -1)
!      DO ig = 1 , 125 !CAC6 IS very freaky c.f. FG's cuprates stuff
!         g_loc = g(:,ig)
!         call cryst_to_cart(1, g_loc, at, -1)
!         DO isym = 1, nsym
!            ism1 = invs(isym)
!            xq_ibk_locr(1) = s (1, 1, isym) * (xq_ibk_loc(1) + g_loc(1)) + s (1, 2, isym) * (xq_ibk_loc(2) + g_loc(2))  + s (1, 3, isym) * (xq_ibk_loc(3) + g_loc(3))
!            xq_ibk_locr(2) = s (2, 1, isym) * (xq_ibk_loc(1) + g_loc(1)) + s (2, 2, isym) * (xq_ibk_loc(2) + g_loc(2))  + s (2, 3, isym) * (xq_ibk_loc(3) + g_loc(3))
!            xq_ibk_locr(3) = s (3, 1, isym) * (xq_ibk_loc(1) + g_loc(1)) + s (3, 2, isym) * (xq_ibk_loc(2) + g_loc(2))  + s (3, 3, isym) * (xq_ibk_loc(3) + g_loc(3))
!            found_q  = (abs(x_q_loc(1) - xq_ibk_locr(1)).le.eps).and. &
!                       (abs(x_q_loc(2) - xq_ibk_locr(2)).le.eps).and. & 
!                       (abs(x_q_loc(3) - xq_ibk_locr(3)).le.eps) 
!            if (found_q) then
!                nig0 = ig
!                inv_q = .true.
!                return
!            endif
!         END DO
!      END DO
!   END DO
!ENDIF
RETURN
END SUBROUTINE find_qg_ibz

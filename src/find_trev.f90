  !-----------------------------------------------------------------------!
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino             !
  ! This file is distributed under the terms of the GNU General Public    !     
  ! License. See the file `LICENSE' in the root directory of the          !     
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .        !
  !-----------------------------------------------------------------------!
subroutine find_trev(xq_ibk, s, invs, iqtr, isymcoul, trev)
  use kinds,         only : DP
  use io_global,     only : stdout, ionode_id, ionode, meta_ionode
  use cell_base,     only : at, bg
  use disp,          only : nqs, nq1, nq2, nq3, wq, x_q
  use lr_symm_base,  only : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  use symm_base,     only : t_rev, irt, ftau, nrot, nsym, &
                            time_reversal, sname, d1, d2, d3, &
                            copy_sym, s_axis_to_cart
  implicit none

  real(DP)                :: xq_ibz(3)
  real(DP)                :: x_q_loc(3)
  real(DP)                :: xq_ibk_locr(3)
  real(DP), intent(in)    :: xq_ibk(3)
  real(DP)                :: xq_ibk_loc(3)
  real(DP), parameter     :: eps=1.0d-5
!_loc so i don't operate on arrays passed to subroutine
!that would require rotating the vectors twice each time. 
  logical                 :: found_q, s_minus_q, inv_q, trev
!variable to acknowledge we've rotated q back the IBZ
!dummy variable for routine to find Sq = -q + G
  integer                 :: iqtr, isymcoul
  integer                 :: s(3,3,48), invs(48), ism1
  integer                 :: i, j, k, invsym, iq, isym, ni

!This routine finds index of 
!a q-point related via time 
!reversal symmetry.
   xq_ibk_loc(:) = xq_ibk(:)
   call cryst_to_cart(1, xq_ibk_loc(:), at, -1)
   found_q=.false.
   isymcoul = 1
!write(1000, '(5x,i3, 4f14.9)') iq, xq_ibk_loc(1), xq_ibk_loc(2), xq_ibk_loc(3)
   do iq = 1, nqs
      x_q_loc(:) = x_q(:,iq)
      call cryst_to_cart(1, x_q_loc(:), at, -1)
      found_q  = ( abs(xq_ibk_loc(1) - x_q_loc(1)).le.eps).and. &
                 ( abs(xq_ibk_loc(2) - x_q_loc(2)).le.eps).and. & 
                 ( abs(xq_ibk_loc(3) - x_q_loc(3)).le.eps) 
      if (found_q) then 
          iqtr = iq
          trev = .false.
          isymcoul = 1
          exit
      endif
   end do

   if(.not. found_q) then
      xq_ibk_loc(:) = -xq_ibk 
      call cryst_to_cart(1, xq_ibk_loc(:), at, -1)
      write(6, '(5x, "xk_ibk ", 3f14.9)') xq_ibk_loc(1), xq_ibk_loc(2), xq_ibk_loc(3)
      do iq = 1, nqs
         x_q_loc(:) = x_q(:,iq)
         call cryst_to_cart(1, x_q_loc(:), at, -1)
         do isym = 1, nsym
           xq_ibk_locr(1)= s(1,1,isym)*(xq_ibk_loc(1))+s(1,2,isym)*(xq_ibk_loc(2))+s(1,3,isym)*(xq_ibk_loc(3))
           xq_ibk_locr(2)= s(2,1,isym)*(xq_ibk_loc(1))+s(2,2,isym)*(xq_ibk_loc(2))+s(2,3,isym)*(xq_ibk_loc(3))
           xq_ibk_locr(3)= s(3,1,isym)*(xq_ibk_loc(1))+s(3,2,isym)*(xq_ibk_loc(2))+s(3,3,isym)*(xq_ibk_loc(3))
           do i = -1, 1
            do j = -1, 1
             do k = -1, 1
                found_q  = (abs(x_q_loc(1)-(xq_ibk_locr(1)+float(i))).le.eps).and. &
                           (abs(x_q_loc(2)-(xq_ibk_locr(2)+float(j))).le.eps).and. & 
                           (abs(x_q_loc(3)-(xq_ibk_locr(3)+float(k))).le.eps) 
                if(found_q) then
                   iqtr = iq
                   isymcoul = isym
                   trev = .true.
 !                 write(6, '(5x, "TR", 3i3, 3f14.9)') i,j,k, isym, x_q_loc(1), x_q_loc(2), x_q_loc(3)
                   GOTO 125
                endif
             enddo
            enddo
           enddo
        enddo
      enddo
  125 CONTINUE
  endif
if (.not.found_q) CALL errore( 'find_trev', 'cant find qpoint in IBZ', 1 )
!if (.not.found_q) write(stdout, '("Cant find qpoint")')
return
end subroutine

!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
subroutine find_xk_unique(xq, xk_kpoints, xk1q, nsq, numsq, wgt)
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE disp,          ONLY : x_q,nqs
  USE symm_base,     ONLY : t_rev, irt, ftau, nrot, nsym, &
                            time_reversal, sname, d1, d2, d3, &
                            copy_sym, s_axis_to_cart, invs,s
  implicit none

  REAL(DP) :: xk_kpoints(3,1), xq(3), xktmp(3)
  real(DP), parameter :: accep=1.e-5_dp
! a zero vector: used in eqvect
  real(DP) :: saq (3, 48), aq (3), raq (3), zero (3)
  real(DP) :: xk1q(3, nsym)
  INTEGER :: numsq, isym, i, counter
  INTEGER :: nsq(48)
  INTEGER :: equiv(nsym)
  REAL(DP) :: wgt(nsym)

  logical, external :: eqvect
  LOGICAL :: unique
  zero(:) = 0.d0
!generate all the xkpoints
     numsq = 1 
     xk1q(:,:) = 0.0
     wgt(:) = 1.0
     equiv(:) = 0
     do isym = 1, nsym
       CALL rotate(xq, aq, s, nsym, invs(isym))
       xktmp  = xk_kpoints(:,1) - aq(:)
       if(numsq.eq.1.and.isym.eq.1) xk1q(:,numsq) = xktmp
       if(numsq.eq.1.and.isym.eq.1) nsq(numsq) = isym

       if(isym.gt.1) then 
         do i = 1, numsq
            if (eqvect (xk1q(:,i), xktmp(1), zero, accep) ) then
              equiv(i) = 1 
              wgt(i)   = wgt(i) + 1.0
            else
              equiv(i) = 0 
             endif
         enddo
         unique = .true.
         do i= 1, numsq
           unique = unique.and.(equiv(i).eq.0)
         enddo
         if(unique) then
            numsq = numsq +1
            xk1q(:, numsq) = xktmp
            nsq(numsq) = isym
         endif
       endif
     enddo
end subroutine find_xk_unique

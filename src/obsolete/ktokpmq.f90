!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 Jesse Noffsinger, Brad Malone,
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
  subroutine ktokpmq ( xk0, xq0, sign, nkq)
  !--------------------------------------------------------
  !
  !   For a given k point in cart coord, find the index 
  !   of the corresponding (k + sign*q) point
  !
  !--------------------------------------------------------

  USE kinds,         ONLY : DP
  USE constants,     ONLY :  eps8
  USE disp,          ONLY : nq1, nq2, nq3, nqs
  USE cell_base,     ONLY : at

  implicit none

  !
  real(kind=DP) :: xk0 (3), xq0 (3)
  ! input: coordinates of k points and q points
  integer :: sign, ipool, nkq, nkq_abs
  ! input: +1 for searching k+q, -1 for k-q
  ! output: in the parallel case, the pool hosting the k+-q point    
  ! output: the index of k+sign*q
  ! output: the absolute index of k+sign*q (in the full k grid)
  ! work variables
  !

  real(kind=DP) :: xxk (3), xxq (3)
  integer :: nkl, nkbl, nkr, iks, ik, i, j, k, n, jpool
  real(kind=DP) :: xx, yy, zz
  logical :: in_the_list
  !
  !if (abs(sign).ne.1) call error ('ktokpmq','sign must be +1 or -1',1)
  !
  ! bring k and q in crystal coordinates
  !
  xxk = xk0
  xxq = xq0
  call cryst_to_cart (1, xxk, at, -1)
  call cryst_to_cart (1, xxq, at, -1)
  !
  !  check that k is actually on a uniform mesh centered at gamma
  !
  xx = xxk(1)*nq1
  yy = xxk(2)*nq2
  zz = xxk(3)*nq3

! HL switching from eps = ^-10 to ^-8.

  in_the_list = abs(xx-nint(xx)).le.eps8 .and. &
                abs(yy-nint(yy)).le.eps8 .and. &
                abs(zz-nint(zz)).le.eps8
 ! if (.not.in_the_list) call error ('ktokpmq','is this a uniform k-mesh?',1)
 
 ! HL: now add the coulomb perturbing wavevector and check that k+q falls again on the k grid.
 
  xxk = xxk + float(sign) * xxq
 !
  xx = xxk(1)*nq1
  yy = xxk(2)*nq2
  zz = xxk(3)*nq3
  in_the_list = abs(xx-nint(xx)).le.eps8 .and. &
                abs(yy-nint(yy)).le.eps8 .and. &
                abs(zz-nint(zz)).le.eps8
! if (.not.in_the_list) call error ('ktokpmq','k+q does not fall on k-grid',1)
  !
  !  find the index of this k+q in the k-grid
  !
  i = mod ( nint ( xx + 2*nq1), nq1 )
  j = mod ( nint ( yy + 2*nq2), nq2 )
  k = mod ( nint ( zz + 2*nq3), nq3 )
  n = i*nq2*nq3 + j*nq3 + k + 1
  !
  nkq = n
  !  now n represents the index of k+sign*q in the original k grid.
  end subroutine ktokpmq
  !--------------------------------------------------------

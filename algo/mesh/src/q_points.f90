!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
subroutine q_points ( )
!---------------------------------------------------------

  use cell_base,  only : bg
  use disp,       only : nq1, nq2, nq3, x_q, nqs, wq, iq1, iq2, iq3
  use io_global,  only : stdout
  use kinds,      only : dp
  use mp,         only : mp_bcast
  use symm_base,  only : nsym, s, time_reversal, t_rev

  implicit none

  integer :: i, iq
  logical :: exist_gamma, skip_equivalence=.FALSE.
  logical, external :: check_q_points_sym
  real(DP), allocatable :: xq(:,:)

  INTEGER :: nqmax
  !
  !  calculate the Monkhorst-Pack grid
  !

  if( nq1 <= 0 .or. nq2 <= 0 .or. nq3 <= 0 ) &
       call errore('q_points','nq1 or nq2 or nq3 <= 0',1)

  nqmax= nq1 * nq2 * nq3

  allocate (wq(nqmax))
  allocate (xq(3,nqmax))
  call kpoint_grid(nsym, time_reversal, skip_equivalence, s, t_rev,& 
                      bg, nqmax, iq1,iq2,iq3, nq1,nq2,nq3, nqs, xq, wq )
                         
  !so that it is equivalent with set_defaults_pw.f90
  !call kpoint_grid( nrot, time_reversal, .false., s, t_rev, bg, nqmax,&
  !                        0,0,0, nq1,nq2,nq3, nqs, xq, wq )
  allocate(x_q(3,nqs))
  x_q(:,:) = xq(:,1:nqs)
  deallocate (xq)
  !
  ! Check if the Gamma point is one of the points and put
  ! it in the first position (it should already be the first)
  !
  exist_gamma = .false.
  do iq = 1, nqs
     if ( abs(x_q(1,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(2,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(3,iq)) .lt. 1.0e-10_dp ) then
        exist_gamma = .true.
        if (iq .ne. 1) then
           do i = 1, 3
              x_q(i,iq) = x_q(i,1)
              x_q(i,1) = 0.0_dp
           end do
        end if
     end if
  end do
  !
  ! Write the q points in the output
  !
  write(stdout, '(//5x,"Coulomb Perturbations for (", 3(i2,","),") &
           & uniform grid of q-points")') nq1, nq2, nq3
  write(stdout, '(5x,"(",i4,"q-points):")') nqs
  write(stdout, '(5x,"  N         xq(1)         xq(2)         xq(3) " )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 3f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
  end do
  !
  IF ( .NOT. exist_gamma) &
     CALL errore('q_points','Gamma is not a q point',1)
!
!  Check that the q point grid is compatible with the symmetry.
!  If this test is not passed, q2r will stop in any case.
!
  return
end subroutine q_points
!

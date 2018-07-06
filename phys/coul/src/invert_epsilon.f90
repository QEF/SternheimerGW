!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
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
SUBROUTINE invert_epsilon(num_g_corr, scrcoul_g_in, lgamma)

USE freq_gw,       ONLY : nfs
USE kinds,         ONLY : DP

IMPLICIT NONE    

!> the number of G vectors in the correlation grid
INTEGER, INTENT(IN) :: num_g_corr

COMPLEX(DP)       :: scrcoul_g_in(num_g_corr, num_g_corr, nfs, 1)
COMPLEX(DP)       :: work(num_g_corr)
INTEGER           :: ig, igp
INTEGER           :: iw
INTEGER           :: iwork(num_g_corr), info
LOGICAL           :: lgamma

!> complex constant of 0
COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

!> complex constant of 1
COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

!at Gamma wings of \Chi are 0.
if(lgamma) then
  do iw = 1, nfs
    do ig = 2, num_g_corr
       scrcoul_g_in(ig,1,iw,1)  = zero 
    enddo
    do igp = 2, num_g_corr
       scrcoul_g_in(1,igp,iw,1) = zero
    enddo
  enddo
endif

!Need block inversion routine if iq is gamma.
do iw = 1, nfs
   call ZGETRF (num_g_corr, num_g_corr,&
   scrcoul_g_in(1:num_g_corr,1:num_g_corr,iw,1), num_g_corr, iwork, info)
   call errore ('invert epsilon', 'factorization', info)
   call ZGETRI (num_g_corr, scrcoul_g_in(1:num_g_corr,1:num_g_corr,iw,1),& 
   num_g_corr, iwork, work, num_g_corr, info)
   call errore ('invert epsilon', 'inversion', info)
enddo

write(6,*)
write(6,'(5x, "Done epsilon inversion.")') 
write(6,'(5x, "")') 

if(lgamma) then
  do iw = 1, nfs
     do ig = 2, num_g_corr
        scrcoul_g_in(ig,1,iw,1) = zero
     enddo
     do igp = 2, num_g_corr
        scrcoul_g_in(1,igp,iw,1) = zero
     enddo
  enddo
endif

!We store epsilon-1 to disk:
do iw = 1, nfs
   do ig = 1, num_g_corr
      scrcoul_g_in(ig,ig,iw,1) = scrcoul_g_in(ig,ig,iw,1) - one
   enddo
enddo

END SUBROUTINE invert_epsilon

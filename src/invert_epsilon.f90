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
SUBROUTINE invert_epsilon(scrcoul_g_in, lgamma)

USE freq_gw,       ONLY : nfs
USE gwsigma,       ONLY : gcutcorr
USE kinds,         ONLY : DP

IMPLICIT NONE    

COMPLEX(DP)       :: scrcoul_g_in(gcutcorr, gcutcorr, nfs, 1)
COMPLEX(DP)       :: work(gcutcorr)
INTEGER           :: ig, igp
INTEGER           :: iw
INTEGER           :: iwork(gcutcorr), info
LOGICAL           :: lgamma

!> complex constant of 0
COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

!> complex constant of 1
COMPLEX(dp), PARAMETER :: one = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

!at Gamma wings of \Chi are 0.
if(lgamma) then
  do iw = 1, nfs
    do ig = 2, gcutcorr
       scrcoul_g_in(ig,1,iw,1)  = zero 
    enddo
    do igp = 2, gcutcorr
       scrcoul_g_in(1,igp,iw,1) = zero
    enddo
  enddo
endif

!Need block inversion routine if iq is gamma.
do iw = 1, nfs
   call ZGETRF (gcutcorr, gcutcorr,&
   scrcoul_g_in(1:gcutcorr,1:gcutcorr,iw,1), gcutcorr, iwork, info)
   call errore ('invert epsilon', 'factorization', info)
   call ZGETRI (gcutcorr, scrcoul_g_in(1:gcutcorr,1:gcutcorr,iw,1),& 
   gcutcorr, iwork, work, gcutcorr, info)
   call errore ('invert epsilon', 'inversion', info)
enddo

write(6,*)
write(6,'(5x, "Done epsilon inversion.")') 
write(6,'(5x, "")') 

if(lgamma) then
  do iw = 1, nfs
     do ig = 2, gcutcorr
        scrcoul_g_in(ig,1,iw,1) = zero
     enddo
     do igp = 2, gcutcorr
        scrcoul_g_in(1,igp,iw,1) = zero
     enddo
  enddo
endif

!We store epsilon-1 to disk:
do iw = 1, nfs
   do ig = 1, gcutcorr
      scrcoul_g_in(ig,ig,iw,1) = scrcoul_g_in(ig,ig,iw,1) - one
   enddo
enddo

END SUBROUTINE invert_epsilon

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
SUBROUTINE para_pool(ngmunique, igstart, igstop)

USE mp_world ,    ONLY :  nproc, mpime
USE io_global,    ONLY : stdout
USE mp_pools,     ONLY : inter_pool_comm, npool, kunit, my_pool_id

IMPLICIT NONE

INTEGER :: ngmunique, ngpool, ngr, igs, igstart, igstop

#ifdef __PARA
      if (npool.gt.1) then
      ! number of g-vec per pool and reminder
        ngpool = ngmunique / npool
        ngr = ngmunique - ngpool * npool
      ! the remainder goes to the first ngr pools
        if ( my_pool_id < ngr ) ngpool = ngpool + 1
        igs = ngpool * my_pool_id + 1
        if ( my_pool_id >= ngr ) igs = igs + ngr
      ! the index of the first and the last g vec in this pool
        igstart = igs
        igstop = igs - 1 + ngpool
      else
#endif
       igstart = 1
       igstop = ngmunique
#ifdef __PARA
      endif
#endif
END SUBROUTINE para_pool

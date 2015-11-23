  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
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

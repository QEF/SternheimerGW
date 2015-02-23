SUBROUTINE para_img(ngmunique, igstart, igstop)

USE mp_world ,   ONLY :  nproc, mpime
USE io_global,   ONLY : stdout
USE mp_global,   ONLY : inter_image_comm, intra_image_comm, &
                        my_image_id, nimage, root_image

IMPLICIT NONE

INTEGER :: npool, ngmunique, ngpool, ngr, igs, igstart, igstop

#ifdef __PARA
      !npool = nproc
      npool = nimage
      if (npool.gt.1) then
      ! number of g-vec per pool and reminder
        ngpool = ngmunique / npool
        ngr = ngmunique - ngpool * npool
      ! the remainder goes to the first ngr pools
        if ( my_image_id < ngr ) ngpool = ngpool + 1
        igs = ngpool * my_image_id + 1
        if ( my_image_id >= ngr ) igs = igs + ngr
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
END SUBROUTINE

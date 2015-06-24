SUBROUTINE sigma_grids()
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : at, bg, tpiba2, tpiba
  USE fft_custom,       ONLY : fft_cus, set_custom_grid, ggent, gvec_init
  USE gvect,            ONLY : nl, ngm, g, nlm, gstart, gl, igtongl
  USE control_flags,    ONLY : gamma_only
  USE grid_subroutines, ONLY : realspace_grid_init
  !USE grid_subroutines, ONLY : realspace_grid_init_custom
  USE fft_base,         ONLY : dfftp
  USE klist,            ONLY : xk, nks
  USE gvect,            ONLY : gcutm
  USE stick_set,        ONLY : pstickset_custom
  USE mp,               ONLY : mp_sum, mp_max,mp_barrier
  USE mp_pools,         ONLY : inter_pool_comm, nproc_pool
  USE mp_bands,         ONLY : me_bgrp, nproc_bgrp, inter_bgrp_comm, &
                               intra_bgrp_comm, root_bgrp, ntask_groups
  USE fft_types,        ONLY : fft_dlay_descriptor, fft_dlay_allocate, &
                               fft_dlay_set, fft_dlay_scalar
  USE io_global,        ONLY :  stdout, ionode, ionode_id
  USE gwsigma,          ONLY : sigma_x_st, sigma_c_st
  USE gwsigma,          ONLY : ecutsex, ecutsco
  USE grid_subroutines, ONLY : realspace_grids_info 

  IMPLICIT NONE

  INTEGER    :: ng
  INTEGER :: ik, ngm_, ngs_, ngw_ , nogrp, kpoint
  REAL (DP) :: gkcut
  INTEGER :: me, nproc, inter_comm, intra_comm, root


!Define COMMUNICATORS here we're using the band group (bgrp) comm
!this is a sub division of the pools. 
  me = me_bgrp
  nproc = nproc_bgrp
  inter_comm = inter_bgrp_comm
  intra_comm = intra_bgrp_comm
  root = root_bgrp
  nogrp = ntask_groups

  IF (nks == 0) THEN
     !
     ! if k-points are automatically generated (which happens later)
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * MAX ( &
            &SQRT (SUM(bg (1:3, 1)**2) ), &
            &SQRT (SUM(bg (1:3, 2)**2) ), &
            &SQRT (SUM(bg (1:3, 3)**2) ) )
  ELSE
     gkcut = 0.0d0
     DO kpoint = 1, nks
        gkcut = MAX (gkcut, SQRT ( SUM(xk (1:3, kpoint)**2) ) )
     ENDDO
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! EXCHANGE GRID !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  sigma_x_st%ecutt  = ecutsex
  sigma_x_st%gcutmt  = ecutsex/tpiba2
  gkcut = (SQRT (sigma_x_st%ecutt) / tpiba + gkcut)**2
!Generate auxilliary exchange grid.
  do ng = 1, ngm
     if ( gl( igtongl (ng) ) .le. sigma_x_st%gcutmt ) sigma_x_st%ngmt = ng
     if ( gl( igtongl (ng) ) .le. sigma_x_st%gcutmt ) sigma_x_st%ngmt_g = ng
  enddo
  CALL set_custom_grid(sigma_x_st)
  !CALL realspace_grid_init_custom(sigma_x_st%dfftt, at, bg, sigma_x_st%gcutmt)
  CALL realspace_grid_init(sigma_x_st%dfftt, at, bg, sigma_x_st%gcutmt)
  CALL pstickset_custom( gamma_only, bg, sigma_x_st%gcutmt, gkcut, sigma_x_st%gcutmt, &
                  dfftp, sigma_x_st%dfftt, ngw_ , ngm_, ngs_, me, root, nproc, &
                  intra_comm, nogrp )
  CALL gvec_init(sigma_x_st, sigma_x_st%ngmt, intra_comm)
  sigma_x_st%initalized = .true.
  CALL ggent(sigma_x_st)
  WRITE(stdout, '(//5x,"Exchange Grid:")')
  WRITE(stdout, '(5x, "nr1, nr2, nr3")')
  WRITE(stdout, '(5x, 3i4)') sigma_x_st%nr1t, sigma_x_st%nr2t, sigma_x_st%nr3t
  WRITE(stdout, '(5x, "G-Vects Exx:")')
  WRITE(stdout, '(5x, f6.2, i5)') sigma_x_st%ecutt, sigma_x_st%ngmt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! CORRELATION GRID !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  sigma_c_st%ecutt   = ecutsco
  sigma_c_st%gcutmt  = ecutsco/tpiba2
  gkcut = (SQRT (sigma_c_st%ecutt) / tpiba + gkcut)**2
!Generate auxilliary correlation grid.
  do ng = 1, ngm
     if ( gl( igtongl (ng) ) .le. sigma_c_st%gcutmt ) sigma_c_st%ngmt = ng
     if ( gl( igtongl (ng) ) .le. sigma_c_st%gcutmt ) sigma_c_st%ngmt_g = ng
  enddo
  CALL set_custom_grid(sigma_c_st)
  CALL realspace_grid_init(sigma_c_st%dfftt, at, bg, sigma_c_st%gcutmt)
  CALL pstickset_custom( gamma_only, bg, sigma_c_st%gcutmt, gkcut, sigma_c_st%gcutmt, &
                  dfftp, sigma_c_st%dfftt, ngw_ , ngm_, ngs_, me, root, nproc, &
                  intra_comm, nogrp )
  CALL gvec_init(sigma_c_st, sigma_c_st%ngmt, intra_comm)
  sigma_c_st%initalized = .true.
  CALL ggent(sigma_c_st)
  WRITE(stdout, '(//5x,"Correlation Grid:")')
  WRITE(stdout, '(5x, "nr1, nr2, nr3")')
  WRITE(stdout, '(5x, 3i4)') sigma_c_st%nr1t, sigma_c_st%nr2t, sigma_c_st%nr3t
  WRITE(stdout, '(5x, "Ecut, Ngmsco:")')
  WRITE(stdout, '(5x, f6.2, i5, //)') sigma_c_st%ecutt, sigma_c_st%ngmt
  CALL realspace_grids_info(sigma_c_st%dfftt, sigma_c_st%dfftt, nproc_pool)
END SUBROUTINE

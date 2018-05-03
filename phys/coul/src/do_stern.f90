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
!> Solve the Sternheimer equation to obtain the screened Coulomb interaction.
!!
!! This routine is the driver routine for the solvers. Depending on the choice
!! in the input file either a direct or an iterative solver is used.
SUBROUTINE do_stern(config, grid, freq)

  USE constants,            ONLY: eps6
  USE control_gw,           ONLY: do_q0_only, solve_direct, do_epsil, plot_coul, &
                                  model_coul, tr2_gw
  USE disp,                 ONLY: nqs, num_k_pts, w_of_q_start, x_q, xk_kpoints
  USE freq_gw,              ONLY: nfs
  USE freqbins_module,      ONLY: freqbins_type
  USE gwsymm,               ONLY: ngmunique, ig_unique, use_symm, sym_friend, sym_ig
  USE io_global,            ONLY: stdout, meta_ionode
  USE kinds,                ONLY: dp
  USE klist,                ONLY: lgauss
  USE mp,                   ONLY: mp_sum, mp_barrier
  USE mp_images,            ONLY: inter_image_comm, my_image_id
  USE mp_world,             ONLY: mpime
  USE parallel_module,      ONLY: parallel_task, mp_gatherv
  USE plot_coulomb_module,  ONLY: plot_coulomb
  USE run_nscf_module,      ONLY: run_nscf
  USE select_solver_module, ONLY: select_solver_type
  USE sigma_grid_module,    ONLY: sigma_grid_type
  USE timing_module,        ONLY: time_coulomb, time_coul_nscf, time_coul_invert, &
                                  time_coul_io, time_coul_symm, time_coul_unfold, &
                                  time_coul_comm
  USE units_gw,             ONLY: lrcoul, iuncoul

IMPLICIT NONE

  !> stores the configuration of the linear solver for the screened Coulomb interaction
  TYPE(select_solver_type), INTENT(IN) :: config

  !> the FFT grid used to represent correlation and exchange
  TYPE(sigma_grid_type),    INTENT(IN) :: grid

  !> the definition of the frequency grid
  TYPE(freqbins_type),      INTENT(IN) :: freq

  !> the number of G vectors in the correlation grid
  INTEGER :: num_g_corr

  !> the number of tasks done on any process
  INTEGER, ALLOCATABLE :: num_task(:)

  !> the number of tasks done on this process
  INTEGER :: num_task_loc

  !> the part of the screened Coulomb interaction calculated on this process
  COMPLEX(dp), ALLOCATABLE :: scrcoul_loc(:,:,:)

  !> the screened Coulomb interaction gathered on the root process
  COMPLEX(dp), ALLOCATABLE :: scrcoul_root(:,:,:)

  !> the unfolded Coulomb interaction
  COMPLEX(dp), ALLOCATABLE :: scrcoul_g(:,:,:)

  !> the dielectric constant at the q + G = 0 point
  COMPLEX(dp), ALLOCATABLE :: eps_m(:)

  !> loop variables
  INTEGER :: iq, ig, igstart, igstop, ios, iq1, iq2
  LOGICAL :: do_band, do_iq, setup_pw, do_matel

  !> we need a special treatment for the case q + G = 0
  !! this flag indicates whether we are at this point
  LOGICAL :: lgamma

  !> is this process the root of the image
  LOGICAL :: is_root

  !> id of the image on which we collect the data
  INTEGER,     PARAMETER :: root_id = 0

  !> constant to initialize arrays to zero
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0, 0.0, KIND=dp)

  CALL start_clock(time_coulomb)

  ! set helper variable
  num_g_corr = grid%corr_fft%ngm

  ! some tasks are only done by the root process
  is_root = my_image_id == root_id

  IF (meta_ionode) THEN
    ALLOCATE(scrcoul_g(num_g_corr, num_g_corr, nfs))
  ELSE
    ALLOCATE(scrcoul_g(1,1,1))
  END IF

  ! allocate arrays for symmetry routine
  ALLOCATE(ig_unique(num_g_corr))
  ALLOCATE(sym_ig(num_g_corr))
  ALLOCATE(sym_friend(num_g_corr))
  ALLOCATE(eps_m(nfs))

  do_iq    = .TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  do_matel = .TRUE.

  IF(lgauss) WRITE(stdout, '(//5x,"SYSTEM IS METALLIC")')
  IF(.NOT.do_epsil) THEN
    iq1 = w_of_q_start
    iq2 = nqs
  ELSE
    ! In case we want to trace a line through the brillouin zone
    ! or get the screening for a particular grid q points (i.e. coulomb matel).
    iq1 = w_of_q_start
    iq2 = num_k_pts
  ENDIF
    
  DO iq = iq1, iq2
    ! Perform head of dielectric matrix calculation.
    CALL start_clock ('epsilq')
    do_matel = .FALSE.

    ! we must evaluate the q + G = 0 case differently, because we need to
    ! shift the vector by a small delta q so that the solver has a
    ! nonvanishing solution
    ! the gamma point is evaluated at the root process
    IF (.NOT.do_epsil) THEN
      lgamma = ALL(ABS(x_q(:,iq)) < eps6)
    ELSE
      lgamma = ALL(ABS(xk_kpoints(:,iq)) < eps6)
    END IF
    IF (lgamma .AND. is_root) THEN

      ! create and initialize array for dielectric constant at q + G = 0
      eps_m = zero

      CALL prepare_q0(do_band, do_iq, setup_pw, iq)

      ! nscf calculation to obtain the wave functions
      CALL start_clock(time_coul_nscf)
      CALL run_nscf(do_band, do_matel, iq)
      CALL stop_clock(time_coul_nscf)

      CALL initialize_gw(.TRUE.)
      CALL coulomb_q0G0(config, eps_m)
      WRITE(stdout,'(5x, "epsM(0) = ", f12.7)') eps_m(1)
      WRITE(stdout,'(5x, "epsM(iwp) = ", f12.7)') eps_m(2)
      CALL clean_pw_gw(.FALSE.)

    END IF ! gamma & root

    !
    ! now the general case for any q + G /= 0
    !
    CALL prepare_q(do_band, do_iq, setup_pw, iq)
 
    ! nscf calculation to obtain the wave functions
    CALL start_clock(time_coul_nscf)
    CALL run_nscf(do_band, do_matel, iq)
    CALL stop_clock(time_coul_nscf)

    CALL initialize_gw(.TRUE.)

    ! symmetrize the G vectors -> only unique ones are calculated
    IF (use_symm) THEN
      WRITE(stdout,'("")')
      WRITE(stdout,'(5x, "SYMMETRIZING COULOMB Perturbations")')
      WRITE(stdout,'("")')
      CALL start_clock(time_coul_symm)
      CALL stern_symm(num_g_corr)
      CALL stop_clock(time_coul_symm)
    ELSE
      ngmunique = num_g_corr
      DO ig = 1, num_g_corr
         ig_unique(ig) = ig
      END DO
    END IF ! use_symm

    ! distribute the number of tasks over processes
    CALL parallel_task(inter_image_comm, ngmunique, igstart, igstop, num_task)
    WRITE(stdout, '(5x, "iq ",i4, " igstart ", i4, " igstop ", i4)') iq, igstart, igstop

    ! note: my_image_id + 1 to convert from C to Fortran notation
    num_task_loc = num_task(my_image_id + 1)

    ! allocate array for distributed part of Coulomb on all processes
    ALLOCATE(scrcoul_loc(num_g_corr, nfs, num_task_loc))

    ! evaluate screened Coulomb interaction and collect on root
    CALL coulomb(config, igstart, num_g_corr, num_task_loc, scrcoul_loc)
    CALL start_clock(time_coul_comm)
    CALL mp_gatherv(inter_image_comm, root_id, num_task, scrcoul_loc, scrcoul_root)
    CALL stop_clock(time_coul_comm)

    ! Only the root of the image should write to file
    IF (meta_ionode) THEN

      ! unfold W from reduced array to full array
      ! also reorder the indices
      CALL start_clock(time_coul_unfold)
      CALL unfold_w(num_g_corr, scrcoul_root, scrcoul_g)
      CALL stop_clock(time_coul_unfold)

      ! set the special |q + G| = 0 element
      IF (lgamma) scrcoul_g(1,1,:) = eps_m

      ! for the direct solver W = eps^-1
      IF (solve_direct) THEN
        WRITE(1000+mpime, '("UNFOLDING, INVERTING, WRITING W")')
        CALL start_clock(time_coul_invert)
        CALL invert_epsilon(num_g_corr, scrcoul_g, lgamma)
        CALL stop_clock(time_coul_invert)
      END IF

      ! write to file
      CALL start_clock(time_coul_io)
      CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, +1, ios)
      CALL stop_clock(time_coul_io)

      IF (plot_coul) CALL plot_coulomb(model_coul, tr2_gw, grid, freq, scrcoul_g)

    END IF ! root

    DEALLOCATE(scrcoul_loc)
    IF (meta_ionode) DEALLOCATE(scrcoul_root)
    IF (ALLOCATED(num_task)) DEALLOCATE(num_task)

    CALL mp_barrier(inter_image_comm)
    CALL clean_pw_gw(.FALSE.)
    CALL print_clock ('epsilq')
    CALL stop_clock ('epsilq')
    IF(do_q0_only) EXIT

  END DO ! iq

  IF (meta_ionode) CLOSE(iuncoul)
  WRITE(stdout, '("Finished Calculating Screened Coulomb")')
  IF (ALLOCATED(scrcoul_g)) DEALLOCATE(scrcoul_g)
  DEALLOCATE(eps_m)
  DEALLOCATE(ig_unique)
  DEALLOCATE(sym_ig)
  DEALLOCATE(sym_friend)

  CALL stop_clock(time_coulomb)

end SUBROUTINE do_stern

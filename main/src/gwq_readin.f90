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
SUBROUTINE gwq_readin(config_coul, config_green, freq, vcut, debug)
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program GW.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE analytic_module,      ONLY : aaa_approx, aaa_pole, godby_needs, pade_approx, pade_robust
  USE cell_base,            ONLY : at, alat
  USE constants,            ONLY : RYTOEV, eps12
  USE control_flags,        ONLY : restart, lkpoint_dir, iverbosity, twfcollect
  USE control_gw,           ONLY : maxter, alpha_mix, reduce_io, tr2_gw, niter_gw, &
                                   lmax_gw, tr2_green, lmax_green, nmix_gw, ldisp, &
                                   tmp_dir_gw, tmp_dir_coul, eta, do_coulomb, do_sigma_c, &
                                   do_sigma_exx, do_sigma_matel, do_q0_only, maxter_green, &
                                   maxter_coul, model_coul, &
                                   solve_direct, do_epsil, set_alpha_pv, &
                                   do_imag, newgrid, double_grid, output_t => output, &
                                   plot_coul, method_truncation => truncation
  USE control_lr,           ONLY : lgamma, lrpa, alpha_pv
  USE debug_module,         ONLY : debug_type
  USE disp,                 ONLY : nq1, nq2, nq3, iq1, iq2, iq3, xk_kpoints, num_k_pts, & 
                                   w_of_q_start, w_of_k_start, w_of_k_stop
  USE freq_gw,              ONLY : fiu, nfs, wsigmamin, wsigmamax, nwsigma, wcoulmax, nwcoul, &
                                   wsig_wind_min, wsig_wind_max, nwsigwin
  USE freqbins_module,      ONLY : freqbins_type, no_symmetry
  USE gw_input_module,      ONLY : gw_input_type, gw_output_type, gw_input_read, gw_input_bcast
  USE gwsigma,              ONLY : nbnd_sig, ecutsex, ecutsco, corr_conv
  USE gwsymm,               ONLY : use_symm
  USE input_parameters,     ONLY : max_seconds, nk1, nk2, nk3, k1, k2, k3, force_symmorphic
  USE io_files,             ONLY : tmp_dir, prefix, check_tempdir
  USE io_global,            ONLY : meta_ionode, meta_ionode_id, stdout
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : nproc_pool_file, nproc_image_file
  USE mp_images,            ONLY : my_image_id, nproc_image
  USE mp_pools,             ONLY : nproc_pool
  USE mp_world,             ONLY : world_comm
  USE output_mod,           ONLY : filsigx, filsigc, filcoul
  USE qpoint,               ONLY : nksq
  USE run_info,             ONLY : title
  USE save_gw,              ONLY : tmp_dir_save
  USE select_solver_module, ONLY : select_solver_type
  USE start_k,              ONLY : reset_grid
  USE truncation_module
  USE wrappers,             ONLY : f_mkdir_safe
  !
  !
  IMPLICIT NONE
  !
  !> We store the configuration for the linear solver for W.
  TYPE(select_solver_type), INTENT(OUT) :: config_coul
  !
  !> We store the configuration for the linear solver for G.
  TYPE(select_solver_type), INTENT(OUT) :: config_green
  !
  !> We store the frequency for the Coulomb solver in this type.
  TYPE(freqbins_type), INTENT(OUT) :: freq
  !
  !> We store the truncated Coulomb potential in this type.
  TYPE(vcut_type),     INTENT(OUT) :: vcut
  !
  !> we store the debug information in this type
  TYPE(debug_type),    INTENT(OUT) :: debug
  !
  !> user input from file
  TYPE(gw_input_type) :: input
  !
  !> user specification for filenames
  TYPE(gw_output_type) :: output
  !
  !> size of the Wigner-Seitz cell
  REAL(dp) atws(3,3)
  !
  !> the cutoff used for the truncated potential
  REAL(dp) ecut_vcut
  !
  !> constant indicating unset priority
  INTEGER, PARAMETER :: no_solver = 0
  !
  !> the priority of the various solvers for the screened Coulomb interaction
  INTEGER priority_coul(10)
  !
  !> the priority of the various solvers for the Green's function
  INTEGER priority_green(10)
  !
  !> number of nontrivial priorities
  INTEGER num_priority
  !
  !> counter on the nontrivial priorities
  INTEGER ipriority
  !
  !> function to convert to lower case
  CHARACTER(LEN=1), EXTERNAL :: lowercase

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, iter, ierr
  ! integer variable for I/O control
  ! counter on polarizations
  ! counter on iterations
  !
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=1), EXTERNAL :: capital
  CHARACTER(LEN=6) :: int_to_char
  INTEGER                    :: i
  LOGICAL :: exst, parallelfs
  REAL(DP)           :: ar, ai
  !

  ! truncation method
  CHARACTER(LEN=trunc_length) :: truncation

  ! read the user input and store it in the input type
  IF (meta_ionode) THEN
    CALL gw_input_read(input, output)
  END IF
  ! broadcast the user input to all CPU
  CALL gw_input_bcast(input, output)
  title = input%title
  prefix = input%prefix
  tmp_dir = trimcheck(input%outdir)
  nbnd_sig = input%num_band
  do_imag = input%int_imag_axis
  nk1 = input%kpt_grid(1)
  nk2 = input%kpt_grid(2)
  nk3 = input%kpt_grid(3)
  nq1 = input%qpt_grid(1)
  nq2 = input%qpt_grid(2)
  nq3 = input%qpt_grid(3)
  do_coulomb = input%do_epsil .OR. input%do_coul
  do_epsil = input%do_epsil
  tr2_gw = input%thres_coul
  maxter_coul = input%max_iter_coul
  priority_coul = input%priority_coul
  lmax_gw = input%lmax_coul
  alpha_pv = input%shift_proj
  niter_gw = input%num_iter_coul
  nmix_gw = input%num_mix_coul
  use_symm = input%use_symm_coul
  solve_direct = (input%solve_coul == 'direct')
  IF (.NOT.solve_direct) THEN
    IF (input%solve_coul /= 'iter' .AND. input%solve_coul /= 'iterative') THEN
      CALL errore(__FILE__, 'unknown solver method for coulomb', 1)
    END IF
  END IF
  wcoulmax = input%max_freq_coul
  nwcoul = input%num_freq_coul
  plot_coul = input%plot_coul
  ! disregard case of input
  DO i = 1, LEN_TRIM(input%model_coul)
    input%model_coul(i:i) = lowercase(input%model_coul(i:i))
  END DO
  SELECT CASE (TRIM(input%model_coul))
  CASE ('gn', 'pp', 'godby-needs')
    model_coul = godby_needs
  CASE ('pade')
    model_coul = pade_approx
  CASE ('pade robust')
    model_coul = pade_robust
  CASE ('aaa')
    model_coul = aaa_approx
  CASE ('aaa pole')
    model_coul = aaa_pole
  CASE DEFAULT
    CALL errore(__FILE__, 'unknown screening model' // TRIM(input%model_coul), 1)
  END SELECT ! input%model_coul
  tr2_green = input%thres_green
  maxter_green = input%max_iter_green
  priority_green = input%priority_green
  lmax_green = input%lmax_green
  do_sigma_c = .NOT.input%do_epsil .AND. input%do_corr
  ecutsco = input%ecut_corr
  wsigmamin = input%min_freq_corr
  wsigmamax = input%max_freq_corr
  nwsigma = input%num_freq_corr
  do_sigma_exx = .NOT.input%do_epsil .AND. input%do_exch
  ecutsex = input%ecut_exch
  do_sigma_matel = .NOT.input%do_epsil .AND. input%do_matrix_el
  wsig_wind_min = input%min_freq_wind
  wsig_wind_max = input%max_freq_wind
  nwsigwin = input%num_freq_wind
  truncation = input%truncation
  do_q0_only = input%only_one_qpt
  w_of_q_start = input%first_qpt
  w_of_k_start = input%first_kpt
  w_of_k_stop = input%last_kpt
  IF (input%verbosity == 'low') THEN
    iverbosity = 0
  ELSE
    iverbosity = 1
  END IF
  debug = input%debug

  ! set Quantum ESPRESSO module variables
  lrpa        = .TRUE.
  ldisp       = .TRUE.
  double_grid = .FALSE.
  max_seconds =  1.E+7_DP
  restart     = .FALSE.

  ! set defaults currently for variables currently not in use

  !for bulk systems alpha_mix = 0.7 is standard
  !for slab systems more rapid convergence can
  !be obtained with alpha_mix = 0.3.
  alpha_mix(:) = 0.D0
  alpha_mix(1) = 0.7D0
  reduce_io    = .FALSE.
  k1           = 0
  k2           = 0
  k3           = 0
  iq1          = 0
  iq2          = 0
  iq3          = 0

  ! interpret the truncation scheme
  SELECT CASE (truncation)
  CASE (NO_TRUNCATION_1, NO_TRUNCATION_2, NO_TRUNCATION_3, NO_TRUNCATION_4, NO_TRUNCATION_5)
    method_truncation = NO_TRUNCATION

  CASE (SPHERICAL_TRUNCATION_1, SPHERICAL_TRUNCATION_2, SPHERICAL_TRUNCATION_3, &
        SPHERICAL_TRUNCATION_4, SPHERICAL_TRUNCATION_5)
    method_truncation = SPHERICAL_TRUNCATION

  CASE (FILM_TRUNCATION_1, FILM_TRUNCATION_2, FILM_TRUNCATION_3, FILM_TRUNCATION_4)
    method_truncation = FILM_TRUNCATION

  CASE (VCUT_SPHERICAL_TRUNCATION_1, VCUT_SPHERICAL_TRUNCATION_2, &
        VCUT_SPHERICAL_TRUNCATION_3, VCUT_SPHERICAL_TRUNCATION_4)
    method_truncation = VCUT_SPHERICAL_TRUNCATION

  CASE (VCUT_WIGNER_SEITZ_TRUNCATION_1, VCUT_WIGNER_SEITZ_TRUNCATION_2, &
        VCUT_WIGNER_SEITZ_TRUNCATION_3, VCUT_WIGNER_SEITZ_TRUNCATION_4)
    method_truncation = VCUT_WIGNER_SEITZ_TRUNCATION

  END SELECT ! truncation

  ! convert frequencies to Ry
  wsigmamin = wsigmamin / RYTOEV
  wsigmamax = wsigmamax / RYTOEV
  wcoulmax  = wcoulmax  / RYTOEV
  wsig_wind_max = wsig_wind_max / RYTOEV
  wsig_wind_min = wsig_wind_min / RYTOEV

  ! copy read data to output type
  filsigx                         = output%file_exch
  filsigc                         = output%file_corr
  filcoul                         = output%file_coul
  output_t%directory              = output%directory
  output_t%pp_dft%filename        = output%file_dft 
  output_t%pp_gw%filename         = output%file_gw
  output_t%pp_vxc%filename        = output%file_vxc
  output_t%pp_exchange%filename   = output%file_hf
  output_t%pp_renorm%filename     = output%file_renorm
  output_t%pp_re_corr%filename    = ''
  output_t%pp_re_corr_iw%filename = ''
  output_t%pp_im_corr%filename    = ''
  output_t%pp_im_corr_iw%filename = ''
  output_t%pp_spec%filename       = ''
  output_t%pp_spec_iw%filename    = ''
  output_t%file_sigma             = output%file_sigma

! if corr_conv not set in input file default to the full
! correlation cutoff.
  if(ABS(corr_conv) < eps12) corr_conv = ecutsco

  !
  ! ... Check all namelist variables
  !
  IF (tr2_gw <= 0.D0) CALL errore (' gwq_readin', ' Wrong tr2_gw ', 1)
  IF (tr2_green <= 0.D0) CALL errore (' gwq_readin', ' Wrong tr2_green ', 1)

  DO iter = 1, maxter
     IF (alpha_mix (iter) .LT.0.D0.OR.alpha_mix (iter) .GT.1.D0) CALL &
          errore ('gwq_readin', ' Wrong alpha_mix ', iter)
  ENDDO
  IF (niter_gw.LT.1.OR.niter_gw.GT.maxter) CALL errore ('gwq_readin', &
       ' Wrong niter_gw ', 1)
  IF (nmix_gw.LT.1.OR.nmix_gw.GT.5) CALL errore ('gwq_readin', ' Wrong &
       &nmix_gw ', 1)
  !
  IF (iverbosity.NE.0.AND.iverbosity.NE.1) CALL errore ('gwq_readin', &
       &' Wrong  iverbosity ', 1)
  IF (max_seconds.LT.0.1D0) CALL errore ('gwq_readin', ' Wrong max_seconds', 1)

! HL here we can just use this to readin the list of frequencies that we want to calculate
! Stored in array  fiu(:), of size nfs.
! reads the frequencies ( just if fpol = .true. )
  nfs=0
  IF (meta_ionode) THEN
     READ (5, *, iostat = ios) card
     IF ( TRIM(card)=='FREQUENCIES'.OR. &
          TRIM(card)=='frequencies'.OR. &
          TRIM(card)=='Frequencies') THEN
        READ (5, *, iostat = ios) nfs
     ENDIF
  ENDIF

  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore ('gwq_readin', 'reading number of FREQUENCIES', ABS(ios) )
  CALL mp_bcast(nfs, meta_ionode_id, world_comm )

  if (nfs < 1) call errore('gwq_readin','Too few frequencies',1)
  ALLOCATE(fiu(nfs), freq%solver(nfs))

  IF (meta_ionode) THEN
     IF ( TRIM(card) == 'FREQUENCIES' .OR. &
          TRIM(card) == 'frequencies' .OR. &
          TRIM(card) == 'Frequencies' ) THEN
        DO i = 1, nfs
           !HL Need to convert frequencies from electron volts into Rydbergs
           READ (5, *, iostat = ios) ar, ai 
           freq%solver(i) = CMPLX(ar, ai, KIND=dp) / RYTOEV
        END DO
     END IF
  END IF

  ! set the small shift into the complex plane
  IF (do_imag) THEN
    !
    ! if we are already in the complex plane, we don't need to shift
    freq%eta = 0.0_dp
    eta = input%eta / RYTOEV
    !
  ELSE
    !
    ! if we are on the real axis, we shift by a small amount into the
    ! complex plane for a numerically stable treatment of the poles
    freq%eta = input%eta / RYTOEV
    eta = input%eta / RYTOEV
    !
  END IF

  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore ('gwq_readin', 'reading FREQUENCIES card', ABS(ios) )
  CALL mp_bcast(freq%solver, meta_ionode_id, world_comm )
  fiu = freq%solver

  ! use symmetry for the frequencies (only for Pade or AAA approximation)
  IF (model_coul == pade_approx .OR. model_coul == aaa_approx .OR. model_coul == aaa_pole) THEN
    freq%freq_symm_coul = input%freq_symm_coul
  ELSE
    ! symmetry not implemented for robust Pade and Godby-Needs
    freq%freq_symm_coul = no_symmetry
  END IF

  num_k_pts = 0
  IF (meta_ionode) THEN
     READ(5, *, IOSTAT = ios) card
     READ(5, *, IOSTAT = ios) card
     IF (TRIM(card) == 'K_POINTS' .OR. &
         TRIM(card) == 'k_points' .OR. &
         TRIM(card) == 'K_points') THEN
       READ(5, *, IOSTAT = ios) num_k_pts
     END IF
  END IF
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore(__FILE__, 'reading number of kpoints', ABS(ios))
  CALL mp_bcast(num_k_pts, meta_ionode_id, world_comm)
  IF (num_k_pts > 2000) CALL errore(__FILE__,'Too many k-points', 1)
  IF (num_k_pts < 1) CALL errore(__FILE__,'Too few kpoints', 1)
  IF (meta_ionode) THEN
    IF (TRIM(card) == 'K_POINTS' .OR. &
        TRIM(card) == 'k_points' .OR. &
        TRIM(card) == 'K_points') THEN
      DO i = 1, num_k_pts
        !should be in units of 2pi/a0 cartesian co-ordinates
        READ(5, *, IOSTAT = ios) xk_kpoints(:,i)
      END DO
    END IF
  END IF
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore(__FILE__, 'reading KPOINTS card', ABS(ios))
  CALL mp_bcast(xk_kpoints, meta_ionode_id, world_comm)

 IF (.NOT.do_epsil .AND. w_of_k_stop == -2) THEN
    w_of_k_stop = num_k_pts
 END IF

 ! sanity check on k-point window
 IF (w_of_k_stop < w_of_k_start) THEN
   IF (.NOT.do_epsil) CALL errore (__FILE__, 'w_of_k_stop less than w_of_k_start', 1)
 END IF
 IF (w_of_k_start < 1 .OR. w_of_k_stop > num_k_pts) THEN
   CALL errore(__FILE__, 'w_of_k_start or w_of_k_stop out of bounds [1:number k-points]', 1)
 END IF
   

  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  tmp_dir_save=tmp_dir
  tmp_dir_gw= TRIM (tmp_dir) //'_gw'//trim(int_to_char(my_image_id))//'/'
  tmp_dir_coul= TRIM (tmp_dir) //'_gw0'//'/'

  ! set output directory if not defined
  IF (output_t%directory == '') THEN
    output_t%directory = trimcheck(tmp_dir_gw)
  ELSE
    output_t%directory = trimcheck(output%directory)
  END IF
  output_t%prefix = prefix

  ! create directory (if it doesn't exist)
  ierr = f_mkdir_safe(output_t%directory)
  IF (ierr > 0) CALL errore(__FILE__, "error when opening/creating directory for output", ierr)

  ! augment sigma file with output directory
  output_t%file_sigma = TRIM(output_t%directory) // output_t%file_sigma

  CALL check_tempdir ( tmp_dir_gw, exst, parallelfs )

  CALL read_file ( )
  force_symmorphic = .true.
  CALL mp_bcast(force_symmorphic, meta_ionode_id, world_comm )
  IF(.not.force_symmorphic) then
      CALL errore( 'FORCE_SYMMORPHIC must be true in GROUND STATE CALCULATIONS!', 'gwq_readin', 1)
  ENDIF

  newgrid = reset_grid (nk1, nk2, nk3, k1, k2, k3)
  tmp_dir=tmp_dir_save

  IF (nproc_image /= nproc_image_file .and. .not. twfcollect)  &
     CALL errore('gwq_readin',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('gwq_readin',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)
  !
  ! If a band structure calculation needs to be done do not open a file 
  ! for k point
  !
  lkpoint_dir=.FALSE.
  !
  IF (.NOT.ldisp) THEN
     IF (lgamma) THEN
        nksq = nks
     ELSE
        nksq = nks / 2
     ENDIF
  ENDIF
  IF (.NOT.do_epsil) THEN
    IF (nq1 <= 0) CALL errore(__FILE__,'1st component of qpt_grid must be positive', 1)
    IF (nq2 <= 0) CALL errore(__FILE__,'2nd component of qpt_grid must be positive', 1)
    IF (nq3 <= 0) CALL errore(__FILE__,'3rd component of qpt_grid must be positive', 1)
  END IF

  ! if alpha_pv was not set in the input, we determine it automatically
  set_alpha_pv = (alpha_pv < 0)

  !
  ! setup the truncation
  !
  ! note: this step is computationally expensive, so we only do it if necessary
  IF (method_truncation == VCUT_SPHERICAL_TRUNCATION .OR. &
      method_truncation == VCUT_WIGNER_SEITZ_TRUNCATION) THEN
    !
    ! determine supercell
    atws = alat * at
    !
    atws(:,1) = atws(:,1) * nq1
    atws(:,2) = atws(:,2) * nq2
    atws(:,3) = atws(:,3) * nq3
    !
    ! we should use a quarter of the cutoff, because vcut assumes WF cutoff
    ! and converts to density cutoff, but for some reason the scaling is
    ! a bit different then for the custom FFT type, so that we increase the
    ! prefactor to 0.3 to be on the safe side
    ecut_vcut = 0.30_dp * MAX(ecutsco, ecutsex)
    CALL vcut_reinit(vcut, atws, ecut_vcut, tmp_dir_gw)
    CALL vcut_info(stdout, vcut)
    !
  END IF ! vcut truncation methods

  !
  ! setup the linear solver
  !
  config_coul%max_iter = maxter_coul
  config_coul%threshold = tr2_gw
  config_coul%bicg_lmax = lmax_gw
  config_green%max_iter = maxter_green
  config_green%threshold = tr2_green
  config_green%bicg_lmax = lmax_green
  !
  ! setup priority for Coulomb solver
  num_priority = COUNT(priority_coul /= no_solver)
  !
  IF (num_priority == 0) THEN
    CALL errore(__FILE__, "priority for Coulomb solver not specified", 1)
  END IF
  !
  ALLOCATE(config_coul%priority(num_priority))
  ipriority = 0
  !
  DO i = 1, SIZE(priority_coul)
    !
    IF (priority_coul(i) /= no_solver) THEN
      ipriority = ipriority + 1
      config_coul%priority(ipriority) = priority_coul(i)
    END IF
    !
  END DO ! i
  !
  ! setup priority for Green's solver
  num_priority = COUNT(priority_green /= no_solver)
  !
  IF (num_priority == 0) THEN
    CALL errore(__FILE__, "priority for Green solver not specified", 1)
  END IF
  !
  ALLOCATE(config_green%priority(num_priority))
  ipriority = 0
  !
  DO i = 1, SIZE(priority_green)
    !
    IF (priority_green(i) /= no_solver) THEN
      ipriority = ipriority + 1
      config_green%priority(ipriority) = priority_green(i)
    END IF
    !
  END DO ! i

  FLUSH(stdout)

END SUBROUTINE gwq_readin

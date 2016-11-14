!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
SUBROUTINE gwq_readin(config_green, freq, vcut, debug)
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program GW.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE cell_base,            ONLY : at, alat
  USE constants,            ONLY : RYTOEV, eps12
  USE control_flags,        ONLY : restart, lkpoint_dir, iverbosity, modenum, twfcollect
  USE control_gw,           ONLY : maxter, alpha_mix, lgamma, lgamma_gamma, epsil, &
                                   reduce_io, tr2_gw, niter_gw, lmax_gw, tr2_green, lmax_green, &
                                   nmix_gw, ldisp, recover, lrpa, lnoloc, start_irr, &
                                   last_irr, start_q, last_q, tmp_dir_gw, tmp_dir_coul, &
                                   ext_recover, ext_restart, modielec, eta, &
                                   do_coulomb, do_sigma_c, do_sigma_exx, do_green, do_sigma_matel, &
                                   do_q0_only, maxter_green, maxter_coul, godbyneeds, padecont,&
                                   cohsex, multishift, do_sigma_extra, &
                                   solve_direct, w_green_start, tinvert, coul_multishift,&
                                   trunc_2d, do_epsil, &
                                   do_diag_g, do_diag_w, do_imag, do_pade_coul, newgrid,&
                                   high_io, prec_direct, prec_shift, just_corr,&
                                   double_grid, name_length, output, &
                                   method_truncation => truncation
  USE debug_module,         ONLY : debug_type
  USE disp,                 ONLY : nq1, nq2, nq3, iq1, iq2, iq3, &
                                   xk_kpoints, kpoints, num_k_pts, & 
                                   w_of_q_start, w_of_k_start, w_of_k_stop
  USE freq_gw,              ONLY : fiu, nfs, wsigmamin, wsigmamax, nwsigma, wcoulmax, nwcoul, &
                                   wsig_wind_min, wsig_wind_max, nwsigwin
  USE freqbins_module,      ONLY : freqbins_type
  USE gwsigma,              ONLY : nbnd_sig, ecutsex, ecutsco, ecutprec, corr_conv, exch_conv
  USE gwsymm,               ONLY : use_symm
  USE input_parameters,     ONLY : max_seconds, nk1, nk2, nk3, k1, k2, k3, force_symmorphic
  USE io_files,             ONLY : tmp_dir, prefix
  USE io_global,            ONLY : meta_ionode, meta_ionode_id, stdout
  USE ions_base,            ONLY : nat, amass
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, nks, nkstot
  USE lsda_mod,             ONLY : nspin
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : nproc_pool_file, nproc_image_file
  USE mp_images,            ONLY : my_image_id, nproc_image
  USE mp_pools,             ONLY : nproc_pool
  USE mp_world,             ONLY : world_comm
  USE output_mod,           ONLY : fildyn, fildvscf, fildrho, filsigx, filsigc, filcoul
  USE parameters,           ONLY : nsx
  USE partial,              ONLY : atomo, list, nat_todo, nrapp
  USE qpoint,               ONLY : nksq, xq
  USE run_info,             ONLY : title
  USE save_gw,              ONLY : tmp_dir_save
  USE select_solver_module, ONLY : select_solver_type, bicgstab_multi, bicgstab_no_multi, &
                                   sgw_linear_solver
  USE start_k,              ONLY : reset_grid
  USE truncation_module
  USE wrappers,             ONLY : f_mkdir_safe
  !
  !
  IMPLICIT NONE
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
  !> size of the Wigner-Seitz cell
  REAL(dp) atws(3,3)
  !
  !> the cutoff used for the truncated potential
  REAL(dp) ecut_vcut
  !
  !> constant indicating unset priority
  INTEGER, PARAMETER :: no_solver = 0
  !
  !> the priority of the various solvers
  INTEGER priority_green(10)
  !
  !> number of nontrivial priorities
  INTEGER num_priority
  !
  !> counter on the nontrivial priorities
  INTEGER ipriority

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ipol, iter, na, ierr
  ! integer variable for I/O control
  ! counter on polarizations
  ! counter on iterations
  ! counter on atoms
  ! counter on types
  REAL(DP) :: amass_input(nsx)
  ! save masses read from input here
  !
  CHARACTER(LEN=256)         :: outdir
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=1), EXTERNAL :: capital
  CHARACTER(LEN=6) :: int_to_char
  INTEGER                    :: i
  LOGICAL                    :: nogg
  INTEGER, EXTERNAL  :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  LOGICAL, EXTERNAL  :: imatches
  LOGICAL :: exst, parallelfs
  REAL(DP)           :: ar, ai
  !
  ! output configuration
  CHARACTER(LEN=name_length) directory
  CHARACTER(LEN=name_length) file_dft
  CHARACTER(LEN=name_length) file_gw
  CHARACTER(LEN=name_length) file_vxc
  CHARACTER(LEN=name_length) file_exchange
  CHARACTER(LEN=name_length) file_renorm
  CHARACTER(LEN=name_length) file_re_corr
  CHARACTER(LEN=name_length) file_re_corr_iw
  CHARACTER(LEN=name_length) file_im_corr
  CHARACTER(LEN=name_length) file_im_corr_iw
  CHARACTER(LEN=name_length) file_spec
  CHARACTER(LEN=name_length) file_spec_iw
  CHARACTER(LEN=name_length) file_sigma

  ! truncation method
  CHARACTER(LEN=trunc_length) :: truncation

  NAMELIST / INPUTGW / tr2_gw, lmax_gw, amass, alpha_mix, niter_gw, nmix_gw,  &
                       nat_todo, iverbosity, outdir, epsil,  &
                       nrapp, max_seconds, reduce_io, &
                       modenum, prefix, fildyn, fildvscf, fildrho,   &
                       ldisp, nq1, nq2, nq3, iq1, iq2, iq3,   &
                       recover, lrpa, lnoloc, start_irr, last_irr, &
                       start_q, last_q, nogg, modielec, nbnd_sig, eta, kpoints,&
                       ecutsco, ecutsex, corr_conv, exch_conv, ecutprec, do_coulomb, do_sigma_c, do_sigma_exx, do_green,& 
                       do_sigma_matel, tr2_green, lmax_green, do_q0_only, wsigmamin, &
                       wsigmamax, wcoulmax, nwsigma, priority_green, &
                       use_symm, maxter_green, maxter_coul, w_of_q_start, w_of_k_start, w_of_k_stop, godbyneeds,& 
                       padecont, cohsex, multishift, do_sigma_extra,&
                       solve_direct, w_green_start, tinvert, coul_multishift, trunc_2d,&
                       do_epsil, do_diag_g, do_diag_w, do_imag, do_pade_coul, nk1, nk2, nk3, high_io,&
                       prec_direct, tmp_dir, prec_shift, just_corr,& 
                       nwcoul, double_grid, wsig_wind_min, wsig_wind_max, nwsigwin, truncation, &
                       filsigx, filsigc, filcoul, debug
  NAMELIST / OUTPUTGW / file_dft, file_gw, file_vxc, file_exchange, file_renorm, &
                       file_re_corr, file_re_corr_iw, file_im_corr, file_im_corr_iw, &
                       file_spec, file_spec_iw, directory, file_sigma

  ! alpha_mix    : the mixing parameter
  ! niter_gw     : maximum number of iterations
  ! nmix_gw      : number of previous iterations used in mixing
  ! nat_todo     : number of atom to be displaced
  ! iverbosity   : verbosity control
  ! outdir       : directory where input, output, temporary files reside
  ! max_seconds  : maximum cputime for this run
  ! reduce_io    : reduce I/O to the strict minimum
  ! prefix       : the prefix of files produced by pwscf
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! filsigx      : output file containing exchange part of sigma
  ! filsigc      : output file containing correlation part of sigma
  ! filcoul      : output file containing screened coulomb
  ! eth_rps      : threshold for calculation of  Pc R |psi> (Raman)
  ! eth_ns       : threshold for non-scf wavefunction calculation (Raman)
  ! recover      : recover=.true. to restart from an interrupted run
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     : 
  ! nogg         : if .true. lgamma_gamma tricks are not used

  IF (meta_ionode) THEN

  !
  ! flib/inpfile.f90!
  ! Reads in from standar input (5)
  ! 
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
     READ( 5, '(A)', IOSTAT = ios ) title
  !   WRITE(6,*) nsx, maxter
  !
  ENDIF
  !
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore( 'gwq_readin', 'reading title ', ABS( ios ) )
  CALL mp_bcast(title, meta_ionode_id, world_comm  )
  !
  ! Rewind the input if the title is actually the beginning of inputgw namelist
  IF( imatches("&inputgw", title)) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    IF (meta_ionode) REWIND(5, iostat=ios)
    CALL mp_bcast(ios, meta_ionode_id, world_comm  )
    CALL errore('gwq_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  ! ... set default values for variables in namelist
  !
  tr2_gw       = 1.D-4
  tr2_green    = 1.D-3
  lmax_gw      = 4
  lmax_green   = 4
  priority_green(1) = bicgstab_multi
  priority_green(2) = sgw_linear_solver
  priority_green(3:) = no_solver
  amass(:)     = 0.D0
  alpha_mix(:) = 0.D0
  !for bulk systems alpha_mix = 0.7 is standard
  !for slab systems more rapid convergence can
  !be obtained with alpha_mix = 0.3.
  alpha_mix(1) = 0.7D0
  niter_gw     = maxter
  nmix_gw      = 3
  nat_todo     = 0
  modenum      = 0
  nrapp        = 0
  iverbosity   = 0
  lnoloc       = .FALSE.
  epsil        = .FALSE.
  max_seconds  =  1.E+7_DP
  reduce_io    = .FALSE.
  prec_direct  = .FALSE.
  prec_shift  = .FALSE.
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix       = 'pwscf'
  fildyn       = 'matdyn'
  fildrho      = ' '
  fildvscf     = ' '
  filsigx      = 'sigma_x'
  filsigc      = 'sigma_c'
  filcoul      = 'coulomb'
  nq1          = 0
  nq2          = 0
  nq3          = 0
  iq1          = 0
  iq2          = 0
  iq3          = 0
  nogg         = .FALSE.
  recover      = .FALSE.
  start_irr    = 0
  last_irr     =-1000
  start_q      = 1
  last_q       =-1000
  ldisp        =.FALSE.
  lrpa         =.FALSE.
  maxter_coul  = 160
  maxter_green = 220
  w_green_start = 1

  coul_multishift = .FALSE.
  trunc_2d        = .FALSE.
  do_epsil        = .FALSE.
  do_diag_g       = .FALSE.
  do_diag_w       = .FALSE.
  do_imag         = .FALSE.
  do_pade_coul    = .FALSE.
  double_grid     = .TRUE.
  high_io    = .TRUE.
!Sigma cutoff, correlation cutoff, exchange cutoff
!this is in case we want to define different cutoffs for 
!W and G. G cannot exceed sigma.
  ecutsco      = 5.0
  ecutsex      = 5.0
  ecutprec     = 15.0
  nbnd_sig     = 8
  nwcoul       = 35
  nwsigma      = 11
  nwsigwin     = 801
!Should have a catch if no model for screening is chosen...
  modielec     = .FALSE.
  godbyneeds   = .FALSE.
  cohsex       = .FALSE.
  padecont     = .FALSE.
  multishift   = .FALSE.
!Imaginary component added to linear system should be in Rydberg
  eta            =  0.02
  kpoints        = .FALSE.
  do_coulomb     = .FALSE.
  do_sigma_c     = .FALSE.
  do_sigma_exx   = .FALSE.
  do_green       = .FALSE.
  do_sigma_matel = .FALSE.
  do_sigma_extra = .FALSE.
  do_q0_only     = .FALSE.
  solve_direct   = .FALSE.
  tinvert        = .TRUE.
!Frequency variables
  wsigmamin      = 0.0d0
  wsigmamax      = 20.0d0
  wcoulmax       = 80.0d0   
  wsig_wind_min   = -50.0
  wsig_wind_max   =  30.0

!Symmetry Default:yes!, which q, point to start on.
!can be used in conjunction with do_q0_only.
  use_symm       = .TRUE.
  w_of_q_start   = 1
  w_of_k_start   = 1
  w_of_k_stop    = -2
  w_green_start  = 1 
! ...  reading the namelist inputgw
  just_corr = .FALSE.

  ! use the default truncation scheme
  truncation = 'on'

  IF (meta_ionode) THEN
    READ( 5, INPUTGW, ERR=30, IOSTAT = ios )

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

  END IF ! meta_ionode

  ! set defaults for output
  directory       = ''
  file_dft        = ''
  file_gw         = ''
  file_vxc        = ''
  file_exchange   = ''
  file_renorm     = ''
  file_re_corr    = ''
  file_re_corr_iw = ''
  file_im_corr    = ''
  file_im_corr_iw = ''
  file_spec       = ''
  file_spec_iw    = ''
  file_sigma      = 'sigma.xml'

  ! read the output from file
  IF (meta_ionode) READ(5, OUTPUTGW, ERR=30, IOSTAT = ios)
  
  ! copy read data to output type
  output%directory              = directory
  output%pp_dft%filename        = file_dft 
  output%pp_gw%filename         = file_gw
  output%pp_vxc%filename        = file_vxc
  output%pp_exchange%filename   = file_exchange
  output%pp_renorm%filename     = file_renorm
  output%pp_re_corr%filename    = file_re_corr
  output%pp_re_corr_iw%filename = file_re_corr_iw
  output%pp_im_corr%filename    = file_im_corr
  output%pp_im_corr_iw%filename = file_im_corr_iw
  output%pp_spec%filename       = file_spec
  output%pp_spec_iw%filename    = file_spec_iw
  output%file_sigma             = file_sigma

! if corr_conv not set in input file default to the full
! correlation cutoff.
  if(ABS(corr_conv) < eps12) corr_conv = ecutsco
!HL TEST PARA FINE
30 CALL mp_bcast(ios, meta_ionode_id, world_comm )
   CALL errore( 'gwq_readin', 'reading namelist', ABS( ios ) )
  IF (meta_ionode) tmp_dir = trimcheck (outdir)

  CALL bcast_gw_input ( ) 
  CALL mp_bcast(nogg, meta_ionode_id, world_comm  )

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
  IF (fildyn.EQ.' ') CALL errore ('gwq_readin', ' Wrong fildyn ', 1)
  IF (max_seconds.LT.0.1D0) CALL errore ('gwq_readin', ' Wrong max_seconds', 1)

  IF (nat_todo.NE.0.AND.nrapp.NE.0) CALL errore ('gwq_readin', &
       &' incompatible flags', 1)
  IF (modenum < 0) CALL errore ('gwq_readin', ' Wrong modenum ', 1)
  !
  !
  IF (meta_ionode) THEN
    ios = 0 
     IF (.NOT. ldisp) &
        READ (5, *, iostat = ios) (xq (ipol), ipol = 1, 3)
  END IF

  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore ('gwq_readin', 'reading xq', ABS (ios) )
  CALL mp_bcast(xq, meta_ionode_id, world_comm  )


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
  fiu = freq%solver

  ! set the small shift into the complex plane
  IF (do_imag) THEN
    !
    ! if we are already in the complex plane, we don't need to shift
    freq%eta = 0.0_dp
    !
  ELSE
    !
    ! if we are on the real axis, we shift by a small amount into the
    ! complex plane for a numerically stable treatment of the poles
    freq%eta = eta
    !
  END IF

  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore ('gwq_readin', 'reading FREQUENCIES card', ABS(ios) )
  CALL mp_bcast(fiu, meta_ionode_id, world_comm )

 IF (kpoints) then
     num_k_pts = 0
     IF (meta_ionode) THEN
        READ (5, *, iostat = ios) card
        READ (5, *, iostat = ios) card
        IF ( TRIM(card)=='K_POINTS'.OR. &
             TRIM(card)=='k_points'.OR. &
             TRIM(card)=='K_points') THEN
           READ (5, *, iostat = ios) num_k_pts
        ENDIF
     ENDIF
     CALL mp_bcast(ios, meta_ionode_id, world_comm )
     CALL errore ('pwq_readin', 'reading number of kpoints', ABS(ios) )
     CALL mp_bcast(num_k_pts, meta_ionode_id, world_comm )
     if (num_k_pts > 2000) call errore('phq_readin','Too many k-points',1) 
     if (num_k_pts < 1) call errore('phq_readin','Too few kpoints',1) 
     IF (meta_ionode) THEN
        IF ( TRIM(card)=='K_POINTS'.OR. &
             TRIM(card)=='k_points'.OR. &
             TRIM(card)=='K_points') THEN
           DO i = 1, num_k_pts
              !should be in units of 2pi/a0 cartesian co-ordinates
              READ (5, *, iostat = ios) xk_kpoints(1,i),& 
                         xk_kpoints(2,i), xk_kpoints(3,i)
           END DO
        END IF
     END IF
     CALL mp_bcast(ios, meta_ionode_id, world_comm)
     CALL errore ('gwq_readin', 'reading KPOINTS card', ABS(ios) )
     CALL mp_bcast(xk_kpoints, meta_ionode_id, world_comm)
 ELSE
     num_k_pts = 1
 ENDIF

 if (w_of_k_stop==-2) then
    w_of_k_stop = num_k_pts
 endif

 if (w_of_k_stop.lt.w_of_k_start) then
     CALL errore ('gwq_readin', 'w_of_k_stop less than w_of_k_start', ABS(ios) )
 else if ((w_of_k_stop.lt.1) .or. (w_of_k_start.lt.1)) then
     CALL errore ('gwq_readin', 'w_of_k_stop or w_of_k_start cannot be less than 1', ABS(ios) )
 endif
   

  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  amass_input(:)= amass(:)
  !
  tmp_dir_save=tmp_dir
  tmp_dir_gw= TRIM (tmp_dir) //'_gw'//trim(int_to_char(my_image_id))//'/'
  tmp_dir_coul= TRIM (tmp_dir) //'_gw0'//'/'

  ! set output directory if not defined
  IF (output%directory == '') THEN
    output%directory = trimcheck(tmp_dir_gw)
  ELSE
    output%directory = trimcheck(output%directory)
  END IF
  output%prefix = prefix

  ! create directory (if it doesn't exist)
  ierr = f_mkdir_safe(output%directory)
  IF (ierr > 0) CALL errore(__FILE__, "error when opening/creating directory for output", ierr)

  ! augment sigma file with output directory
  output%file_sigma = TRIM(output%directory) // output%file_sigma

  CALL check_tempdir ( tmp_dir_gw, exst, parallelfs )
  ext_restart=.FALSE.
  ext_recover=.FALSE.
  recover=.false.

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
  restart = recover
  !
  lgamma_gamma=.FALSE.
  IF (.NOT.ldisp) THEN
     IF (nkstot==1.OR.(nkstot==2.AND.nspin==2)) THEN
        lgamma_gamma=(lgamma.AND.(ABS(xk(1,1))<1.D-12) &
                            .AND.(ABS(xk(2,1))<1.D-12) &
                            .AND.(ABS(xk(3,1))<1.D-12) )
     ENDIF
     IF (nogg) lgamma_gamma=.FALSE.
     IF ((nat_todo /= 0 .or. nrapp /= 0 ) .and. lgamma_gamma) CALL errore( &
        'gwq_readin', 'gamma_gamma tricks with nat_todo or nrapp &
       & not available. Use nogg=.true.', 1)
     !
     IF (lgamma) THEN
        nksq = nks
     ELSE
        nksq = nks / 2
     ENDIF
  ENDIF
  IF (nat_todo.NE.0) THEN
     IF (meta_ionode) &
     READ (5, *, iostat = ios) (atomo (na), na = 1, &
          nat_todo)
     CALL mp_bcast(ios, meta_ionode_id, world_comm )
     CALL errore ('gwq_readin', 'reading atoms', ABS (ios) )
     CALL mp_bcast(atomo, meta_ionode_id, world_comm )
  ENDIF
  IF (nrapp.LT.0.OR.nrapp.GT.3 * nat) CALL errore ('gwq_readin', &
       'nrapp is wrong', 1)
  IF (nrapp.NE.0) THEN
     IF (meta_ionode) &
     READ (5, *, iostat = ios) (list (na), na = 1, nrapp)
     CALL mp_bcast(ios, meta_ionode_id, world_comm )
     CALL errore ('gwq_readin', 'reading list', ABS (ios) )
     CALL mp_bcast(list, meta_ionode_id, world_comm )
  ENDIF
  IF (ldisp .AND. (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
      CALL errore('gwq_readin','nq1, nq2, and nq3 must be greater than 0',1)

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
  config_green%max_iter = maxter_green
  config_green%threshold = tr2_green
  !
  CALL mp_bcast(priority_green, meta_ionode_id, world_comm)
  num_priority = COUNT(priority_green /= no_solver)
  !
  IF (num_priority == 0) THEN
    CALL errore(__FILE__, "priority for solvers not specified", 1)
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

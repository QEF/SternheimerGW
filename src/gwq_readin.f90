!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE gwq_readin()
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program GW.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE kinds,         ONLY : DP
  USE parameters,    ONLY : nsx
  USE constants,     ONLY : amconv, RYTOEV
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_bcast
  USE input_parameters, ONLY : max_seconds
  USE ions_base,     ONLY : amass, pmass, atm
  USE klist,         ONLY : xk, nks, nkstot, lgauss, two_fermi_energies
  USE control_flags, ONLY : gamma_only, tqr, restart, lkpoint_dir
  USE uspp,          ONLY : okvan
  USE fixed_occ,     ONLY : tfixed_occ
  USE lsda_mod,      ONLY : lsda, nspin
  USE printout_base, ONLY : title
  USE control_gw,    ONLY : maxter, alpha_mix, lgamma, lgamma_gamma, epsil, &
                            zue, zeu,  &
                            trans, reduce_io, tr2_gw, niter_gw, tr2_green, &
                            nmix_gw, ldisp, recover, lrpa, lnoloc, start_irr, &
                            last_irr, start_q, last_q, current_iq, tmp_dir_gw, &
                            ext_recover, ext_restart, u_from_file, modielec, eta, &
                            do_coulomb, do_sigma_c, do_sigma_exx, do_green, do_sigma_matel, &
                            do_q0_only, maxter_green

  USE save_gw,       ONLY : tmp_dir_save
  USE gamma_gamma,   ONLY : asr
  USE qpoint,        ONLY : nksq, xq
  USE partial,       ONLY : atomo, list, nat_todo, nrapp
  USE output,        ONLY : fildyn, fildvscf, fildrho
  USE disp,          ONLY : nq1, nq2, nq3, iq1, iq2, iq3, xk_kpoints, kpoints, num_k_pts, w_of_q_start
  USE io_files,      ONLY : tmp_dir, prefix, trimcheck
  USE noncollin_module, ONLY : i_cons, noncolin
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : iverbosity, modenum
  USE io_global,     ONLY : ionode, stdout
  USE mp_global,     ONLY : nproc, nproc_pool, nproc_file, nproc_pool_file, &
                            nimage, my_image_id,    &
                            nproc_image_file, nproc_image, mp_global_end, mpime, &
                            mp_barrier
  USE control_flags, ONLY : twfcollect
  USE paw_variables, ONLY : okpaw
! HL USE ramanm,        ONLY : eth_rps, eth_ns, lraman, elop, dek

  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, wsigmamin, wsigmamax, deltaw, wcoulmax
  USE gw_restart,    ONLY : gw_readfile
  USE gwsigma,       ONLY : ecutsig, nbnd_sig, ecutsex, ecutsco
  USE gwsymm,        ONLY : use_symm
  !
  !
  IMPLICIT NONE
  !
  INTEGER :: ios, ipol, iter, na, it, ierr
  ! integer variable for I/O control
  ! counter on polarizations
  ! counter on iterations
  ! counter on atoms
  ! counter on types
  REAL(DP) :: amass_input(nsx)
  ! save masses read from input here
  CHARACTER (LEN=256) :: outdir
  !
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=1), EXTERNAL :: capital
  CHARACTER(LEN=6) :: int_to_char
  INTEGER                    :: i
  LOGICAL                    :: nogg
  INTEGER, EXTERNAL  :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  LOGICAL, EXTERNAL  :: imatches
  !
  NAMELIST / INPUTGW / tr2_gw, amass, alpha_mix, niter_gw, nmix_gw,  &
                       nat_todo, iverbosity, outdir, epsil,  &
                       trans, zue, zeu, nrapp, max_seconds, reduce_io, &
                       modenum, prefix, fildyn, fildvscf, fildrho,   &
                       ldisp, nq1, nq2, nq3, iq1, iq2, iq3,   &
                       recover, fpol, asr, lrpa, lnoloc, start_irr, last_irr, &
                       start_q, last_q, nogg, modielec, ecutsig, nbnd_sig, eta, kpoints,&
                       ecutsco, ecutsex, do_coulomb, do_sigma_c, do_sigma_exx, do_green,& 
                       do_sigma_matel, tr2_green, do_q0_only, wsigmamin, wsigmamax, deltaw, wcoulmax,&
                       use_symm, maxter_green, w_of_q_start

  ! HL commented these vars in Namelist: eth_rps, eth_ns, lraman, elop, dek 
  ! tr2_ph       : convergence threshold
  ! amass        : atomic masses
  ! alpha_mix    : the mixing parameter
  ! niter_gw     : maximum number of iterations
  ! nmix_gw      : number of previous iterations used in mixing
  ! nat_todo     : number of atom to be displaced
  ! iverbosity   : verbosity control
  ! outdir       : directory where input, output, temporary files reside
  ! epsil        : if true calculate dielectric constant
  ! trans        : if true calculate phonon
  ! elph         : if true calculate electron-phonon coefficients
  ! zue          : if .true. calculate effective charges ( d force / dE )
  ! zeu          : if .true. calculate effective charges ( d P / du )
  ! lraman       : if true calculate raman tensor
  ! elop         : if true calculate electro-optic tensor
  ! nrapp        : the representations to do
  ! max_seconds  : maximum cputime for this run
  ! reduce_io    : reduce I/O to the strict minimum
  ! modenum      : single mode calculation
  ! prefix       : the prefix of files produced by pwscf
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! eth_rps      : threshold for calculation of  Pc R |psi> (Raman)
  ! eth_ns       : threshold for non-scf wavefunction calculation (Raman)
  ! dek          : delta_xk used for wavefunctions derivation (Raman)
  ! recover      : recover=.true. to restart from an interrupted run
  ! asr          : in the gamma_gamma case apply acoustic sum rule
  ! start_q      : in q list does the q points from start_q to last_q
  ! last_q       : 
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     : 
  ! nogg         : if .true. lgamma_gamma tricks are not used

  IF (ionode) THEN

  !
  ! flib/inpfile.f90!
  ! Reads in from standar input (5)
  ! 
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
     READ( 5, '(A)', IOSTAT = ios ) title
  !
  ENDIF
  !
  CALL mp_bcast(ios, ionode_id )
  CALL errore( 'gwq_readin', 'reading title ', ABS( ios ) )
  !
  ! Rewind the input if the title is actually the beginning of inputgw namelist
  IF( imatches("&inputgw", title)) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    REWIND(5, iostat=ios)
    CALL errore('gwq_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  ! ... set default values for variables in namelist
  !
  tr2_gw       = 1.D-9
  tr2_green    = 1.D-5
  amass(:)     = 0.D0
  alpha_mix(:) = 0.D0
  alpha_mix(1) = 0.7D0
  niter_gw     = maxter
  nmix_gw      = 3
  nat_todo     = 0
  modenum      = 0
  nrapp        = 0
  iverbosity   = 0
  trans        = .TRUE.
  lnoloc       = .FALSE.
  epsil        = .FALSE.
  zeu          = .TRUE.
  zue          = .FALSE.
  fpol         = .FALSE.
  max_seconds  =  1.E+7_DP
  reduce_io    = .FALSE.
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix       = 'pwscf'
  fildyn       = 'matdyn'
  fildrho      = ' '
  fildvscf     = ' '
  nq1          = 0
  nq2          = 0
  nq3          = 0
  iq1          = 0
  iq2          = 0
  iq3          = 0
  nogg         = .FALSE.
  recover      = .FALSE.
  asr          = .FALSE.
  start_irr    = 0
  last_irr     = -1000
  start_q      = 1
  last_q       =-1000
  ldisp        = .FALSE.
  lrpa         = .FALSE.
  maxter_green = 200

!Sigma cutoff, correlation cutoff, exchange cutoff
  ecutsig      = 2.5
  ecutsco      = 2.5 
  ecutsex      = 5.0
  nbnd_sig     = 8
  modielec     = .FALSE.

!imaginary component added to linear system should be in Rydberg
!eta            = 0.6/RYTOEV
!THIS IS half a volt you idiot!
  eta            = 0.04
  kpoints        = .FALSE.
  do_coulomb     = .FALSE.
  do_sigma_c     = .FALSE.
  do_sigma_exx   = .FALSE.
  do_green       = .FALSE.
  do_sigma_matel = .FALSE.
  do_q0_only     = .FALSE.

!Frequency variables
  wsigmamin      = -10.0d0
  wsigmamax      =  10.0d0
  deltaw         =   0.25d0 
  wcoulmax       = 110.0d0   

 !Symmetry Default:yes!, which q, point to start on.
 !can be used in conjunction with do_q0_only.
  use_symm       = .TRUE.
  w_of_q_start   = 1

  

  ! ...  reading the namelist inputgw

  IF (ionode) READ( 5, INPUTGW, IOSTAT = ios )


!HL TEST PARA FINE

   CALL mp_bcast(ios, ionode_id)
  !
   CALL errore( 'gwq_readin', 'reading inputgw namelist', ABS( ios ) )
  !
  IF (ionode) tmp_dir = trimcheck (outdir)
  CALL bcast_gw_input ( ) 
  CALL mp_bcast(nogg, ionode_id )

!HL FINE

  !
  ! ... Check all namelist variables
  !
  IF (tr2_gw <= 0.D0) CALL errore (' gwq_readin', ' Wrong tr2_gw ', 1)

  !HL raman thresholds
  !IF (eth_rps<= 0.D0) CALL errore ( 'gwq_readin', ' Wrong eth_rps', 1)
  !IF (eth_ns <= 0.D0) CALL errore ( 'gwq_readin', ' Wrong eth_ns ', 1)

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

  !  IF (dek <= 0.d0) CALL errore ( 'phq_readin', ' Wrong dek ', 1)
  !  epsil = epsil .OR. lraman .OR. elop
  !  IF ( (lraman.OR.elop) .AND. fildrho == ' ') fildrho = 'drho'
  !
  !  We can calculate  dielectric, raman or elop tensors and no Born effective
  !  charges dF/dE, but we cannot calculate Born effective charges dF/dE
  !  without epsil.

  !
  IF (zeu) zeu = epsil
  !
  !
  IF (ionode) THEN
     IF (.NOT. ldisp) &
        READ (5, *, iostat = ios) (xq (ipol), ipol = 1, 3)
  END IF
  CALL mp_bcast(ios, ionode_id)
  CALL errore ('gwq_readin', 'reading xq', ABS (ios) )
  CALL mp_bcast(xq, ionode_id )
  IF (.NOT.ldisp) THEN
     lgamma = xq (1) .EQ.0.D0.AND.xq (2) .EQ.0.D0.AND.xq (3) .EQ.0.D0
     IF ( (epsil.OR.zue) .AND..NOT.lgamma) CALL errore ('gwq_readin', &
          'gamma is needed for elec.field', 1)
  ENDIF

!  HL - Commenting the born effective charge stuff again
!  IF (zue.AND..NOT.trans) CALL errore ('phq_readin', 'trans must be &
!       &.t. for Zue calc.', 1)
!  IF (trans.AND.(lrpa.OR.lnoloc)) CALL errore('phq_readin', &
!                    'only dielectric constant with lrpa or lnoloc',1)
!  IF (lrpa.or.lnoloc) THEN
!     zeu=.FALSE.
!     lraman=.FALSE.
!     elop = .FALSE.
!  ENDIF

! HL here we can just use this to readin the list of frequencies that we want to calculate
! Stored in array  fiu(:), of size nfs.
! reads the frequencies ( just if fpol = .true. )

  IF ( fpol ) THEN
     nfs=0
     IF (ionode) THEN
        READ (5, *, iostat = ios) card
        IF ( TRIM(card)=='FREQUENCIES'.OR. &
             TRIM(card)=='frequencies'.OR. &
             TRIM(card)=='Frequencies') THEN
           READ (5, *, iostat = ios) nfs
        ENDIF
     ENDIF

     CALL mp_bcast(ios, ionode_id )
     CALL errore ('gwq_readin', 'reading number of FREQUENCIES', ABS(ios) )
     CALL mp_bcast(nfs, ionode_id )

     if (nfs > nfsmax) call errore('gwq_readin','Too many frequencies',1) 
     if (nfs < 1) call errore('gwq_readin','Too few frequencies',1) 

     IF (ionode) THEN
        IF ( TRIM(card) == 'FREQUENCIES' .OR. &
             TRIM(card) == 'frequencies' .OR. &
             TRIM(card) == 'Frequencies' ) THEN
           DO i = 1, nfs
              !HL Need to convert frequencies from electron volts into Rydbergs
              READ (5, *, iostat = ios) fiu(i)
              fiu(i) = fiu(i) / RYTOEV
           END DO
        END IF
     END IF

     CALL mp_bcast(ios, ionode_id)
     CALL errore ('gwq_readin', 'reading FREQUENCIES card', ABS(ios) )
     CALL mp_bcast(fiu, ionode_id )

  ELSE
     nfs=0
     fiu=0.0_DP
  END IF

! Reading in kpoints specified by user.
! Note max number of k-points is 10. 
! Why? Because that's the number I picked. 
! If k-points option is not specified it defaults to Gamma. 

 IF (kpoints) then
     num_k_pts = 0
     IF (ionode) THEN
        READ (5, *, iostat = ios) card
        READ (5, *, iostat = ios) card
        IF ( TRIM(card)=='K_POINTS'.OR. &
             TRIM(card)=='k_points'.OR. &
             TRIM(card)=='K_points') THEN
           READ (5, *, iostat = ios) num_k_pts
        ENDIF
     ENDIF
     CALL mp_bcast(ios, ionode_id )
     CALL errore ('gwq_readin', 'reading number of kpoints', ABS(ios) )
     CALL mp_bcast(num_k_pts, ionode_id )
     if (num_k_pts > 10) call errore('gwq_readin','Too many k-points',1) 
     if (num_k_pts < 1) call errore('gwq_readin','Too few kpoints',1) 
     IF (ionode) THEN
        IF ( TRIM(card)=='K_POINTS'.OR. &
             TRIM(card)=='k_points'.OR. &
             TRIM(card)=='K_points') THEN
           DO i = 1, num_k_pts
           !DO i = 1, 2
              !HL Need to convert frequencies from electron volts into Rydbergs
              READ (5, *, iostat = ios) xk_kpoints(1,i), xk_kpoints(2,i), xk_kpoints(3,i)
              !write(6,'(3f11.7)') xk_kpoints(:,i)
           END DO
        END IF
     END IF
     CALL mp_bcast(ios, ionode_id)
     CALL errore ('gwq_readin', 'reading FREQUENCIES card', ABS(ios) )
     CALL mp_bcast(xk_kpoints, ionode_id )
 ELSE
     num_k_pts = 1
 ENDIF
  !
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  amass_input(:)= amass(:)
  !
  tmp_dir_save=tmp_dir
  tmp_dir_gw= TRIM (tmp_dir) // '_gw' // int_to_char(my_image_id)
  ext_restart=.FALSE.
  ext_recover=.FALSE.

  IF (recover) THEN
     CALL gw_readfile('init',ierr)
     IF (ierr /= 0 ) THEN
        recover=.FALSE.
        goto 1001
     ENDIF
     tmp_dir=tmp_dir_gw
     CALL check_restart_recover(ext_recover, ext_restart)
     tmp_dir=tmp_dir_save
     IF (ldisp) lgamma = (current_iq==1)
!
!  If there is a restart or a recover file gw.x has saved its own data-file 
!  and we read the initial information from that file
!
     IF ((ext_recover.OR.ext_restart).AND..NOT.lgamma) &
                                                      tmp_dir=tmp_dir_gw
     u_from_file=.true.
  ENDIF
1001 CONTINUE

  ! HL !!ATTENZIONE!! This is where the files from the SCF step are read.
  ! QE description:
  ! read_file reads the variables defining the system
  ! in parallel execution, only root proc reads the file
  ! and then broadcast the values to all other procs
  ! distribute across pools k-points and related variables.
  ! nks is defined by divide et_impera.
  ! It also reads in wavefunctions in "distributed form" i.e. split over k-points...


  CALL read_file ( )
  tmp_dir=tmp_dir_save
  !
  IF (modenum > 3*nat) CALL errore ('gwq_readin', ' Wrong modenum ', 2)

  IF (gamma_only) CALL errore('gwq_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)

!HL
!  IF (lda_plus_u) CALL errore('gwq_readin',&
!     'The phonon code with LDA+U is not yet available',1)

!  IF (okpaw.and.(lraman.or.elop.or.elph)) CALL errore('phq_readin',&
!     'The phonon code with paw and raman, elop or elph is not yet available',1)

!  IF (okvan.and.(lraman.or.elop)) CALL errore('phq_readin',&
!     'The phonon code with US-PP and raman or elop not yet available',1)

!  IF (noncolin.and.(lraman.or.elop.or.elph)) CALL errore('phq_readin', &
!      'lraman, elop, or e-ph and noncolin not programed',1)


  IF (nproc_image /= nproc_image_file .and. .not. twfcollect)  &
     CALL errore('gwq_readin',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('gwq_readin',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)


!  IF (elph.and.nimage>1) CALL errore('phq_readin',&
!     'elph with image parallelization is not yet available',1)

!  IF (two_fermi_energies.or.i_cons /= 0) &
!     CALL errore('gwq_readin',&
!     'The phonon code with constrained magnetization is not yet available',1)

  IF (tqr) CALL errore('gwq_readin',&
     'The phonon code with Q in real space not available',1)

  IF (start_irr < 0 ) CALL errore('gwq_readin', 'wrong start_irr',1)
  !
  IF (start_q <= 0 ) CALL errore('gwq_readin', 'wrong start_q',1)
  !
  !
  ! If a band structure calculation needs to be done do not open a file 
  ! for k point
  !
  lkpoint_dir=.FALSE.
  restart = recover
  !
  !  set masses to values read from input, if available;
  !  leave values read from file otherwise
  !
  DO it = 1, ntyp
     IF (amass_input(it) < 0.0_DP) amass_input(it)= &
              atom_weight(atomic_number(TRIM(atm(it))))
     IF (amass_input(it) > 0.D0) amass(it) = amass_input(it)
     IF (amass(it) <= 0.D0) CALL errore ('gwq_readin', 'Wrong masses', it)
     !
     !  convert masses to a.u.
     !
     pmass(it) = amconv * amass(it)
  ENDDO


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

!HL
!    IF (tfixed_occ) &
!     CALL errore('phq_readin','phonon with arbitrary occupations not tested',1)
!  !
!  IF (elph.AND..NOT.lgauss) CALL errore ('phq_readin', 'Electron-&
!       &phonon only for metals', 1)
!  IF (elph.AND.lsda) CALL errore ('phq_readin', 'El-ph and spin not &
!       &implemented', 1)
!  IF (elph.AND.fildvscf.EQ.' ') CALL errore ('phq_readin', 'El-ph needs &
!       &a DeltaVscf file', 1)
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here

  ! HL lets comment this allocate part function
  ! allocate_part := dynamical allocation of arrays for the control of partial computation
  ! of the dynamical matrix
  ! CALL allocate_part ( nat )
  ! HL

  IF ( nat_todo < 0 .OR. nat_todo > nat ) &
     CALL errore ('gwq_readin', 'nat_todo is wrong', 1)
  IF (nat_todo.NE.0) THEN
     IF (ionode) &
     READ (5, *, iostat = ios) (atomo (na), na = 1, &
          nat_todo)
     CALL mp_bcast(ios, ionode_id )
     CALL errore ('gwq_readin', 'reading atoms', ABS (ios) )
     CALL mp_bcast(atomo, ionode_id )
  ENDIF
  IF (nrapp.LT.0.OR.nrapp.GT.3 * nat) CALL errore ('gwq_readin', &
       'nrapp is wrong', 1)
  IF (nrapp.NE.0) THEN
     IF (ionode) &
     READ (5, *, iostat = ios) (list (na), na = 1, nrapp)
     CALL mp_bcast(ios, ionode_id )
     CALL errore ('gwq_readin', 'reading list', ABS (ios) )
     CALL mp_bcast(list, ionode_id )
  ENDIF
  
  !HL Commmenting various epsil/lgaus stuff.
  ! IF (epsil.AND.lgauss) &
  !      CALL errore ('phq_readin', 'no elec. field with metals', 1)
  !IF (modenum > 0) THEN
  !   IF ( ldisp ) &
  !        CALL errore('phq_readin','Dispersion calculation and &
  !        & single mode calculation not possibile !',1)
  !   nrapp = 1
  !   nat_todo = 0
  !   list (1) = modenum
  ! ENDIF
  !IF (modenum > 0 .OR. lraman ) lgamma_gamma=.FALSE.
  IF (.NOT.lgamma_gamma) asr=.FALSE.
  !
   IF (ldisp .AND. (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
       CALL errore('phq_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  RETURN
  !
END SUBROUTINE gwq_readin

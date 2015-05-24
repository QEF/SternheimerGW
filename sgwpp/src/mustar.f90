MODULE units_coulmat
  ! ... the units of the files and the record lengths
  SAVE
  INTEGER :: iuncoulmat, iuncoul
  INTEGER :: lrcoulmat, lrcoul
END MODULE units_coulmat

MODULE dielectric
  USE kinds,    ONLY : DP
  SAVE
  REAL(DP) :: qtf, kf
END MODULE dielectric

MODULE control_coulmat
  USE kinds,    ONLY : DP
  SAVE
  LOGICAL  :: do_coulmat, do_fsavg, do_lind, do_diag, do_plotmuk, do_dosband
  INTEGER  :: nbndmin, nbndmax, ngcoul
  REAL(DP) :: degaussfs, debye_e
END MODULE control_coulmat
!------------------------------
PROGRAM mustar
!------------------------------
  !
#  define DIRECT_IO_FACTOR 8 
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id, meta_ionode 
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  USE iotk_module
  USE xml_io_base
  USE io_files,    ONLY : tmp_dir, prefix, outdir, diropn
  USE constants,   ONLY : RYTOEV, eps8
  USE ener,        ONLY : ef
  USE klist,       ONLY : lgauss
  USE ktetra,      ONLY : ltetra
  USE wvfct,       ONLY : nbnd
  USE lsda_mod,    ONLY : nspin
  USE mp_global,   ONLY : mp_startup
  USE environment, ONLY : environment_start
  USE klist,       ONLY : nks, nkstot, degauss,xk,wk
  USE start_k,     ONLY : nks_start, xk_start, wk_start, &
                          nk1, nk2, nk3, k1, k2, k3
  USE units_coulmat,   ONLY : iuncoulmat, lrcoulmat, lrcoul, iuncoul
  USE dielectric,      ONLY : qtf, kf
  USE control_coulmat, ONLY : do_coulmat, do_fsavg, nbndmin, degaussfs, ngcoul, do_lind, &
                              debye_e, do_diag, do_plotmuk, do_dosband, nbndmax
  USE mp_global,        ONLY : inter_image_comm, intra_image_comm, &
                               my_image_id, nimage, root_image
  IMPLICIT NONE
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! input variables
  !
  INTEGER                 :: nw,nshell,ibndmin,ibndmax
  INTEGER                 :: nk1tmp, nk2tmp, nk3tmp
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift,eta, qmod_par
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc, exst
  !
  NAMELIST / inputpp / prefix, outdir, calculation, nk1, nk2, nk3, qtf, do_coulmat, &
                       do_fsavg, nbndmin, kf, degaussfs, ngcoul, do_lind, debye_e,  &
                       do_diag, do_plotmuk, do_dosband, nbndmax
  NAMELIST / energy_grid / smeartype,intersmear,intrasmear,wmax,wmin,nw,shift,nshell,eta,ibndmin,ibndmax, qmod_par
  !
  ! local variables
  !
  INTEGER :: ios

  COMPLEX(DP), ALLOCATABLE :: vcnknpkp(:,:,:,:)
  character(len=256) :: tempfile, filename
  INTEGER :: recl
  integer*8 :: unf_recl
  

!---------------------------------------------
! program body
!---------------------------------------------
!
  ! initialise environment
  !
#ifdef __MPI
!  CALL mp_startup ()
  CALL mp_startup (start_images=.true.)
#endif
  CALL environment_start ( 'coulmatel' )
  !
  ! Set default values for variables in namelist
  !
  calculation  = 'coulmat'
  prefix       = 'pwscf'
  shift        = 0.0d0
  outdir       = './'
  intersmear   = 0.136
  wmin         = 0.0d0
  wmax         = 30.0d0
  eta          = 0.3
  nbndmin      = 1
  nbndmax      = 8
  ibndmin      = 1
  nshell       = 0
  ibndmax      = 0
  qmod_par     = 0.5
  nw           = 600
  smeartype    = 'gauss'
  intrasmear   = 0.0d0
  metalcalc    = .false.
  do_coulmat   = .true.
  do_fsavg     = .true.
  do_lind      = .false.
  qtf          = 1.0
  kf           = 1.0
  degaussfs    = 0.05
  debye_e      = 0.0023
  do_diag      =.false.

  !SHOULD READ FROM INPUT FILE

  IF(ionode) then
    CALL input_from_file( )
    IF (ionode) WRITE( stdout, "( 2/, 5x, 'Reading input file...' ) " )
    ios = 0
    IF ( ionode ) READ (5, inputpp, IOSTAT=ios)
  ENDIF
  !
  ! read input file
  !

  nk1tmp = nk1
  nk2tmp = nk2
  nk3tmp = nk3
  !
  CALL mp_bcast ( ios, ionode_id, world_comm )
  !
  IF (ios/=0) CALL errore('epsilon', 'reading namelist INPUTPP', abs(ios))

  IF ( ionode ) THEN
     tmp_dir = trimcheck(outdir)
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id, world_comm )
  IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', abs(ios))
  !
  ! ... Broadcast variables
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Broadcasting variables...' ) " )
  CALL mp_bcast( smeartype, ionode_id, world_comm )
  CALL mp_bcast( calculation, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( shift, ionode_id, world_comm )
  CALL mp_bcast( outdir, ionode_id, world_comm )
  CALL mp_bcast( intrasmear, ionode_id, world_comm )
  CALL mp_bcast( intersmear, ionode_id, world_comm )
  CALL mp_bcast( wmax, ionode_id, world_comm )
  CALL mp_bcast( wmin, ionode_id, world_comm )
  CALL mp_bcast( nw, ionode_id, world_comm      )
  CALL mp_bcast( nbndmin, ionode_id, world_comm )
  CALL mp_bcast( nbndmax, ionode_id, world_comm )
  CALL mp_bcast( nshell, ionode_id, world_comm  )
  CALL mp_bcast( eta, ionode_id, world_comm     )
  CALL mp_bcast( do_coulmat, ionode_id, world_comm  )
  CALL mp_bcast( do_fsavg,   ionode_id, world_comm  )
  CALL mp_bcast( do_lind,   ionode_id, world_comm  )
  CALL mp_bcast( kf,   ionode_id, world_comm  )
  CALL mp_bcast( qtf,   ionode_id, world_comm  )
  CALL mp_bcast( do_diag,   ionode_id, world_comm  )
  CALL mp_bcast( do_plotmuk,   ionode_id, world_comm  )
  CALL mp_bcast( do_dosband,   ionode_id, world_comm  )

  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )

  CALL read_file

! CALL read_file_coul

! openfil_pp set twfcollet to false and opens wavefunctions.
! it assumes wave functions are \emph{always} required 
! in distributed format.
  CALL openfil_pp

!nk1 etc. gets over written by readfile.
  nk1 = nk1tmp
  nk2 = nk2tmp
  nk3 = nk3tmp

  CALL mp_bcast( nk1, ionode_id, world_comm )
  CALL mp_bcast( nk2, ionode_id, world_comm )
  CALL mp_bcast( nk3, ionode_id, world_comm )

  IF (ionode) WRITE(stdout,"(2/, 5x, 'Fermi energy [eV] is: ',f8.5)") ef*RYTOEV

  IF (lgauss .or. ltetra) THEN
      metalcalc=.true.
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal...' ) " )
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a dielectric...' ) " )
  ENDIF

!  IF (nbndmax == 0) nbndmax = nbnd
!  IF (ibndmax == 0) ibndmax = nbnd

  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") trim(calculation)


  SELECT CASE ( trim(calculation) )
  !
  CASE ( 'coulmat' )

 ALLOCATE(vcnknpkp(nks,nks,nbnd,nbnd))
 vcnknpkp = (0.0d0,0.d0)

 iuncoulmat = 29
 lrcoulmat  = nks 
 iuncoul    = 28
 lrcoul     = 2 * ngcoul * ngcoul

 CALL mp_bcast( iuncoul, ionode_id, world_comm )
 CALL mp_bcast( lrcoul,  ionode_id, world_comm )
 CALL mp_bcast( ngcoul,  ionode_id, world_comm )
 CALL mp_bcast( iuncoulmat,  ionode_id, world_comm )

!OPEN DIRECTORY FOR COULOMB MATRIX ELEMENTS
 IF(ionode) THEN
    CALL diropn (iuncoulmat, 'coulmat', lrcoulmat, exst)
 ENDIF
!OPEN COULOMB directory
  !CALL diropn (iuncoul, 'coul', lrcoul, exst)
  !tempfile = trim(tmp_dir) // trim('mgb2.coul1')
  tempfile = trim(tmp_dir) // trim(prefix) // trim('.coul1')
  unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
  open (iuncoul, file = trim(adjustl(tempfile)), iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)

  if (ios /= 0) call errore ('diropn', 'error opening '//trim(tempfile), iuncoul)
! ENDIF

 WRITE(stdout, '(\5x, "nk1 nk2 nk3", 3i4)'), nk1, nk2, nk3
 WRITE(stdout, '(\5x, "Fermi Vector", f12.7)'), kf
 WRITE(stdout, '(\5x, "Thomas-Fermi Vector", f12.7)'), qtf
!Calculate Coulomb matrix elements stored in struct vcnknpkp
  CALL start_clock( 'calculation' )
  IF(do_coulmat) THEN
    IF (ionode) WRITE( stdout, "( 5x, 'Calculating Coulomb Matrix Elements' ) " )
        CALL coulmatsym()
!        CALL coulmatstar()
  ELSE IF(do_dosband) THEN
    IF (ionode) WRITE( stdout, "( 5x, 'Calculating Dosband' ) " )
        CALL dosband()
  ENDIF
  
  IF (do_plotmuk) THEN
     IF(meta_ionode)   CALL plotmuk() 
  ENDIF
  call mp_barrier(inter_image_comm)

  CASE DEFAULT
      CALL errore('sgwpp','invalid CALCULATION = '//trim(calculation),1)
  END SELECT
  CALL stop_clock( 'calculation' )
  CALL print_clock( 'calculation' )
  IF ( ionode ) WRITE( stdout, *  )
  CALL stop_pp ()
END PROGRAM mustar



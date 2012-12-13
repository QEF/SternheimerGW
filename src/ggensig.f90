  !
  !----------------------------------------------------------------
  SUBROUTINE ggensig ()
  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------
  ! FIND THE G-VECTORS FOR THE SMALL SIGMA CUTOFF
  !----------------------------------------------------------------
  !
  ! the G^2 cutoff in units of 2pi/a_0
  ! Note that in Ry units the kinetic energy is G^2, not G^2/2
  !
  ! FG
  ! determine G-vectors within the cutoff from the
  ! array already created in ggen
  ! HL basically the same as in SGWI except i've had to change all the ngm(s) etc. into (sig) 
  ! TESTED THIS DOES PRODUCE SAME ORDERING A ggen.f90 4.2.1 if the same cutoff is used. 
  ! nls is already reserved for  the smooth grid in Quantum espresso. 
  ! Also added variable ecutsig which can be defined by the user at input.

  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE gvect,            ONLY : gcutm, ecutwfc, dual, nr1, nr2, nr3, ngm, g, igtongl, gl, nl
  USE gwsigma,          ONLY : gcutmsig, nlsig, ngmsig, nr1sig, nr2sig, nr3sig, nrsig, ecutsig, &
                               nlsex, nlsco
  USE cell_base,        ONLY : at, tpiba2
  USE fft_scalar,       ONLY : allowed
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE qpoint,           ONLY : xq
  USE control_gw,      ONLY : done_bands, reduce_io, recover, tmp_dir_gw, &
                              ext_restart, bands_computed
  USE save_gw,         ONLY : tmp_dir_save
  USE input_parameters,ONLY : pseudo_dir
  USE io_files,        ONLY : prefix, tmp_dir
  USE control_flags,   ONLY : restart
  
  IMPLICIT NONE

  integer :: n1, n2, n3, i, j, k, ipol, ig, igl, ng

  !HL ecutsig defined in punch card
  !cutmsig = 4.D0 * ecutsig / tpiba2

!HL
  gcutmsig = 4.D0 * ecutsig / tpiba2

  nr1sig = 1 + int (2 * sqrt (gcutmsig) * sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) )
  nr2sig = 1 + int (2 * sqrt (gcutmsig) * sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) )
  nr3sig = 1 + int (2 * sqrt (gcutmsig) * sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) )

  do while (.not.allowed(nr1sig))
    nr1sig = nr1sig + 1
  enddo

  do while (.not.allowed(nr2sig))
    nr2sig = nr2sig + 1
  enddo

  do while (.not.allowed(nr3sig))
    nr3sig = nr3sig + 1
  enddo
 
  xq(1) = 0.0d0
  xq(2) = 0.0d0
  xq(3) = 0.0d0 


  ! ... arrays allocated in input.f90, read_file.f90 or setup.f90
  CALL clean_pw( .FALSE. )

  !From now on, work only on the _gw virtual directory

  tmp_dir=tmp_dir_gw

!  write(6,*)tmp_dir_gw

  ! ... Setting the values for the nscf run

  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  restart = ext_restart
  pseudo_dir= TRIM( tmp_dir_save ) // TRIM( prefix ) // '.save'

  CALL setup_nscf(xq)
  CALL init_run()

  !Loop over all g vectors in (order of size) and test that |G|.le.ngmsig.
  !igtongl should only be defined when init_run is called in pwscf. No idea
  !where the array values would come from otherwise...
  !write(6,*)igtongl(:)
  !stop
  
  do ng = 1, ngm
    if ( gl( igtongl (ng) ) .le. gcutmsig ) ngmsig = ng
   !write(6,*) igtongl(ng) 
  enddo

  ALLOCATE ( nlsig(ngmsig) )

!HL upping cutoff of nlsig
  !ALLOCATE ( nlsig(ngmsig) )
  !
  ! Now set nl with the correct fft correspondence
  !

  do ng = 1, ngmsig
     ! n1 is going to be i+1, folded to positive when <= 0
     n1 = nint (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, ng) * at (3, 1) ) + 1
     if (n1.lt.1) n1 = n1 + nr1sig
     ! n2 is going to be j+1, folded to positive when <= 0
     n2 = nint (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, ng) * at (3, 2) ) + 1
     if (n2.lt.1) n2 = n2 + nr2sig
     ! n3 is going to be k+1, folded to positive when <= 0
     n3 = nint (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, ng) * at (3, 3) ) + 1
     if (n3.lt.1) n3 = n3 + nr3sig
     !
     if (n1.le.nr1sig.and.n2.le.nr2sig.and.n3.le.nr3sig) then
       nlsig (ng) = n1 + (n2 - 1) * nr1sig + (n3 - 1) * nr1sig * nr2sig
       !write (6,*)nlsig(ng)
     else
        !call error('ggens','Mesh too small?',ng)
        WRITE(6,'("ERROR: Mesh too small?")')
        STOP
     endif
  enddo

  !total number of real-space grid points
  !write(6,*) nl(:)
  !write(6,*) nlsig(:)

  nrsig = nr1sig * nr2sig * nr3sig

  write(6,'(4x,"")')
  write(6,'(4x,"ngmsig = ",i10)') ngmsig
  write(6,'(4x,"nr1sig = ",i10)') nr1sig
  write(6,'(4x,"nr2sig = ",i10)') nr2sig
  write(6,'(4x,"nr3sig = ",i10)') nr3sig
  write(6,'(4x,"nrsig  = ",i10)') nrsig
  
  !
  END SUBROUTINE ggensig

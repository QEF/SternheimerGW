! Copyright (C) 2004-2009 Andrea Benassi and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!------------------------------
 MODULE grid_module
!------------------------------
  USE kinds,        ONLY : DP
  IMPLICIT NONE
  PRIVATE

  !
  ! general purpose vars
  !
  REAL(DP), ALLOCATABLE  :: focc(:,:), wgrid(:)
  REAL(DP)               :: alpha
  !
  !
  PUBLIC :: grid_build, grid_destroy
  PUBLIC :: focc, wgrid, alpha
  !
CONTAINS

!---------------------------------------------
  SUBROUTINE grid_build(nw, wmax, wmin)
  !-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE wvfct,     ONLY : nbnd, wg
  USE klist,     ONLY : nks, wk, nelec, xk
  USE lsda_mod,  ONLY : nspin
  USE uspp,      ONLY : okvan
  USE io_global, ONLY : stdout, ionode, ionode_id

  !
  IMPLICIT NONE
  !
  ! input vars
  INTEGER,  INTENT(in) :: nw
  REAL(DP), INTENT(in) :: wmax ,wmin
  !
  ! local vars
  INTEGER         :: iw,ik,i,ierr

  !
  ! check on the number of bands: we need to include empty bands in order to allow
  ! to write the transitions
  !
  !print*,  nelec, nbnd 
  IF ( REAL(nbnd, DP)  <= nelec / 2.0_DP ) CALL errore('epsilon', 'bad band number', 1)

  !
  ! spin is not implemented
  !
  IF( nspin > 2 ) CALL errore('grid_build','Non collinear spin  calculation not implemented',1)

  !
  ! USPP are not implemented (dipole matrix elements are not trivial at all)
  !
  IF ( okvan ) CALL errore('grid_build','USPP are not implemented',1)

  ALLOCATE ( focc( nbnd, nks), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', abs(ierr))
  !
  ALLOCATE( wgrid( nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating wgrid', abs(ierr))

  !
  ! check on k point weights, no symmetry operations are allowed
  !
!  DO ik = 2, nks
!     !
!     IF ( abs( wk(1) - wk(ik) ) > 1.0d-8 ) &
!        CALL errore('grid_build','non unifrom kpt grid', ik )
!     !
!  ENDDO
  !
  ! PRINT the k point grid and weights as a test
  !
  !DO ik = 1, nks
  !   !
  !   IF ( ionode ) WRITE (6,*) xk(1,ik), xk(2,ik),xk(3,ik) , wk(ik)
  !   !
  !ENDDO
  !
  ! occupation numbers, to be normalized differently
  ! whether we are spin resolved or not
  !
  IF(nspin==1) THEN
    DO ik = 1,nks
    DO i  = 1,nbnd
         focc(i,ik)= wg(i, ik ) * 2.0_DP / wk( ik )
    ENDDO
    ENDDO
  ELSEIF(nspin==2) THEN
    DO ik = 1,nks
    DO i  = 1,nbnd
         focc(i,ik)= wg(i, ik ) * 1.0_DP / wk( ik )
    ENDDO
    ENDDO
  ENDIF
  !
  ! set the energy grid
  !
  alpha = (wmax - wmin) / REAL(nw, DP)
  !
  DO iw = 1, nw
      wgrid(iw) = wmin + iw * alpha
  ENDDO
  !
END SUBROUTINE grid_build
!
!
!----------------------------------
  SUBROUTINE grid_destroy
  !----------------------------------
  IMPLICIT NONE
  INTEGER :: ierr
  !
  IF ( allocated( focc) ) THEN
      !
      DEALLOCATE ( focc, wgrid, STAT=ierr)
      CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))
      !
  ENDIF
  !
END SUBROUTINE grid_destroy

END MODULE grid_module


!------------------------------
PROGRAM epsilon
!------------------------------
  !
  ! Compute the complex macroscopic dielectric function,
  ! at the RPA level, neglecting local field effects.
  ! Eps is computed both on the real or immaginary axis
  !
  ! Authors: Andrea Benassi, Andrea Ferretti, Carlo Cavazzoni
  !
  ! NOTE: Part of the basic implementation is taken from pw2gw.f90;
  !
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  USE iotk_module
  USE xml_io_base
  USE io_files,    ONLY : tmp_dir, prefix, outdir
  USE constants,   ONLY : RYTOEV
  USE ener,        ONLY : ef
  USE klist,       ONLY : lgauss
  USE ktetra,      ONLY : ltetra
  USE wvfct,       ONLY : nbnd
  USE lsda_mod,    ONLY : nspin
  USE mp_global,   ONLY : mp_startup
  USE environment, ONLY : environment_start
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! input variables
  !
  INTEGER                 :: nw,nbndmin,nbndmax,nshell,ibndmin,ibndmax
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift,eta, qmod_par
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype,intersmear,intrasmear,wmax,wmin,nbndmin,nbndmax,nw,shift,nshell,eta,ibndmin,ibndmax, qmod_par
  !
  ! local variables
  !
  INTEGER :: ios

!---------------------------------------------
! program body
!---------------------------------------------
!
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'epsilon' )
  !
  ! Set default values for variables in namelist
  !
  calculation  = 'eps'
  prefix       = 'pwscf'
  shift        = 0.0d0
  outdir       = './'
  intersmear   = 0.136
  wmin         = 0.0d0
  wmax         = 30.0d0
  eta          = 0.3
  nbndmin      = 1
  ibndmin      = 1
  nshell       = 0
  nbndmax      = 0
  ibndmax      = 0
  qmod_par     = 0.5
  nw           = 600
  smeartype    = 'gauss'
  intrasmear   = 0.0d0
  metalcalc    = .false.

  !
  ! this routine allows the user to redirect the input using -input
  ! instead of <
  !
  CALL input_from_file( )

  !
  ! read input file
  !
  IF (ionode) WRITE( stdout, "( 2/, 5x, 'Reading input file...' ) " )
  ios = 0
  !
  IF ( ionode ) READ (5, inputpp, IOSTAT=ios)
  !
  CALL mp_bcast ( ios, ionode_id )
  !
  IF (ios/=0) CALL errore('epsilon', 'reading namelist INPUTPP', abs(ios))
  !
  IF ( ionode ) THEN
     !
     READ (5, energy_grid, IOSTAT=ios)
     !
     tmp_dir = trimcheck(outdir)
     !
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id )
  IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', abs(ios))
  !
  ! ... Broadcast variables
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Broadcasting variables...' ) " )

  CALL mp_bcast( smeartype, ionode_id )
  CALL mp_bcast( calculation, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( shift, ionode_id )
  CALL mp_bcast( outdir, ionode_id )
  CALL mp_bcast( intrasmear, ionode_id )
  CALL mp_bcast( intersmear, ionode_id)
  CALL mp_bcast( wmax, ionode_id )
  CALL mp_bcast( wmin, ionode_id )
  CALL mp_bcast( nw, ionode_id )
  CALL mp_bcast( ibndmin, ionode_id )
  CALL mp_bcast( ibndmax, ionode_id )
  CALL mp_bcast( qmod_par, ionode_id )
  CALL mp_bcast( nbndmin, ionode_id )
  CALL mp_bcast( nbndmax, ionode_id )
  CALL mp_bcast( nshell, ionode_id )
  CALL mp_bcast( eta, ionode_id )

  !
  ! read PW simulation parameters from prefix.save/data-file.xml
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )

  CALL read_file
  CALL openfil_pp

  !
  ! few conversions
  !

  IF (ionode) WRITE(stdout,"(2/, 5x, 'Fermi energy [eV] is: ',f8.5)") ef *RYTOEV

  IF (lgauss .or. ltetra) THEN
      metalcalc=.true.
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal...' ) " )
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a dielectric...' ) " )
  ENDIF

  IF (nbndmax == 0) nbndmax = nbnd
  IF (ibndmax == 0) ibndmax = nbnd

  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") trim(calculation)

  CALL start_clock( 'calculation' )
  !
  SELECT CASE ( trim(calculation) )
  !
  CASE ( 'gwppa' )
      !
      CALL gwppa_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin,nshell,eta,ibndmin,ibndmax,qmod_par)
      !
  CASE ( 'eps' )
      !
      CALL eps_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin )
      !
  CASE ( 'jdos' )
      !
      CALL jdos_calc ( smeartype,intersmear,nw,wmax,wmin,nbndmin,nbndmax,shift,nspin )
      !
  CASE ( 'offdiag' )
      !
      CALL offdiag_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin )
      !
  CASE ( 'occ' )
      !
      CALL occ_calc ()
      !
  CASE DEFAULT
      !
      CALL errore('epsilon','invalid CALCULATION = '//trim(calculation),1)
      !
  END SELECT
  !
  CALL stop_clock( 'calculation' )

  !
  ! few info about timing
  !
  CALL stop_clock( 'epsilon' )
  !
  IF ( ionode ) WRITE( stdout , "(/)" )
  !
  CALL print_clock( 'epsilon' )
  CALL print_clock( 'calculation' )
  CALL print_clock( 'dipole_calc' )
  !
  IF ( ionode ) WRITE( stdout, *  )

  !
  !
  CALL stop_pp ()

END PROGRAM epsilon

!-----------------------------------------------------------------------------
SUBROUTINE gwppa_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax, shift, &
                      metalcalc , nspin, nshell,eta, ibndmin, ibndmax, qmod_par)
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV, e2
  USE cell_base,            ONLY : tpiba2, omega, at, alat
  USE wvfct,                ONLY : nbnd, et, npw, igk, npwx,  g2kin, ecutwfc
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss,xk,wk
  USE io_global,            ONLY : ionode, stdout
  !
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm,me_bgrp
  USE mp,                   ONLY : mp_sum
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE wavefunctions_module, ONLY : evc !,psic
  USE fft_base,             ONLY : dfftp, dffts
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : nls, nlsm 
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE funct,                ONLY : xc
  USE grid_subroutines,     ONLY : realspace_grids_info

  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin,nshell,ibndmin,ibndmax
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear,intrasmear, shift, eta, qmod_par
  LOGICAL,         INTENT(in) :: metalcalc
  !
  ! local variables
  !
  INTEGER       :: i,j,k, ik, iband1, iband2,is, jbnd, ibnd,ig, jg, kg, ir 
  INTEGER       :: iw, iwp, ierr, nw1, nksigma, sysd
  REAL(DP)      :: etrans, const, w, renorm(3), wplasma, wplasma0, wtilde0, eps0, eps1, q, kTF, qG, qG1
  COMPLEX(DP)   :: wtilde
  REAL(DP)      :: wmin1, wmax1, wksigma, epssum, epsinv, epssum1 
  REAL(DP)      :: vtxc, etxc, GN_freq, epsqdep 
  COMPLEX(DP)   :: vxc,vc1,ovlp
  COMPLEX(DP)               ::   ZDOTC
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  COMPLEX(DP)   :: prefactor
  REAL(DP) :: r(3), g2(3)
  
  !
  REAL(DP), ALLOCATABLE    :: epsr(:,:), epsi(:,:), epsrc(:,:,:), epsic(:,:,:)
  REAL(DP), ALLOCATABLE    :: ieps(:,:), eels(:,:), iepsc(:,:,:), eelsc(:,:,:), ieps_qdep(:,:) 
  REAL(DP), ALLOCATABLE    :: ieps_tmp (:,:)
  REAL(DP), ALLOCATABLE    :: dipole(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: sigmare (:)
!  COMPLEX(DP), sigmax
  REAL(DP), ALLOCATABLE    :: spectrum(:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_aux(:,:,:)
  COMPLEX(DP)              :: tempevc(dffts%nnr), tempevc1(dffts%nnr), tempevc2(dffts%nnr)
  COMPLEX(DP),ALLOCATABLE  :: psi(:),vpsi(:) !, tempevc(:)
  COMPLEX(DP)  :: evc1(npwx,nbnd)
  REAL(DP) :: rhox, arhox, zeta, amag, vs, ex, ec, vx(2), vc(2), rhoneg(2), qmin
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10
  INTEGER  :: npw1, ncount, indexqmin, indexgmin, npar, ipar, ipole
  real(DP), ALLOCATABLE :: par (:), par_tmp (:)
  INTEGER  :: igk1 (npwx)
  REAL(DP) :: absg0vec,ovlpsum
  INTEGER, ALLOCATABLE   :: gmap  (:,:)
  REAL(DP), ALLOCATABLE  :: g0vec (:,:)
  REAL(DP) :: g0 (3), gtmp(3,ngm)
  REAL(DP) :: eta_large(2)
  INTEGER  :: ip, index, index0, ir_end, ng0vec
  CHARACTER*2 label
  CHARACTER*12 filesigma
  CHARACTER*15 filespectrum
!
!--------------------------
! main routine body
!--------------------------

  eta_large(1) = 0.d0
  eta_large(2) = 0.d0
  
  
  CALL grid_build(nw, wmax, wmin)

IF (nspin == 1) THEN

   ! Read PPA paramenter from file (named ppa.in) 
!   npar = 2 
   open (22,file="ppa.in",status="old",action="read")

   read(22,*) npar
   ALLOCATE(par(npar))
   ALLOCATE(par_tmp (npar))
   do ipar = 1, npar
     read(22,*) par(ipar)
   enddo
   close(22)
!   par(2) = sqrt(par(2))
!   par(4) = sqrt(par(4))
!   par(1) = -par(1)
!   par(3) = -par(3)
   !print *, par


!  !----------------------------------------------------------------------------------
  !print out the k-point grid 
  if(ionode)then
    open(44,file='kgrid.dat')
      write(44,'(A44)') '--------------------------------------------------------------------------------------'
      write(44,'(A4,A10,A10,A10,A10)') 'ik', 'xk(1,ik)','xk(2,ik)','xk(3,ik)','wk(ik)' 
      write(44,'(A44)') '--------------------------------------------------------------------------------------'
      do ik =1 ,nks  
        write(44,'(I4,F10.6,F10.6,F10.6,F10.6)') ik, xk(1,ik), xk(2,ik), xk(3,ik),wk(ik)
      enddo
    close(44)
  endif
!  !----------------------------------------------------------------------------------

  !construct the frequency grid ------------------------------------
  CALL grid_destroy()
  wmin1 = -40.0d0 ! eV
  wmax1 = 30.0d0  ! eV
  nw1 = 500
  CALL grid_build(nw1, wmax1, wmin1)
  !-----------------------------------------------------------------

 ! print*, wgrid 
  call print_epsilon (1, par, wgrid, nw1 )
      
  ALLOCATE( vpsi ( npwx ), STAT=ierr )
  IF (ierr/=0) CALL errore('gwppa','allocating vpsi', abs(ierr))
  ALLOCATE( psi  ( npwx ), STAT=ierr )
  IF (ierr/=0) CALL errore('gwppa','allocating psi', abs(ierr))
  ALLOCATE( sigmare (nw1), STAT=ierr )
  IF (ierr/=0) CALL errore('gwppa','allocating sigmare', abs(ierr))
  ALLOCATE( spectrum (nw1), STAT=ierr )
  IF (ierr/=0) CALL errore('gwppa','allocating spectrum', abs(ierr))

  nksigma     = 1 
  wksigma     = 2.0d0 / (nks - nksigma)
  !-----------------------------------------------------------------


! EXTERNAL LOOP OVER BANDS FOR WHICH THE SELF-ENERGY MUST BE CALCULATED
!  DO ibnd = 4,4 

!   IF(ionode) PRINT*, "Smallest band index: ",        ibndmin
!   IF(ionode) PRINT*, "tpiba2: ",      tpiba2,alat 
!   IF(ionode) PRINT*, "Biggest  band index: ",        ibndmax
!   IF(ionode) PRINT*, "Parameter for q-dependence: ", qmod_par 
  ALLOCATE( ieps_tmp(3,nw) ) 

  DO ibnd = ibndmin, ibndmax
     IF(ionode) PRINT*, "Evaluating sigma for band: ", ibnd 

     !calculate expectation value of v_xc ----------------------------
     ik=1
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     
     IF (ibnd .LT.10) THEN
       WRITE(label,'(A,I1)') "0", ibnd
     ELSE
       WRITE(label,'(I2)') ibnd 
     ENDIF 
    
     filespectrum = "spectrum_"//label//".dat"
     filesigma    = "sigma_"//label//".dat"
     
     tempevc (:) = (0.0d0,0.0d0)
     vpsi    (:) = (0.0d0,0.0d0)
     vxc         = (0.0d0,0.0d0) 
     vc1         = (0.0d0,0.0d0) 
     
!     !construct the exchange-correlation potential ------------------
     v%of_r(:,:) = (0.0d0, 0.0d0)
     tempevc(:)  = (0.0d0, 0.0d0)
     vpsi    (:) = (0.0d0, 0.0d0)
     CALL davcio (evc, nwordwfc, iunwfc, ik, -1)
     DO ig = 1 , npw
        tempevc(nls(igk(ig))) = evc(ig,ibnd)
     ENDDO
     
     CALL invfft ('Wave', tempevc(:), dffts)
     
     !evaluate matrix elements of v_xc
     !v%of_r(:,:) = (0.0d0, 0.0d0)
     CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
     DO ir = 1, dfftp%nnr
         tempevc (ir) = v%of_r (ir, 1) * tempevc (ir)
     ENDDO
     CALL fwfft ('Wave', tempevc(:), dffts)
     DO ig = 1, npw
        vpsi(ig) = tempevc (nls(igk(ig)))
     ENDDO
     vxc = ZDOTC (npw, evc (1, ibnd), 1, vpsi, 1)
     CALL mp_sum( vxc, intra_pool_comm )
!     if (ionode) PRINT * , 'VXC  = ',  real(vxc)*RYTOEV, 'VC  = ',  real(vc1)*RYTOEV
!     !--------------------------------------------

     !construct the exchange-correlation potential ------------------
     v%of_r(:,:) = (0.0d0, 0.0d0)
     tempevc(:)  = (0.0d0, 0.0d0)
     vpsi    (:) = (0.0d0, 0.0d0)
     CALL davcio (evc, nwordwfc, iunwfc, ik, -1)
     DO ig = 1 , npw
        tempevc(nls(igk(ig))) = evc(ig,ibnd)
     ENDDO

     CALL invfft ('Wave', tempevc(:), dffts)

     !evaluate matrix elements of v_xc
     !v%of_r(:,:) = (0.0d0, 0.0d0)
     CALL get_vc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
     DO ir = 1, dfftp%nnr
         tempevc (ir) = v%of_r (ir, 1) * tempevc (ir)
     ENDDO
     CALL fwfft ('Wave', tempevc(:), dffts)
     DO ig = 1, npw
        vpsi(ig) = tempevc (nls(igk(ig)))
     ENDDO
     vc1 = ZDOTC (npw, evc (1, ibnd), 1, vpsi, 1)
     !CALL mp_sum( vxc, intra_pool_comm )
     CALL mp_sum( vc1, intra_pool_comm )
     if (ionode) PRINT * , 'VXC  = ',  real(vxc)*RYTOEV, 'VC  = ',  real(vc1)*RYTOEV
     !--------------------------------------------

     
     IF ( ionode ) WRITE(6,*) 'EVHERE:', et(ibnd,1) * RYTOEV, efermi*RYTOEV 
     
     
     !evaluate the super-approximated gwppa self-energy ----------------------------
     sigmare (:) = (0.0d0,0.0d0)
     vxc         = (0.0d0,0.0d0)
     !vc1         = (0.0d0,0.0d0)
      
     !nshell = 1
     ng0vec = (2*nshell+1)**3
     IF(.NOT.ALLOCATED(gmap))THEN
       ALLOCATE (gmap(ngm,ng0vec))
     ENDIF
     IF(.NOT.ALLOCATED(g0vec))THEN
       ALLOCATE (g0vec(3,ng0vec))
     ENDIF
     gmap(:,:)=0.d0
     CALL refold (gmap,nshell,g0vec)


     !determine the minimum q-point, for the integration of the q=0 singularity 
     qmin = 100.d0
     DO ik = nksigma+1, nks
       DO ig = 1 , ng0vec 
         q =sum((xk(:,ik)-xk(:,1)+g0vec(:,ig))**2.d0)
         IF(q .LT. qmin)THEN
           qmin = q 
           indexqmin = ik 
           indexgmin = ig 
         ENDIF
       ENDDO
     ENDDO
  !-----------------------------------------------------------------


     CALL gk_sort (xk (1, 1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     CALL davcio (evc, nwordwfc, iunwfc, 1, -1)
     tempevc(:)=(0.d0,0.d0)
     tempevc(nls(igk(1:npw))) = evc(1:npw,ibnd)

!     !test gmap
!     DO ig=1,1
!       DO jg =1,(2*nshell+1)**3
!         IF(ionode.and..not. gmap(ig,jg).eq.0)THEN 
!           PRINT*, '--------------------------'
!           PRINT*, ig, jg ,igk(ig), gmap(ig,jg)
!           PRINT*, g(:,ig)
!           PRINT*, g(:,gmap(ig,jg)) 
!         ENDIF 
!       ENDDO
!     ENDDO

     !DO ik = 1, nks !always exclude the gamma (q=0) point from the summation 
     DO ik = nksigma+1, nks !always exclude the gamma (q=0) point from the summation 
       
!       IF (ionode) PRINT*, ' ik = ' , ik , ' / ', nks

       !construct the overlap --------------------------------------
       CALL gk_sort (xk (1,ik), ngm, g, ecutwfc / tpiba2, npw1, igk1, g2kin)
       CALL davcio (evc1, nwordwfc, iunwfc, ik, -1)
       
       DO ig = 1 , ng0vec !loop over G-vectors
     
         DO jbnd = 1, nbnd

!!!   !!!              !calculate an overlap by hand --------------------
!!!   !!!    !          IF(ik.EQ.1)THEN
!!!   !!!                ovlp = (0.d0,0.d0)
!!!   !!!                g0(:) = g0vec(:,ig)
!!!   !!!                gtmp(:,:) = g(:,:)
!!!   !!!                
!!!   !!!                call cryst_to_cart (ngm, gtmp, at, -1)
!!!   !!!                call cryst_to_cart (1, g0, at, -1)
!!!   !!!   
!!!   !!!                DO kg = 1, npw 
!!!   !!!                  g2(:) = gtmp(:,igk(kg))+g0(:)
!!!   !!!                  DO jg = 1, npw1 
!!!   !!!                    IF( (nint(gtmp(1,igk1(jg))) .EQ. nint(g2(1))) .AND. &
!!!   !!!                        (nint(gtmp(2,igk1(jg))) .EQ. nint(g2(2))) .AND. &
!!!   !!!                        (nint(gtmp(3,igk1(jg))) .EQ. nint(g2(3))) )THEN
!!!   !!!                      ovlp = ovlp + evc(kg,ibnd) * conjg(evc1(jg,jbnd)) 
!!!   !!!                    ENDIF
!!!   !!!                  ENDDO 
!!!   !!!                ENDDO
!!!   !!!   
!!!   !!!                IF(ionode)PRINT*, 'ibnd = ', ibnd, 'jbnd = ', jbnd
!!!   !!!                IF(ionode)PRINT*, 'ovlp1 = ', ovlp!,  ncount
!!!   !!!    !          ENDIF
!!!   !!!              !-------------------------------------------------

           ovlp   = (0.d0,0.d0)
           tempevc1(:)=(0.d0,0.d0)
           DO jg = 1, npw1 !loop over G'-vectors
               IF(gmap(igk1(jg),ig) .GT. 0 )THEN
                 tempevc1(nls(gmap(igk1(jg),ig))) = evc1(jg,jbnd)
               ENDIF
           ENDDO

           ovlp = ZDOTC (dffts%nnr, tempevc1, 1, tempevc, 1)
           CALL mp_sum( ovlp, intra_pool_comm )
           
       !------------------------------------------------------------
            qG = sum((xk(:,ik)-xk(:,1)+g0vec(:,ig))**2) *tpiba2 
            qG1 = sum((xk(:,ik) +g0vec(:,ig))**2)

            DO ipole = 1, npar/2 
               
              wtilde  = (par( ipole * 2 ) + (0.d0,1.d0) * eta_large(ipole) ) * dsqrt(1.d0 + qG / 1.2 )
              !if( ionode .and. ibnd == 1 .and. jbnd == 1 ) print *, qG, qG / tpiba2, wtilde
              !wtilde  = par( 2 + ( ipole - 1 ) * 2 ) * dsqrt(1.d0 + qG1 / qmod_par)

!              if( ionode .and. ibnd == 1 .and. jbnd == 1 ) print * , tpiba2
!              if( ionode .and. ibnd == 1 .and. jbnd == 1 )  print * , qG1, dsqrt(1.d0 + qG1 /0.6), ik
              wplasma = par( ipole * 2 - 1 )
 
              prefactor = ovlp * conjg (ovlp) * &
                 wplasma**2 / (2.0d0*wtilde ) * 4.0d0 * PI 
      
              IF (.NOT. (ik.EQ.indexqmin .AND. ig .EQ. indexgmin)) THEN
                 DO iw = 1, nw1
                    sigmare(iw) = sigmare(iw) + &
                     prefactor / (wgrid(iw) + SIGN(real(wtilde),efermi - et(jbnd,ik)) - et(jbnd,ik)  &
                     * RYTOEV  - (0.0d0,1.0d0) * SIGN(eta,efermi - et(jbnd,ik))) &
                     / qG * wksigma / omega *RYTOEV !* 1.d0/ (1.d0 + qG /kTF/tpiba2)
   
                 ENDDO
              ELSE 
                 !IF (ionode) WRITE(6,*) "I am neglecting the q=0 singularity!"
                 DO iw = 1, nw1
  
                    sigmare(iw) = sigmare(iw) +  &
                     prefactor / ( wgrid(iw) + SIGN(real(wtilde),efermi - et(jbnd,ik)) - et(jbnd,ik)  &
                     * RYTOEV  - (0.0d0,1.0d0) * SIGN(eta,efermi - et(jbnd,ik))) * &
                     ( 3.d0 * wksigma / 4.d0 / PI / omega )**(1.d0/3.d0)/PI * RYTOEV
  
                 ENDDO !loop over frequency points
              ENDIF ! (.NOT.(ik.EQ.indexqmin .AND. absg0vec .LT. 1e-3)) 
            ENDDO
          ENDDO !loop over all bands
       ENDDO !loop over G-vectors shells 
     ENDDO !loop over k-points
     CALL mp_sum( sigmare, inter_pool_comm )
     !------------------------------------------------------------------------------
     
     
     !evaluate the spectral function ----------------------------------------------- 
     DO iw = 1, nw1 
        spectrum(iw) =  abs(aimag (sigmare(iw))) /       &
              ((wgrid(iw) - (et(ibnd,1)-real(vc1)) * RYTOEV &
              - real(sigmare(iw)))**2 + aimag (sigmare(iw))**2)
     ENDDO 
     !------------------------------------------------------------------------------ 
     IF( ionode ) PRINT*, ' v_c [eV]: ', real(vc1)  * RYTOEV  
     
     
     !print to file ---------------------------------------------------------------- 
     IF ( ionode ) THEN
        OPEN (30, FILE=filesigma)
        OPEN (31, FILE=filespectrum)
           DO iw = 1, nw1
             WRITE(30,*) wgrid(iw), real(sigmare(iw)), aimag (sigmare(iw))
             WRITE(31,*) wgrid(iw)-efermi*RYTOEV, spectrum(iw)
             !WRITE(31,*) wgrid(iw)-efermi*RYTOEV, spectrum(iw)
           ENDDO
        CLOSE(30)
        CLOSE(31)
     ENDIF
!     DEALLOCATE (sigmare)
  !------------------------------------------------------------------------------ 

  ENDDO ! LOOP OVER OCCUPIED BANDS

ELSEIF (nspin == 2 ) THEN
  IF ( ionode ) WRITE(6,*) 'GWPPA works only for nspin = 1 !'
ENDIF 

IF(ionode ) WRITE(6,*) " DONE !!! "

  !--------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------
  !
  ! local cleaning
  !
  CALL grid_destroy()
  !
  DEALLOCATE (  dipole, dipole_aux )


END SUBROUTINE gwppa_calc

!-----------------------------------------------------------------------------
SUBROUTINE eps_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax, shift, &
                      metalcalc , nspin)
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss
  USE io_global,            ONLY : ionode, stdout
  !
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(in) :: metalcalc
  !
  ! local variables
  !
  INTEGER       :: i, ik, iband1, iband2,is
  INTEGER       :: iw, iwp, ierr
  REAL(DP)      :: etrans, const, w, renorm(3)
  !
  REAL(DP), ALLOCATABLE    :: epsr(:,:), epsi(:,:), epsrc(:,:,:), epsic(:,:,:)
  REAL(DP), ALLOCATABLE    :: ieps(:,:), eels(:,:), iepsc(:,:,:), eelsc(:,:,:)
  REAL(DP), ALLOCATABLE    :: dipole(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_aux(:,:,:)
!
!--------------------------
! main routine body
!--------------------------
!
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin)
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( dipole(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole', abs(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', abs(ierr) )

!
! spin unresolved calculation
!
IF (nspin == 1) THEN
  !
  ALLOCATE( epsr( 3, nw), epsi( 3, nw), eels( 3, nw), ieps(3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', abs(ierr))

  !
  ! initialize response functions
  !
  epsr(:,:)  = 0.0_DP
  epsi(:,:)  = 0.0_DP
  ieps(:,:)  = 0.0_DP

  !
  ! main kpt loop
  !
  kpt_loop: &
  DO ik = 1, nks
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g
     !                           recover g parallelism getting the total dipole matrix
     !
     CALL dipole_calc( ik, dipole_aux, metalcalc , nbndmin, nbndmax)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * conjg(dipole_aux(:,:,:)), DP )

     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     !Interband
     !
     DO iband2 = nbndmin,nbndmax
         !
         IF ( focc(iband2,ik) < 2.0d0) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF (iband1==iband2) CYCLE
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
         IF (abs(focc(iband2,ik)-focc(iband1,ik))< 1e-3) CYCLE
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epsi(:,iw) = epsi(:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )

                   epsr(:,iw) = epsr(:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )
               ENDDO

         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF ( focc(iband1,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                  epsi(:,iw) = epsi(:,iw) +  dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**2 * (exp((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )

                  epsr(:,iw) = epsr(:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                                            (exp((et(iband1,ik)-efermi)/degauss )) * w**2 / &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loop

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epsr, inter_pool_comm )
  CALL mp_sum( epsi, inter_pool_comm )

  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epsr(:,:) = 1.0_DP + epsr(:,:) * const
  epsi(:,:) =          epsi(:,:) * const

  !
  ! Calculation of eels spectrum
  !
  DO iw = 1, nw
      !
      eels(:,iw) = epsi(:,iw) / ( epsr(:,iw)**2 + epsi(:,iw)**2 )
      !
  ENDDO

  !
  !  calculation of dielectric function on the immaginary frequency axe
  !

  DO iw = 1, nw
  DO iwp = 2, nw
      !
      ieps(:,iw) = ieps(:,iw) + wgrid(iwp) * epsi(:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2)
      !
  ENDDO
  ENDDO

  ieps(:,:) = 1.0d0 + 2 / PI * ieps(:,:) * alpha

  !
  ! check  dielectric function  normalizzation via sumrule
  !
 DO i=1,3
     renorm(i) = alpha * sum( epsi(i,:) * wgrid(:) )
 ENDDO
  !
  IF ( ionode ) THEN
      !
      WRITE(stdout,"(/,5x, 'The bulk xx plasmon frequency [eV] is: ',f15.9 )")  sqrt(renorm(1) * 2.0d0 / PI)
      WRITE(stdout,"(5x, 'The bulk yy plasmon frequency [eV] is: ',f15.9 )")  sqrt(renorm(2) * 2.0d0 / PI)
      WRITE(stdout,"(5x, 'The bulk zz plasmon frequency [eV] is: ',f15.9 )")  sqrt(renorm(3) * 2.0d0 / PI)
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      ! write results on data files
      !

      OPEN (30, FILE='epsr.dat', FORM='FORMATTED' )
      OPEN (40, FILE='epsi.dat', FORM='FORMATTED' )
      OPEN (41, FILE='eels.dat', FORM='FORMATTED' )
      OPEN (42, FILE='ieps.dat', FORM='FORMATTED' )
      !
      WRITE(30, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(40, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(41, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(42, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      !
      DO iw =1, nw
          !
          WRITE(30,"(4f15.6)") wgrid(iw), epsr(1:3, iw)
          WRITE(40,"(4f15.6)") wgrid(iw), epsi(1:3, iw)
          WRITE(41,"(4f15.6)") wgrid(iw), eels(1:3, iw)
          WRITE(42,"(4f15.6)") wgrid(iw), ieps(1:3, iw)
          !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      !
  ENDIF

  DEALLOCATE ( epsr, epsi, eels, ieps)
!
! collinear spin calculation
!
ELSEIF (nspin == 2 ) THEN
  !
  ALLOCATE( epsrc( 0:1, 3, nw), epsic( 0:1,3, nw), eelsc( 0:1,3, nw), iepsc(0:1,3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', abs(ierr))

  !
  ! initialize response functions
  !
  epsrc(:,:,:)  = 0.0_DP
  epsic(:,:,:)  = 0.0_DP
  iepsc(:,:,:)  = 0.0_DP

  !
  ! main kpt loop
  !

spin_loop: &
DO is=0,1
  kpt_loopspin: &
! if nspin=2 the number of nks must be even (even if the calculation
! is performed at gamma point only), so nks must be always a multiple of 2
  DO ik = 1 + is * int(nks/2), int(nks/2) +  is * int(nks/2)
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g
     !                           recover g parallelism getting the total dipole matrix
     !
     CALL dipole_calc( ik, dipole_aux, metalcalc , nbndmin, nbndmax)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * conjg(dipole_aux(:,:,:)), DP )

     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     !Interband
     !
     DO iband2 = nbndmin,nbndmax
         !
         IF ( focc(iband2,ik) < 1.0d0) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF (iband1==iband2) CYCLE
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
         IF (abs(focc(iband2,ik)-focc(iband1,ik))< 1e-3) CYCLE
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epsic(is,:,iw) = epsic(is,:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )

                   epsrc(is,:,iw) = epsrc(is,:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )
               ENDDO

         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF ( focc(iband1,ik) < 1.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                  epsic(is,:,iw) = epsic(is,:,iw) +  dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**2 * (exp((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )

                  epsrc(is,:,iw) = epsrc(is,:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                                            (exp((et(iband1,ik)-efermi)/degauss )) * w**2 / &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loopspin
ENDDO spin_loop
  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epsr, inter_pool_comm )
  CALL mp_sum( epsi, inter_pool_comm )

  !
  ! impose the correct normalization
  !
  const = 128.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epsrc(:,:,:) = 1.0_DP + epsrc(:,:,:) * const
  epsic(:,:,:) =          epsic(:,:,:) * const

  !
  ! Calculation of eels spectrum
  !
  DO iw = 1, nw
      !
      eelsc(:,:,iw) = epsic(:,:,iw) / ( epsrc(:,:,iw)**2 + epsic(:,:,iw)**2 )
      !
  ENDDO

  !
  !  calculation of dielectric function on the immaginary frequency axe
  !

  DO iw = 1, nw
  DO iwp = 2, nw
      !
      iepsc(:,:,iw) = iepsc(:,:,iw) + wgrid(iwp) * epsic(:,:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2)
      !
  ENDDO
  ENDDO

  iepsc(:,:,:) = 1.0d0 + 2.0_DP / PI * iepsc(:,:,:) * alpha

  IF (ionode) THEN
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      ! write results on data files
      !

      OPEN (30, FILE='uepsr.dat', FORM='FORMATTED' )
      OPEN (40, FILE='uepsi.dat', FORM='FORMATTED' )
      OPEN (41, FILE='ueels.dat', FORM='FORMATTED' )
      OPEN (42, FILE='uieps.dat', FORM='FORMATTED' )
      OPEN (43, FILE='depsr.dat', FORM='FORMATTED' )
      OPEN (44, FILE='depsi.dat', FORM='FORMATTED' )
      OPEN (45, FILE='deels.dat', FORM='FORMATTED' )
      OPEN (46, FILE='dieps.dat', FORM='FORMATTED' )
      OPEN (47, FILE='epsr.dat', FORM='FORMATTED' )
      OPEN (48, FILE='epsi.dat', FORM='FORMATTED' )
      OPEN (49, FILE='eels.dat', FORM='FORMATTED' )
      OPEN (50, FILE='ieps.dat', FORM='FORMATTED' )
      !
      WRITE(30, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(40, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(41, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(42, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      WRITE(43, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(44, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(45, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(46, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      WRITE(47, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(48, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(49, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(50, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      !
      DO iw =1, nw
          !
          WRITE(30,"(4f15.6)") wgrid(iw), epsrc(0,1:3, iw)
          WRITE(40,"(4f15.6)") wgrid(iw), epsic(0,1:3, iw)
          WRITE(41,"(4f15.6)") wgrid(iw), eelsc(0,1:3, iw)
          WRITE(42,"(4f15.6)") wgrid(iw), iepsc(0,1:3, iw)
          WRITE(43,"(4f15.6)") wgrid(iw), epsrc(1,1:3, iw)
          WRITE(44,"(4f15.6)") wgrid(iw), epsic(1,1:3, iw)
          WRITE(45,"(4f15.6)") wgrid(iw), eelsc(1,1:3, iw)
          WRITE(46,"(4f15.6)") wgrid(iw), iepsc(1,1:3, iw)
          WRITE(47,"(4f15.6)") wgrid(iw), epsrc(1,1:3, iw)+epsrc(0,1:3, iw)
          WRITE(48,"(4f15.6)") wgrid(iw), epsic(1,1:3, iw)+epsic(0,1:3, iw)
          WRITE(49,"(4f15.6)") wgrid(iw), eelsc(1,1:3, iw)+eelsc(0,1:3, iw)
          WRITE(50,"(4f15.6)") wgrid(iw), iepsc(1,1:3, iw)+iepsc(0,1:3, iw)
          !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
      CLOSE(46)
      CLOSE(47)
      CLOSE(48)
      CLOSE(49)
      CLOSE(50)
      !
  ENDIF
  DEALLOCATE ( epsrc, epsic, eelsc, iepsc)
ENDIF
  !
  ! local cleaning
  !
  CALL grid_destroy()
  !
  DEALLOCATE (  dipole, dipole_aux )

END SUBROUTINE eps_calc

!----------------------------------------------------------------------------------------
SUBROUTINE jdos_calc ( smeartype,intersmear,nw,wmax,wmin,nbndmin,nbndmax,shift,nspin )
  !--------------------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE wvfct,                ONLY : nbnd, et
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy
  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear, shift
  CHARACTER(*),    INTENT(in) :: smeartype
  !
  ! local variables
  !
  INTEGER       :: ik, is, iband1, iband2
  INTEGER       :: iw, ierr
  REAL(DP)      :: etrans, w, renorm, count, srcount(0:1), renormzero,renormuno
  !
  REAL(DP), ALLOCATABLE    :: jdos(:),srjdos(:,:)
!
!--------------------------
! main routine body
!--------------------------
!
! No wavefunctions are needed in order to compute jdos, only eigenvalues,
! they are distributed to each task so
! no mpi calls are necessary in this routine
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin )

!
! spin unresolved calculation
!
IF (nspin == 1) THEN
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( jdos(nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating jdos',abs(ierr))
  !
  ! initialize jdos
  !
  jdos(:)=0.0_DP

  ! Initialising a counter for the number of transition
  count=0.0_DP

  !
  ! main kpt loop
  !

  IF (smeartype=='lorentz') THEN

    kpt_lor: &
    DO ik = 1, nks
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
           IF ( focc(iband2,ik) <  2.0d0) THEN
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 count = count + (focc(iband1,ik)-focc(iband2,ik))
                 !
                 ! loop over frequencies
                 !
                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
       ENDDO
           ENDIF
       ENDDO

    ENDDO kpt_lor

  ELSEIF (smeartype=='gauss') THEN

    kpt_gauss: &
    DO ik = 1, nks

       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !

                 count=count+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                                exp(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * sqrt(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO kpt_gauss

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//trim(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !

  jdos(:)=jdos(:)/count

  !
  ! check jdos normalization
  !

  renorm = alpha * sum( jdos(:) )
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Integration over JDOS gives: ',f15.9,' instead of 1.0d0' )") renorm
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")

     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     JDOS [1/eV] ')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), jdos(iw)
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF
  !
  ! local cleaning
  !
  DEALLOCATE ( jdos )

!
! collinear spin calculation
!
ELSEIF(nspin==2) THEN
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( srjdos(0:1,nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating spin resolved jdos',abs(ierr))
  !
  ! initialize jdos
  !
  srjdos(:,:)=0.0_DP

  ! Initialising a counter for the number of transition
  srcount(:)=0.0_DP

  !
  ! main kpt loop
  !

  IF (smeartype=='lorentz') THEN

  DO is=0,1
    ! if nspin=2 the number of nks must be even (even if the calculation
    ! is performed at gamma point only), so nks must be always a multiple of 2
    DO ik = 1 + is * int(nks/2), int(nks/2) +  is * int(nks/2)
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
           IF ( focc(iband2,ik) <  2.0d0) THEN
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !
                 srcount(is)=srcount(is)+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     srjdos(is,iw) = srjdos(is,iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
       ENDDO
           ENDIF
       ENDDO

    ENDDO
 ENDDO

  ELSEIF (smeartype=='gauss') THEN

  DO is=0,1
    ! if nspin=2 the number of nks must be even (even if the calculation
    ! is performed at gamma point only), so nks must be always a multiple of 2
    DO ik = 1 + is * int(nks/2), int(nks/2) +  is * int(nks/2)
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !

                 srcount(is)=srcount(is)+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     srjdos(is,iw) = srjdos(is,iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                                exp(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * sqrt(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO
 ENDDO

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//trim(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !
  DO is = 0,1
    srjdos(is,:)=srjdos(is,:)/srcount(is)
  ENDDO
  !
  ! check jdos normalization
  !

  renormzero = alpha * sum( srjdos(0,:) )
  renormuno = alpha * sum( srjdos(1,:) )
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Integration over spin UP JDOS gives: ',f15.9,' instead of 1.0d0' )") renormzero
     WRITE(stdout,"(/,5x, 'Integration over spin DOWN JDOS gives: ',f15.9,' instead of 1.0d0' )") renormuno
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")

     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     UJDOS [1/eV]      DJDOS[1:eV]')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), srjdos(0,iw), srjdos(1,iw)
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF

  DEALLOCATE ( srjdos )
ENDIF
  !
  ! local cleaning
  !
  CALL grid_destroy()

END SUBROUTINE jdos_calc

!-----------------------------------------------------------------------------
SUBROUTINE offdiag_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax,&
                          shift, metalcalc, nspin )
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss
  USE grid_module,          ONLY : focc, wgrid, grid_build, grid_destroy
  USE io_global,            ONLY : ionode, stdout
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum

  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(in) :: metalcalc
  !
  ! local variables
  !
  INTEGER       :: ik, iband1, iband2
  INTEGER       :: iw, ierr, it1, it2
  REAL(DP)      :: etrans, const, w
  !
  COMPLEX(DP), ALLOCATABLE :: dipole_aux(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: epstot(:,:,:),dipoletot(:,:,:,:)
  !
  !--------------------------
  ! main routine body
  !--------------------------
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin )
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( dipoletot(3,3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipoletot', abs(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', abs(ierr) )
  !
  ALLOCATE(epstot( 3,3, nw),STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating epstot', abs(ierr))

   !
   ! initialize response functions
   !
   epstot  = (0.0_DP,0.0_DP)
   !
   ! main kpt loop
   !
   DO ik = 1, nks
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g
     !                           recover g parallelism getting the total dipole matrix
     !
     CALL dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax)
     !
     DO it2 = 1, 3
        DO it1 = 1, 3
           dipoletot(it1,it2,:,:) = tpiba2 * dipole_aux(it1,:,:) * conjg( dipole_aux(it2,:,:) )
        ENDDO
     ENDDO
     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     DO iband2 = 1,nbnd
         IF ( focc(iband2,ik) <  2.0d0) THEN
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
             !
             ! transition energy
             !
             etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
             !
             IF (abs(focc(iband2,ik)-focc(iband1,ik))< 1e-4) CYCLE
             !
             ! loop over frequencies
             !
             DO iw = 1, nw
                  !
                  w = wgrid(iw)
                  !
                  epstot(:,:,iw) = epstot(:,:,iw) + dipoletot(:,:,iband1,iband2)*RYTOEV**3/(etrans) *&
                                   focc(iband1,ik)/(etrans**2 - w**2 - (0,1)*intersmear*w)
             ENDDO
             !
         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband1,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epstot(:,:,iw) = epstot(:,:,iw) - dipoletot(:,:,iband1,iband1)* &
                                RYTOEV**2 * (exp((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**2 + (0,1)*intrasmear*w)*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epstot, inter_pool_comm )
  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epstot(:,:,:) = epstot(:,:,:) * const
  !
  ! add diagonal term
  !
  epstot(1,1,:) = 1.0_DP + epstot(1,1,:)
  epstot(2,2,:) = 1.0_DP + epstot(2,2,:)
  epstot(3,3,:) = 1.0_DP + epstot(3,3,:)
  !
  ! write results on data files
  !
  IF (ionode) THEN
      !
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      OPEN (41, FILE='epsxx.dat', FORM='FORMATTED' )
      OPEN (42, FILE='epsxy.dat', FORM='FORMATTED' )
      OPEN (43, FILE='epsxz.dat', FORM='FORMATTED' )
      OPEN (44, FILE='epsyx.dat', FORM='FORMATTED' )
      OPEN (45, FILE='epsyy.dat', FORM='FORMATTED' )
      OPEN (46, FILE='epsyz.dat', FORM='FORMATTED' )
      OPEN (47, FILE='epszx.dat', FORM='FORMATTED' )
      OPEN (48, FILE='epszy.dat', FORM='FORMATTED' )
      OPEN (49, FILE='epszz.dat', FORM='FORMATTED' )
      !
      WRITE(41, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(42, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(43, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(44, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(45, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(46, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(47, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(48, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(49, "(2x,'# energy grid [eV]     epsr     epsi')" )
      !
      DO iw =1, nw
         !
         WRITE(41,"(4f15.6)") wgrid(iw), REAL(epstot(1,1, iw)), aimag(epstot(1,1, iw))
         WRITE(42,"(4f15.6)") wgrid(iw), REAL(epstot(1,2, iw)), aimag(epstot(1,2, iw))
         WRITE(43,"(4f15.6)") wgrid(iw), REAL(epstot(1,3, iw)), aimag(epstot(1,3, iw))
         WRITE(44,"(4f15.6)") wgrid(iw), REAL(epstot(2,1, iw)), aimag(epstot(2,1, iw))
         WRITE(45,"(4f15.6)") wgrid(iw), REAL(epstot(2,2, iw)), aimag(epstot(2,2, iw))
         WRITE(46,"(4f15.6)") wgrid(iw), REAL(epstot(2,3, iw)), aimag(epstot(2,3, iw))
         WRITE(47,"(4f15.6)") wgrid(iw), REAL(epstot(3,1, iw)), aimag(epstot(3,1, iw))
         WRITE(48,"(4f15.6)") wgrid(iw), REAL(epstot(3,2, iw)), aimag(epstot(3,2, iw))
         WRITE(49,"(4f15.6)") wgrid(iw), REAL(epstot(3,3, iw)), aimag(epstot(3,3, iw))
         !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      !
  ENDIF

  !
  ! local cleaning
  !
  CALL grid_destroy()
  DEALLOCATE ( dipoletot, dipole_aux, epstot )

END SUBROUTINE offdiag_calc

!-------------------------------------------------
SUBROUTINE occ_calc ()
  !-------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE klist,     ONLY : nkstot, wk, degauss
  USE wvfct,     ONLY : nbnd, wg, et
  USE ener,      ONLY : ef
  USE mp_global, ONLY : me_pool
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE  :: focc(:,:),foccp(:,:)
  CHARACTER(25)          :: filename
  INTEGER                :: ierr, i, ik
  !
  ALLOCATE ( focc( nbnd, nkstot), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', abs(ierr))
  !
  ALLOCATE ( foccp( nbnd, nkstot), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating foccp', abs(ierr))

  IF (me_pool==0) THEN
      !
      filename = 'occupations.dat'
!      WRITE(filename,"(I3,'.occupation.dat')")me_pool
      OPEN (unit=50, file=trim(filename))
      WRITE(50,*) '#energy (Ry)      occupation factor       derivative'

      DO ik = 1,nkstot
      DO i  = 1,nbnd
           focc(i,ik)= wg(i, ik ) * 2.0_DP/wk( ik )
           foccp(i,ik)= 2* exp((et(i,ik)-ef)/degauss)/((1+exp((et(i,ik)-ef)/degauss))**2*degauss)
           WRITE(50,*)et(i,ik),focc(i,ik),foccp(i,ik)
      ENDDO
      ENDDO

      CLOSE (50)
      !
  ENDIF
  !
  DEALLOCATE ( focc, STAT=ierr)
  CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))
  !
  DEALLOCATE ( foccp, STAT=ierr)
  CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))

END SUBROUTINE occ_calc

!--------------------------------------------------------------------
SUBROUTINE dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax )
  !------------------------------------------------------------------
  !
  ! Evaluate overlaps (n, k-q | n', k) where k is fixed 
  ! (it should be the first k-point in the pwscf input file),
  ! and q (labeled by ik) varies in the BZ
  !
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npw, nbnd, igk, g2kin, ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE grid_module,          ONLY : focc
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(in)        :: ik,nbndmin,nbndmax
  COMPLEX(DP), INTENT(inout) :: dipole_aux(3,nbnd,nbnd)
  LOGICAL, INTENT(in)        :: metalcalc
  !
  ! local variables
  INTEGER :: iband1,iband2,ig
  COMPLEX(DP)   :: caux

  !
  ! Routine Body
  !
  CALL start_clock( 'dipole_calc' )

  !
  ! setup k+G grids for each kpt
  !
  CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  !
  ! read wfc for the given kpt
  !
  CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  DO iband2 = nbndmin,nbndmax
      IF ( focc(iband2,ik) <  2.0d0) THEN
  DO iband1 = nbndmin,nbndmax
      !
      IF ( iband1==iband2 ) CYCLE
      IF ( focc(iband1,ik) >= 1e-4 ) THEN
            !
            DO  ig=1,npw
                 !
                 caux= conjg(evc(ig,iband1))*evc(ig,iband2)
                 !
                 dipole_aux(:,iband1,iband2) = dipole_aux(:,iband1,iband2) + &
                       ( g(:,igk(ig)) ) * caux
                 !
            ENDDO
      ENDIF
      !
  ENDDO
      ENDIF
  ENDDO
  !
  ! The diagonal terms are taken into account only if the system is treated like a metal, not
  ! in the intraband therm. Because of this we can recalculate the diagonal component of the dipole
  ! tensor directly as we need it for the intraband therm, without interference with interband one.
  !
  IF (metalcalc) THEN
     !
     DO iband1 = nbndmin,nbndmax
        DO  ig=1,npw
          !
          caux= conjg(evc(ig,iband1))*evc(ig,iband1)
          !
          dipole_aux(:,iband1,iband1) = dipole_aux(:,iband1,iband1) + &
                                        ( g(:,igk(ig))+ xk(:,ik) ) * caux
          !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! recover over G parallelization (intra_pool)
  !
  CALL mp_sum( dipole_aux, intra_pool_comm )
  !
  CALL stop_clock( 'dipole_calc' )
  !
END SUBROUTINE dipole_calc


  !
  !----------------------------------------------------------------
  SUBROUTINE refold ( gmap, nshell, g0vec )
  !----------------------------------------------------------------
  !
  ! the map is defined as follows
  ! g(:,gmap(ig,ig0)) = g(:,ig) - g0vec(:,ig0)
  !
  ! at the exit, the folding vectors are in cartesian coordinates
  !
  !----------------------------------------------------------------
  !
  USE parameters
  USE gvect !,     ONLY: g, ngm
  USE cell_base, ONLY: bg, at 
  USE io_global, ONLY: ionode
  implicit none
  !
  integer   ::  notfound, g2(3), count
  integer   :: ig0, ig, igp, ig1, ig2, ig3
  integer   :: nshell, ng0vec, npw
  real(DP)  :: g1(3) 
  real(DP)  :: g0vec(3,(nshell*2+1)**3)
  integer   :: gmap(ngm,(nshell*2+1)**3)
  real(DP)  :: gtmp(3,ngm)
  
  !
  ! the 3^3 possible translations (in very anisotropic materials
  ! this could go to 5^3 or even more - see EPW code)
  !
  ! crystal coordinates
  !
  ng0vec = 0
  do ig1 = -nshell,nshell
    do ig2 = -nshell,nshell
      do ig3 = -nshell,nshell
        ng0vec = ng0vec + 1
        !print*, ng0vec
        g0vec(1,ng0vec) = ig1
        g0vec(2,ng0vec) = ig2
        g0vec(3,ng0vec) = ig3
      enddo
    enddo
  enddo
  !
  ! bring all the G-vectors in crystal coordinates
  ! (in real codes we use the array of the miller indices mill_ )
  !  
!  gtmp(:,:) = g(:,:)
!  call cryst_to_cart (ngm, gtmp, at, -1)
  call cryst_to_cart (ngm, g, at, -1)
  ! 
  count = 0
  notfound = 0
  do ig = 1, ngm
    !
    do ig0 = 1, ng0vec
      !
      gmap (ig,ig0) = 0 
      count = count + 1
      !
      g2 = nint ( g(:,ig) - g0vec(:,ig0) )
      !
      do igp = 1, ngm

!        if ( ( g2(1) .eq. nint(gtmp(1,igp)) ) .and. &
!             ( g2(2) .eq. nint(gtmp(2,igp)) ) .and. &
!             ( g2(3) .eq. nint(gtmp(3,igp)) ) )     &
        if ( ( g2(1) .eq. nint(g(1,igp)) ) .and. &
             ( g2(2) .eq. nint(g(2,igp)) ) .and. &
             ( g2(3) .eq. nint(g(3,igp)) ) )     &
              gmap (ig,ig0) = igp
!              if(ionode)print*,' '
!              if(ionode)print*,'  --- ', ig, ig0 ,igp, gmap(ig,ig0),' ---  '
!              if(ionode)print*,g2(:)
!              if(ionode)print*,g(:,igp)
!              if(ionode)print*,g(:,gmap(ig,ig0))
      enddo
      if (gmap (ig,ig0).eq.0) notfound = notfound + 1
      !
    enddo
    !
  enddo
  write(6,'(4x,"refold: notfound = ",i6," out of ",i6)') notfound,count
  !if(ionode) print*, gmap 
  !
  ! back to cartesian
  !  
  
  call cryst_to_cart (ngm, g, bg, 1)
  call cryst_to_cart (ng0vec, g0vec, bg, 1)
  !
  END SUBROUTINE refold
  !
  !
  !----------------------------------------------------------------
  SUBROUTINE refold2 ( gmap, nshell, g0vec )
  !----------------------------------------------------------------
  !
  ! the map is defined as follows
  ! g(:,gmap(ig,ig0)) = g(:,ig) - g0vec(:,ig0)
  !
  ! at the exit, the folding vectors are in cartesian coordinates
  !
  !----------------------------------------------------------------
  !
  USE parameters
  USE gvect !,     ONLY: g, ngm
  USE cell_base, ONLY: bg, at 
  USE io_global, ONLY: ionode
  implicit none
  !
  integer   ::  notfound, g2(3), count
  integer   :: ig0, ig, igp, ig1, ig2, ig3
  integer   :: nshell, ng0vec, npw
!  parameter :: ng0vec
  real(DP)  :: g1(3) 
  real(DP)  :: g0vec(3,(nshell*2+1)**3)
  integer   :: gmap(ngm,(nshell*2+1)**3)
  
  !
  ! the 3^3 possible translations (in very anisotropic materials
  ! this could go to 5^3 or even more - see EPW code)
  !
  ! crystal coordinates
  !
  ng0vec = 0
  do ig1 = -nshell,nshell
    do ig2 = -nshell,nshell
      do ig3 = -nshell,nshell
        ng0vec = ng0vec + 1
        !print*, ng0vec
        g0vec(1,ng0vec) = ig1
        g0vec(2,ng0vec) = ig2
        g0vec(3,ng0vec) = ig3
      enddo
    enddo
  enddo
  !
  ! bring all the G-vectors in crystal coordinates
  ! (in real codes we use the array of the miller indices mill_ )
  !  
  ! call cryst_to_cart (ngm, g, at, -1)
  ! 
  count = 0
  notfound = 0
  do ig = 1, ngm
!    if(ionode)print*,g(:,ig)
    !
    do ig0 = 1, ng0vec
      !
      gmap (ig,ig0) = 0 
      count = count + 1
      !
      g2 = nint( g(:,ig) - g0vec(:,ig0) )
      !
      do igp = 1, ngm
        if ( ( g2(1) .eq. nint(g(1,igp)) ) .and. &
             ( g2(2) .eq. nint(g(2,igp)) ) .and. &
             ( g2(3) .eq. nint(g(3,igp)) ) )  then 
              gmap (ig,ig0) = igp
!              if(ionode)print*,' '
!              if(ionode)print*,'  --- ', ig, ig0 ,igp, gmap(ig,ig0),' ---  '
!              if(ionode)print*,g2(:)
!              if(ionode)print*,g(:,igp)
!              if(ionode)print*,g(:,gmap(ig,ig0))
        endif 
              
      enddo
      if (gmap (ig,ig0).eq.0) notfound = notfound + 1
      !
    enddo
    !
  enddo
  write(6,'(4x,"refold: notfound = ",i6," out of ",i6)') notfound,count
  !if(ionode) print*, gmap 
  !
  ! back to cartesian
  !  
!  call cryst_to_cart (ngm, g, bg, 1)
!  call cryst_to_cart (ng0vec, g0vec, bg, 1)
  !
  END SUBROUTINE refold2
  !
SUBROUTINE get_vc( rho, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !
  ! ... Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout, ionode
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin
  USE scf,              ONLY : scf_type
  USE mp_global,        ONLY : intra_pool_comm, intra_bgrp_comm, mpime
  USE mp,               ONLY : mp_sum


  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin), vtxc, etxc
    ! V_xc potential
    ! integral V_xc * rho
    ! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhox, arhox, zeta, amag, vs, ex, ec, vx(2), vc(2), rhoneg(2)
    ! the total charge in each point
    ! the absolute value of the charge
    ! the absolute value of the charge
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, ipol
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20
  !
  !
  CALL start_clock( 'v_xc' )
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     !
     ! ... spin-unpolarized case
     !
!$omp parallel do private( rhox, arhox, ex, ec, vx, vc ), &
!$omp             reduction(+:etxc,vtxc), reduction(-:rhoneg)
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           CALL xc( arhox, ex, ec, vx(1), vc(1) )
           !
           v(ir,1) = e2*( vc(1)  )
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           !
           vtxc = vtxc + v(ir,1) * rho%of_r(ir,1)
           !
        ENDIF
        !
        IF ( rho%of_r(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        !
     END DO
!$omp end parallel do
     !
  ELSE   IF ( nspin == 2 ) THEN
!     IF(ionode) THEN 
!       PRINT*, 'This is supposed to work only with ispin = 1 '   
!     ENDIF 
!     stop
!  ENDIF
     !
     ! ... spin-polarized case
     !
!$omp parallel do private( rhox, arhox, zeta, ex, ec, vx, vc ), &
!$omp             reduction(+:etxc,vtxc), reduction(-:rhoneg)
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho%of_r(ir,2) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = ( rho%of_r(ir,1) - rho%of_r(ir,2) ) / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) zeta = SIGN( 1.D0, zeta )
           !
           IF ( rho%of_r(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
           IF ( rho%of_r(ir,2) < 0.D0 ) rhoneg(2) = rhoneg(2) - rho%of_r(ir,2)
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           v(ir,:) = e2*( vx(:) + vc(:) )
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           !
           vtxc = vtxc + v(ir,1) * rho%of_r(ir,1) + v(ir,2) * rho%of_r(ir,2)
           !
        END IF
        !
     END DO
!$omp end parallel do
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     ! ... noncolinear case
     !
     DO ir = 1,dfftp%nnr
        !
        amag = SQRT( rho%of_r(ir,2)**2 + rho%of_r(ir,3)**2 + rho%of_r(ir,4)**2 )
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        !
        IF ( rho%of_r(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = amag / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) THEN
              !
              rhoneg(2) = rhoneg(2) + 1.D0 / omega
              !
              zeta = SIGN( 1.D0, zeta )
              !
           END IF
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           vs = 0.5D0*( vx(1) + vc(1) - vx(2) - vc(2) )
           !
           v(ir,1) = e2*( 0.5D0*( vx(1) + vc(1) + vx(2) + vc(2 ) ) )
           !
           IF ( amag > vanishing_mag ) THEN
              !
              DO ipol = 2, 4
                 !
                 v(ir,ipol) = e2 * vs * rho%of_r(ir,ipol) / amag
                 !
                 vtxc = vtxc + v(ir,ipol) * rho%of_r(ir,ipol)
                 !
              END DO
              !
           END IF
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           vtxc = vtxc + v(ir,1) * rho%of_r(ir,1)
           !
        END IF
        !
     END DO
     !
  END IF
  !
  CALL mp_sum(  rhoneg , intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2E10.3)') rhoneg
  !
  ! ... energy terms, local-density contribution
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  ! ... add gradient corrections (if any)
  !
!  CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, v )
 
  !
  ! ... add non local corrections (if any)
  !
!  CALL nonloccorr(rho%of_r, rho_core, etxc, vtxc, v)
  !
  CALL mp_sum(  vtxc , intra_bgrp_comm )
  CALL mp_sum(  etxc , intra_bgrp_comm )
  !
  CALL stop_clock( 'v_xc' )
  !
  RETURN
  !
END SUBROUTINE get_vc

SUBROUTINE GN_ppa (nw,eta, wgrid, ieps, wtilde, wplasma , sysd)

  USE kinds,                ONLY : DP
  USE io_global,        ONLY : stdout, ionode

  integer   nw   !  frequency index for the fitting of rht PPA parameter 
  integer   sysd !  dimensionality of the system 
  real(DP)      wgrid (nw) ! (imaginary) frequency for the PPA parameter
  REAL(DP) :: ieps(3,nw) 
  real(DP)      wtilde, wplasma ! PPA parameters
  real(DP)      eta

  !---- internal variables ------  
  real(DP)      eps0 ,eps1
  integer   iw

  !use Godby Needs plasmon-pole model

  !for 3D
  if (sysd == 3)then
    iw = nw/2
    eps0 =1.d0/(sum(ieps(:,1))/3.d0) 
    eps1 =1.d0/(sum(ieps(:,iw))/3.d0) 
  elseif(sysd ==2 )then
    iw = nw/2
    eps0 =1.d0/((ieps(1,1) +ieps(1,1) )/2.d0) 
    eps1 =1.d0/((ieps(1,iw)+ieps(2,iw))/2.d0) 
  endif
  
  wtilde = (( eps0 - 1.d0 )/(eps0-eps1) - 1.d0) * wgrid(iw)**2
  if (wtilde.gt.0)then
     wtilde = sqrt(wtilde)
  else
     if(ionode)print*, 'wtilde is imaginary! Stop!'
     stop
  endif
  wplasma = (1.d0 - eps0)*wtilde**2
  if (wplasma.gt.0)then
     wplasma = sqrt(wplasma)
  else
     if(ionode)print*, 'wplasma is imaginary! Stop!'
     stop
  endif

  IF(ionode) THEN 

     PRINT *
     PRINT *, ' Godby-Needs Plasmon-pole parameters :  '
     PRINT *, ' wtilde  = ',wtilde 
     PRINT *, ' wplasma = ',wplasma
     PRINT *

     OPEN (30, FILE='epsppa_GN.dat')
     DO jw = 1, nw
       WRITE(30,*) wgrid(jw), & 
         real(1.0d0/(1.0d0+ wplasma**2/(wgrid(jw)**2-(wtilde-(0.0d0,1.0d0)*eta/100)**2))) , &
         aimag(1.0d0/(1.0d0+ wplasma**2/(wgrid(jw)**2-(wtilde-(0.0d0,1.0d0)*eta/100)**2))) !, &
         !epsr(1,jw)
     ENDDO
     CLOSE(30)

     OPEN (30, FILE='epsppa_GN_inv.dat')
     DO jw = 1, nw
       WRITE(30,*) wgrid(jw), & 
         real((1.0d0+ wplasma**2/(wgrid(jw)**2-(wtilde-(0.0d0,1.0d0)*eta/100)**2))) , &
         aimag((1.0d0+ wplasma**2/(wgrid(jw)**2-(wtilde-(0.0d0,1.0d0)*eta/100)**2))) !, &
         !epsr(1,jw)
     ENDDO
     CLOSE(30)


  ENDIF 

END SUBROUTINE GN_ppa

SUBROUTINE print_epsilon (npoles, par, wgrid, nw )

    USE io_global,        ONLY : stdout, ionode
    USE kinds,                ONLY : DP
 
    integer npoles
    real(DP) par(npoles*2)
    integer nw
    real(DP) wgrid(nw)

    integer iw, ipole
    real eta
    complex(DP) epsinv

    eta = 0.01
    
    IF(ionode)THEN
 
    PRINT*, ' I recreating epsilon, and printing it to file'
    PRINT*, par 
    
 
    OPEN (30, FILE='epsilon_new.dat')
    
    DO iw = 1, nw, 1

      epsinv = (1.d0 , 0.d0)
      DO ipole = 1, npoles 
        epsinv = epsinv + par(ipole*2)**2 / ( wgrid(iw)**2-(par(ipole*2-1) - (0.0d0,1.0d0)*eta)**2 )
      ENDDO 

      WRITE(30,*) wgrid(iw), real (1.d0 / epsinv  )
           !real(1.0d0/(1.0d0+ wplasma0**2/(wgrid(iw)**2-(wtilde-(0.0d0,1.0d0)*eta)**2))), &
    ENDDO

    CLOSE(30)
 
    ENDIF

END SUBROUTINE print_epsilon 


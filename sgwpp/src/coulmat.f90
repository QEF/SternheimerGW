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
  USE klist,     ONLY : nks, wk, nelec
  USE lsda_mod,  ONLY : nspin
  USE uspp,      ONLY : okvan
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
  DO ik = 2, nks
     !
     IF ( abs( wk(1) - wk(ik) ) > 1.0d-8 ) &
        CALL errore('grid_build','non unifrom kpt grid', ik )
     !
  ENDDO
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
  alpha = (wmax - wmin) / REAL(nw-1, KIND=DP)
  !
  DO iw = 1, nw
      wgrid(iw) = wmin + (iw-1) * alpha
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

END MODULE

MODULE units_coulmat
  ! ... the units of the files and the record lengths
  SAVE
  INTEGER :: iuncoulmat
  INTEGER :: lrcoulmat
END MODULE units_coulmat

MODULE dielectric
  USE kinds,    ONLY : DP
  SAVE
  REAL(DP) :: qtf, kf
END MODULE dielectric

MODULE control_coulmat
  USE kinds,    ONLY : DP
  SAVE
  LOGICAL  :: do_coulmat, do_fsavg
  INTEGER  :: nbndmin, nbndmax
  REAL(DP) :: degaussfs
END MODULE control_coulmat

!------------------------------
PROGRAM mustar
!------------------------------
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id
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
  USE start_k,     only : nks_start, xk_start, wk_start, &
                          nk1, nk2, nk3, k1, k2, k3
  USE units_coulmat,  ONLY : iuncoulmat, lrcoulmat
  USE dielectric,     ONLY : qtf, kf
  USE control_coulmat, ONLY : do_coulmat, do_fsavg, nbndmin, degaussfs

  IMPLICIT NONE

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! input variables
  !
  INTEGER                 :: nw,nbndmax,nshell,ibndmin,ibndmax
  INTEGER                 :: nk1tmp, nk2tmp, nk3tmp
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift,eta, qmod_par
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc, exst
  !
  NAMELIST / inputpp / prefix, outdir, calculation, nk1, nk2, nk3, qtf, do_coulmat, do_fsavg, nbndmin, kf, degaussfs
  NAMELIST / energy_grid / smeartype,intersmear,intrasmear,wmax,wmin,nbndmax,nw,shift,nshell,eta,ibndmin,ibndmax, qmod_par
  !
  ! local variables
  !
  INTEGER :: ios

  COMPLEX(DP), ALLOCATABLE :: vcnknpkp(:,:,:,:)

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
  calculation  = 'coulmat'
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
  do_coulmat   = .true.
  do_fsavg     = .true.
  qtf          = 1.0
  kf           = 1.0
  degaussfs    = 0.05
 !SHOULD READ FROM INPUT FILE
  CALL input_from_file( )
  !
  ! read input file
  !
  IF (ionode) WRITE( stdout, "( 2/, 5x, 'Reading input file...' ) " )
  ios = 0

  IF ( ionode ) READ (5, inputpp, IOSTAT=ios)

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
  CALL mp_bcast( qtf, ionode_id, world_comm     )
  CALL mp_bcast( do_coulmat, ionode_id, world_comm  )
  CALL mp_bcast( do_fsavg,   ionode_id, world_comm  )

  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )

  CALL read_file
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

  IF (nbndmax == 0) nbndmax = nbnd
  IF (ibndmax == 0) ibndmax = nbnd

  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") trim(calculation)

  CALL start_clock( 'calculation' )

  SELECT CASE ( trim(calculation) )
  !
  CASE ( 'coulmat' )

    ALLOCATE(vcnknpkp(nks,nks,nbnd,nbnd))
    vcnknpkp = (0.0d0,0.d0)

!OPEN DIRECTORY FOR COULOMB MATRIX ELEMENTS
  iuncoulmat = 29
  lrcoulmat  = 2 * nks * nks * nbnd * nbnd
  CALL diropn (29, 'coulmat', lrcoulmat, exst)

  print *,"nk1 nk2 nk3", nk1, nk2, nk3
  print *,"Thomas-Fermi Vector", qtf

  !Calculate coulomb matrix elements stored in struct vcnknpkp
  IF(do_coulmat) THEN
    IF (ionode) WRITE( stdout, "( 5x, 'Calculating Coulomb Matrix Elements' ) " )
        CALL coulmats(vcnknpkp)
  ENDIF


  !Perform fermi surface averaging of matrix elements in vcnknpkp
  IF(do_fsavg) THEN
    IF (ionode) WRITE( stdout, "( 5x, 'Performing Fermi Surface average' ) " )
      CALL fsaverage()
  ENDIF

  CASE DEFAULT
      CALL errore('sgwpp','invalid CALCULATION = '//trim(calculation),1)
  END SELECT
  !
  CALL stop_clock( 'calculation' )
  !
  IF ( ionode ) WRITE( stdout, *  )

  CALL stop_pp ()
END PROGRAM mustar

SUBROUTINE coulmats(vcnknpkp)
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE cell_base,            ONLY : tpiba2, omega, at, alat
  USE wvfct,                ONLY : nbnd, et, npw, igk, npwx,  g2kin, ecutwfc
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss,xk,wk
  USE io_global,            ONLY : ionode, stdout
  USE start_k,              ONLY : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm,me_bgrp
  USE mp,                   ONLY : mp_sum
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE wavefunctions_module, ONLY : evc 
  USE fft_base,             ONLY : dfftp, dffts
  USE gvect,                ONLY : ngm, g, nl
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE io_files,             ONLY : nwordwfc, iunwfc, diropn
  USE funct,                ONLY : xc
  USE grid_subroutines,     ONLY : realspace_grids_info
  USE buffers,              ONLY : save_buffer, get_buffer
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat
  USE control_coulmat,      ONLY : nbndmin

IMPLICIT NONE

REAL        :: lind_eps
COMPLEX(DP) :: norm
REAL(DP)    :: qg2, xq(3), xkp(3)
COMPLEX(DP) :: vcnknpkp(nks,nks,nbnd,nbnd), fnknpkp(dffts%nnr)
COMPLEX(DP) :: psink(dffts%nnr,nbnd), psinpkp(dffts%nnr,nbnd), psi_temp(dffts%nnr)
COMPLEX(DP), ALLOCATABLE :: vc(:,:)
INTEGER     :: ik, ikp, ibnd, jbnd, ig, igp
INTEGER     :: nkp, nkp_abs, ipool
INTEGER     :: ngcoul, ir, iq 
LOGICAL     :: limq, exst
COMPLEX(DP) :: evc1(npwx,nbnd)

ngcoul = 200
ALLOCATE(vc(ngcoul,ngcoul))
!For testing use lindhard function
vc    = (0.0d0,0.0d0)
vcnknpkp = (0.0d0,0.0d0)
ipool = 0
!SHOULD READ FROM INPUT
do iq = 1, nks
!Screened Coulomb interaction
   xq(:) = xk(:,iq)
   print*, "xq:  ", xq, iq
   do ig =1, ngcoul
    qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
    limq = (qg2.lt.eps8) 
    if(.not.limq) then
      vc(ig,ig) = (1.0d0/lind_eps(qg2))*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
    endif
   enddo
   do ik = 1, nks
      CALL ktokpmq (xk(1,ik), xq(1), -1, ipool, nkp, nkp_abs)
      ikp = nkp
      xkp(:) = xk(:,ik)-xq(:)
!read in wavefunctions at k: 
      psink(:,:)   = (0.0d0,0.0d0)
      psinpkp(:,:) = (0.0d0,0.0d0)

      CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
!HL SHOULD ONLY STORE BANDS AROUND FERMI LEVEL
      psink(nls(igk(1:npw)), :) = evc(1:npw,:)
!read in wavefunctions at k' = k-q: 
      CALL gk_sort (xkp(1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL davcio (evc, 2*nwordwfc, iunwfc, ikp, -1 )
      psinpkp(nls(igk(1:npw)),:) = evc(1:npw,:)
!LOOP over bands calculated f functions etc.
      do ibnd = nbndmin, nbnd
         CALL invfft ('Wave', psink(:,ibnd), dffts)
         do jbnd = nbndmin, nbnd
!calculate f_{nk,npkp}(\G)
            psi_temp = (0.0d0, 0.0d0)
            psi_temp = psinpkp(:,jbnd)
            CALL invfft ('Wave', psi_temp(:), dffts)
            fnknpkp = (0.0d0,0.0d0)
            do ir = 1, dffts%nnr  
               fnknpkp(ir) = fnknpkp(ir) +  conjg(psi_temp(ir))*psink(ir,ibnd)
            enddo
            CALL fwfft ('Wave', fnknpkp(:), dffts)
            do ig = 1, ngcoul 
!nolocal fields do igp = 1, ngcoul
               vcnknpkp(ik, ikp, ibnd, jbnd) = vcnknpkp(ik, ikp, ibnd, jbnd) + &
                                               conjg(fnknpkp(nls(ig)))*vc(ig,ig)*fnknpkp(nls(ig))
!               enddo
            enddo
        enddo
     enddo!ibnd
   enddo
enddo

CALL davcio (vcnknpkp, lrcoulmat, iuncoulmat, 1,+1)

END SUBROUTINE coulmats

SUBROUTINE fsaverage()
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE cell_base,            ONLY : tpiba2, omega, at, alat
  USE ener,                 ONLY : ef
  USE io_global,            ONLY : ionode, stdout
  USE start_k,              only : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE klist,                ONLY : nks, nkstot, ngauss, degauss, xk, wk
  USE wavefunctions_module, ONLY : evc 
  USE wvfct,      ONLY : nbnd, et
  USE lsda_mod,   ONLY : nspin
  USE ktetra,     ONLY : ntetra, tetra, ltetra
  USE start_k,              only : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE io_files,             ONLY : nwordwfc, iunwfc, diropn
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat
  USE control_coulmat,      ONlY : degaussfs, nbndmin

IMPLICIT NONE

  COMPLEX(DP) :: vcnknpkp(nks,nks,nbnd,nbnd)
  REAL(DP)    :: enk, enpkp
  REAL(DP)    :: norm, nqs
  REAL(DP)    :: En, DOSofE(2), N0
  REAL(DP)    :: degaussw0, w0g1, w0g2
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  INTEGER     :: ik, ikp, ibnd, jbnd, ig, igp
!IO
  LOGICAL     :: exst
!to put in coul struct:
  REAL(DP)    :: mu, mustar, debye_e


  CALL davcio (vcnknpkp, lrcoulmat, iuncoulmat, 1,-1)

!First calculate dosef 
   !dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
   !   N(Ef) in the equation for lambda is the DOS per spin
   !dosef = dosef / two
   tetra=.true.
    WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
        &        "ngauss,degauss=",i4,f12.6/)') ngauss, degauss

!CALCULATE Dos at Fermi Level gaussian:
   CALL dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, ef, DOSofE)
!use tetrahedron method:
  !CALL dos_t(et, nspin, nbnd, nks, ntetra, tetra, ef, DOSofE)
  !want density of states per spin.
   N0 = DOSofE(1)/2.d0
   !degaussw0 = 0.025 !eV
   !degaussw0 = degaussw0/RYTOEV !eV
  if (degaussfs.eq.0.0d0) then
     degaussw0 = degauss
  else
     degaussw0 = degaussfs
  endif

   En      = ef
   nqs     = dble(nk1*nk2*nk3)
   print*, ef*rytoev, N0/rytoev
   print*, "degauss fs"
  mu = 0.0d0
  do ik = 1, nks
   if (mod(ik,20).eq.0) print*, "xk:  ", ik
   do ikp = 1, nks
!  do ik = ikstart, ikstop
!   do ikp = ikstart, ikstop
    do ibnd = nbndmin, nbnd
     do jbnd = nbndmin, nbnd
        enk = (et(ibnd, ik) - ef)
        enpkp = (et(jbnd, ikp) - ef)
        w0g1 = w0gauss ( enk / degaussw0, 0) / degaussw0
        w0g2 = w0gauss ( enpkp / degaussw0, 0) / degaussw0
!       mu = mu + (1.0d0/N0)*vcnknpkp(ik,ikp,ibnd,jbnd)*w0g1*w0g2*wk(ik)*wk(ikp)
!       Again we want per spin.
        mu = mu + (1.0d0/N0)*vcnknpkp(ik,ikp,ibnd,jbnd)*w0g1*w0g2*(wk(ik)/2.0)
     enddo
    enddo
   enddo
  enddo
!CALL mp_sum(mu, inter_pool_comm)
!for aluminum we take debye temperature is 428 K convert to rydberg.
  debye_e = 0.036
!Factors not included when we calculate V^{c}_{nkn'k'}.
  mu = mu/omega/nqs
  print*, "\mu", mu
  print*, "\mu^{*}", mu/(1+mu*log((12.43/RYTOEV)/debye_e))
END SUBROUTINE fsaverage

REAL(DP) FUNCTION lind_eps(qg2)
  USE kinds,       ONLY : DP
  USE dielectric,  ONLY : qtf, kf
  IMPLICIT NONE
  REAL(DP) :: qg2, x
!qtf tommy fermi wave vector for Al... (according to my calculations)

  x = qg2/(2*kf)
  lind_eps = 1.0d0 + (qtf)**2.0/(2.0*qg2**2)*(1.0d0+(1.0d0/(2.0*x))*(1-x**2)*log(abs((1+x)/(1-x))))
  !lind_eps = 1.0d0 
  RETURN
END FUNCTION
!Construct screened coulomb interaction
!rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
!do ig = 1, ngmpol
!  qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
!  limq = (qg2.lt.eps8) 
!  if(.not.limq) then
!     do igp = 1, ngmpol
!        vc(ig, igp) = scrcoul(ig,igp,1)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
!     enddo
!  endif
!  qg = sqrt(qg2)
!  spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
!  if(.not.limq) then
!      do igp = 1, ngmpol
!         vc = scrcoul(ig,igp,1)*dcmplx(spal, 0.0d0)
!      enddo
!  else
!      scrcoul_g_R(ig,igp,1) = scrcoul_g_R(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
!  endif
!enddo
!         do ir = 1, dffts%nnr
!            !fnknpkp(ir) = fnknpkp(ir) + conjg(psi_temp(ir))*psi_temp(ir)
!             fnknpkp(ir) = fnknpkp(ir) + conjg(psink(ir,ibnd))*psink(ir,ibnd)
!             norm   = norm + conjg(psink(ir,ibnd))*psink(ir,ibnd)/(dffts%nnr)
!         enddo
!         CALL fwfft ('Wave', fnknpkp(:), dffts)
!         print*, sum(fnknpkp(nls(igk(1:npw)))), fnknpkp(1)
          !print*, sum(fnknpkp(nls(igk(1:npw)))), fnknpkp(1)

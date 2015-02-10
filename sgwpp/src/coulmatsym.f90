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
  LOGICAL  :: do_coulmat, do_fsavg, do_lind
  INTEGER  :: nbndmin, nbndmax, ngcoul
  REAL(DP) :: degaussfs, debye_e
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
  USE start_k,     ONLY : nks_start, xk_start, wk_start, &
                          nk1, nk2, nk3, k1, k2, k3
  USE units_coulmat,   ONLY : iuncoulmat, lrcoulmat, lrcoul, iuncoul
  USE dielectric,      ONLY : qtf, kf
  USE control_coulmat, ONLY : do_coulmat, do_fsavg, nbndmin, degaussfs, ngcoul, do_lind, &
                              debye_e
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
  NAMELIST / inputpp / prefix, outdir, calculation, nk1, nk2, nk3, qtf, do_coulmat, do_fsavg, nbndmin, kf, degaussfs, ngcoul, do_lind, debye_e
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
  do_lind      = .false.
  qtf          = 1.0
  kf           = 1.0
  degaussfs    = 0.05
  debye_e      = 0.0026

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

!OPEN COULOMB directory
  iuncoul = 28
  lrcoul = 2 * ngcoul * ngcoul
  CALL diropn (iuncoul, 'coul', lrcoul, exst)

  print *,"nk1 nk2 nk3", nk1, nk2, nk3
  print *,"Thomas-Fermi Vector", qtf

  !Calculate Coulomb matrix elements stored in struct vcnknpkp
  IF(do_coulmat) THEN
    IF (ionode) WRITE( stdout, "( 5x, 'Calculating Coulomb Matrix Elements' ) " )
        CALL coulmats(vcnknpkp)
  ENDIF

  !Perform Fermi surface averaging of matrix elements in vcnknpkp
  IF(do_fsavg) THEN
    IF (ionode) WRITE( stdout, "( 5x, 'Performing Fermi Surface average' ) " )
      CALL fsaverage()
  ENDIF

  CASE DEFAULT
      CALL errore('sgwpp','invalid CALCULATION = '//trim(calculation),1)
  END SELECT
  !
  CALL stop_clock( 'calculation' )

  CALL print_clock( 'calculation' )
  !
  IF ( ionode ) WRITE( stdout, *  )

  CALL stop_pp ()
END PROGRAM mustar

SUBROUTINE load_coul()

IMPLICIT

   vcnknpkp = (0.0d0,0.0d0)
   ipool = 0
   xq(:) = xk(:,iq)
   vc    = (0.0d0,0.0d0)
   if (mod(iq,30).eq.0) print*, "xq:  ", iq
   if(.not.do_lind) then 
      print*, "reading xq:  ", xq, iq
      CALL davcio (vc, lrcoul, iuncoul, iq, -1)
!what is stored is: eps^-1(q)-1
      do ig = 1, ngcoul
         vc(ig,ig) = vc(ig,ig) + 1.0d0
      enddo 
   endif
   do ig = 1, ngcoul
      qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
      limq = (qg2.lt.eps8) 
      if (limq) cycle
      if (.not.do_lind) then
         do igp = 1, ngcoul
            vc(ig,igp) = vc(ig,igp)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
         enddo
      else
            vc(ig,ig) = (1.0d0/lind_eps(qg2))*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
      endif
   enddo
END SUBROUTINE load_coul

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
  USE control_coulmat,      ONlY : degaussfs, nbndmin, debye_e

IMPLICIT NONE

  COMPLEX(DP) :: vcnknpkp(nks,nks,nbnd,nbnd)
  REAL(DP)    :: enk, enpkp
  REAL(DP)    :: norm, nqs
  REAL(DP)    :: En, DOSofE(2), N0
  REAL(DP)    :: degaussw0, w0g1, w0g2
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  COMPLEX(DP) :: psink(dffts%nnr,nbnd), psinpkp(dffts%nnr,nbnd), psi_temp(dffts%nnr), fnknpkp(dffts%nnr)
  INTEGER     :: ik, ikp, ibnd, jbnd, ig, igp
!IO
  LOGICAL     :: exst
!to put in coul struct:
  REAL(DP)    :: mu, mustar
  CALL davcio (vcnknpkp, lrcoulmat, iuncoulmat, 1,-1)
  !First calculate dosef: 
  !dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
  !N(Ef) in the equation for lambda is the DOS per spin
  !dosef = dosef / two
  WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
       &        "ngauss,degauss=",i4,f12.6/)') ngauss, degauss
  !CALCULATE Dos at Fermi Level gaussian:
   CALL dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, ef, DOSofE)
  !tetra=.true.
  !use tetrahedron method:
  !CALL dos_t(et, nspin, nbnd, nks, ntetra, tetra, ef, DOSofE)
  !want density of states per spin.
  N0 = DOSofE(1)/2.d0
  if (degaussfs.eq.0.0d0) then
     degaussw0 = degauss
  else
     degaussw0 = degaussfs
  endif
  En      = ef
  nqs     = dble(nk1*nk2*nk3)
  mu = 0.0d0

!Loop over IBZ on kpoints:
  do ik = 1, nks
   if (mod(ik,2).eq.0) print*, "xk:  ", ik
   do iq = 1, nqs
!do S = nsymops:
!find index for ikp
!need to modify this:
      CALL ktokpmq (xk(1,ik), xq(1), -1, ipool, nkp, nkp_abs)
!on a gen'l grid ktokpmq should also return the symmetry operation which rotates
!xkp to the brillouin zone. then we keep track of that symmetry else and use it
!to invert one of kpoints in the IBZ. Keeping with the idea that we only need
!IBZk and IBZq.
      ikp = nkp
      xkp(:) = xk(:,ik)-xq(:)
!zero wavefunctions at k and k': 
      psink(:,:)   = (0.0d0,0.0d0)
      psinpkp(:,:) = (0.0d0,0.0d0)

      CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
      psink(nls(igk(1:npw)), :) = evc(1:npw,:)

      CALL gk_sort (xkp(1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL davcio (evc, 2*nwordwfc, iunwfc, ikp, -1 )
      psinpkp(nls(igk(1:npw)),:) = evc(1:npw,:)

!where to parallelize? I think over q.
      do ibnd = nbndmin, nbnd
        CALL invfft ('Wave', psink(:,ibnd), dffts)
        do jbnd = nbndmin, nbnd

            psi_temp = (0.0d0, 0.0d0)
            psi_temp = psinpkp(:,jbnd)
            CALL invfft ('Wave', psi_temp(:), dffts)

            fnknpkp = (0.0d0,0.0d0)
            do ir = 1, dffts%nnr  
               fnknpkp(ir) = fnknpkp(ir) +  conjg(psi_temp(ir))*psink(ir,ibnd)
            enddo

!calculate f_{nk,npkp}(\G)
            CALL fwfft ('Wave', fnknpkp(:), dffts)

           enk = (et(ibnd, ik) - ef)
           enpkp = (et(jbnd, ikp) - ef)

           w0g1 = w0gauss ( enk / degaussw0, 0) / degaussw0
           w0g2 = w0gauss ( enpkp / degaussw0, 0) / degaussw0
!       Again we want per spin.
            if(.not.do_lind) then
              do ig = 1, ngcoul 
                 do igp = 1, ngcoul
                   vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,igp)*fnknpkp(nls(ig))
                 enddo
              enddo
            else
              do ig = 1, ngcoul 
                   vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,ig)*fnknpkp(nls(ig))
              enddo
            endif
            mu = mu + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)
        enddo!jbnd
      enddo!ibnd
    enddo!iq
  enddo!ik

!CALL mp_sum(mu, inter_pool_comm)
!Factors not included when we calculate V^{c}_{nkn'k'}.

  mu = mu/omega/nqs
  print*, "Ef", ef*rytoev, "N(0)", N0/rytoev
  print*, "debye temp Ry", debye_e
  print*, "\mu", mu
  print*, "\mu^{*}", mu/(1+mu*log((ef)/debye_e))

END SUBROUTINE fsaverage

REAL(DP) FUNCTION lind_eps(qg2)
  USE kinds,       ONLY : DP
  USE dielectric,  ONLY : qtf, kf
  IMPLICIT NONE
  REAL(DP) :: qg2, x
  x = qg2/(2*kf)
  lind_eps = 1.0d0 + (qtf)**2.0/(2.0*qg2**2)*(1.0d0+(1.0d0/(2.0*x))*(1-x**2)*log(abs((1+x)/(1-x))))
  RETURN
END FUNCTION


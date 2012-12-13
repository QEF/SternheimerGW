!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the GW program
!  
MODULE modes
  USE kinds,  ONLY : DP
  ! ... The variables needed to describe the modes and the small group of q
  SAVE
  !
  INTEGER :: irgq(48), nsymq, irotmq, nirr, nmodes
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  ! number of irreducible representations contained in the dynamical matrix
  ! number of modes
  ! number of crystal sym.ops. for q=0 
  INTEGER, ALLOCATABLE, TARGET :: npert(:) !3 * nat )
  ! the number of perturbations per IR
  INTEGER :: npertx
  ! max number of perturbations per IR
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  COMPLEX (DP), POINTER :: &
       u(:,:),                     &!  3 * nat, 3 * nat),
       ubar(:),                    &!  3 * nat), &
       t(:,:,:,:),                 &! npertx, npertx, 48,3 * nat),
       tmq(:,:,:)                   ! npertx, npertx, 3 * nat)
  ! the transformation modes patterns
  ! the mode for deltarho
  ! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa
  LOGICAL :: &
       minus_q,  &    !  if .TRUE. there is the symmetry sending q<->-q
       invsymq        !  if .TRUE. the small group of q has inversion

  CHARACTER(15), ALLOCATABLE :: name_rap_mode(:) ! symmetry type of each mode
  !     
END MODULE modes
 
MODULE dynmat
  USE kinds, ONLY :  DP
  !
  ! ... The dynamical matrix 
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE :: &
       dyn00(:,:),           &! 3 * nat, 3 * nat),
       dyn(:,:),             &! 3 * nat, 3 * nat)
       dyn_rec(:,:)           ! 3 * nat, 3 * nat)
  ! the initial dynamical matrix
  ! the dynamical matrix
  ! the contribution of each representation to the dynamical matrix
  REAL (DP), ALLOCATABLE :: &
       w2(:)                  ! 3 * nat)
  ! omega^2
  !
END MODULE dynmat
!

! MODULE gwconst
!  USE kinds, ONLY : DP
!  SAVE
!  COMPLEX(DP),PARAMETER :: ci = (0.000000D0, 1.0000000D0)
!  COMPLEX(DP),PARAMETER :: czero = (0.000000D0, 0.0000000D0)
! END MODULE gwconst

MODULE qpoint
  USE kinds, ONLY :  DP
  USE parameters, ONLY : npk
  !
  ! ... The q point
  !
  SAVE
  !
  INTEGER, POINTER :: igkq(:)     ! npwx
  ! correspondence k+q+G <-> G
  INTEGER :: nksq, npwq
  ! the real number of k points
  ! the number of plane waves for q
  INTEGER, ALLOCATABLE :: ikks(:), ikqs(:)
  ! the index of k point in the list of k
  ! the index of k+q point in the list of k
  REAL (DP) :: xq(3)
  ! the coordinates of the q point
  COMPLEX (DP), ALLOCATABLE :: eigqts(:) ! nat)
  ! the phases associated to the q

  !
END MODULE qpoint
!
!
MODULE eqv
  USE kinds, ONLY :  DP
  !
  ! ... The wavefunctions at point k+q 
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  !
  ! ... The variable describing the linear response problem 
  !
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), drhoscfs (:,:), dpsip(:,:), dpsim(:,:)
  ! the product of dV psi
  ! dpsi the change of the wavefunction
  ! dpsip change of wavefunction for positive frequency
  ! dpsim change of wavefunction for negative frequency

  COMPLEX (DP), ALLOCATABLE :: dvbare(:)

  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin),
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the derivative of the xc potential
  ! the local potential at q+G
  REAL (DP), ALLOCATABLE :: eprec(:,:) ! needed for preconditioning
  !
END MODULE eqv
!
!
MODULE efield_mod
  USE kinds, ONLY :  DP
  !
  ! ... the variables for the electric field perturbation
  !  
  SAVE
  !
  REAL (DP) :: epsilon (3, 3)
  REAL (DP), ALLOCATABLE :: &
       zstareu(:,:,:),       &! 3, 3, nat),
       zstarue(:,:,:)         ! 3, nat, 3)
  ! the dielectric constant
  ! the effective charges Z(E,Us) (E=scf,Us=bare)
  ! the effective charges Z(Us,E) (Us=scf,E=bare)
  COMPLEX (DP), ALLOCATABLE :: &
       zstareu0(:,:),        &! 3, 3 * nat),
       zstarue0(:,:),        &! 3 * nat, 3)
       zstarue0_rec(:,:)      ! 3 * nat, 3)
  ! the effective charges
  !
END MODULE efield_mod
!
!
MODULE nlcc_gw
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for non-linear core correction
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE, TARGET :: drc(:,:) ! ngm, ntyp)
  ! contain the rhoc (without structure fac) for all atomic types
  LOGICAL :: nlcc_any
  ! .T. if any atom-type has nlcc
  !
END MODULE nlcc_gw
!
!
MODULE gc_gw
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),              &! 3, nrxx, nspin),
       gmag(:,:,:),              &! 3, nrxx, nspin),
       vsgga(:),                 &! nrxx
       segni(:),                 &! nrxx
       dvxc_rr(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_sr(:,:,:),           &! nrxx, nspin, nspin),
       dvxc_ss(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_s(:,:,:)              ! nrxx, nspin, nspin)
  !
  ! in the noncollinear case gmag contains the gradient of the magnetization
  ! grho the gradient of rho+ and of rho-, the eigenvalues of the spin density
  ! vsgga= 0.5* (V_up-V_down) to be used in the calculation of the change
  ! of the exchange and correlation magnetic field.
  ! gradient of the unpert. density
  !
  ! derivatives of the E_xc functiona
  ! r=rho and s=|grad(rho)|
  !
END MODULE gc_gw
!
!
MODULE gwus
  USE kinds, ONLY :  DP
  USE becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... program with the US pseudopotentials
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       alphasum(:,:,:,:),   &! nhm*(nhm+1)/2,3,nat,nspin)
                             ! used to compute modes
       dpqq(:,:,:,:)         ! (nhm, nhm, 3, ntyp) 
  ! alphasum contains \sum_i <psi_i| d/du (|\beta_n><beta_m|) | psi_i> + (m-n)
  ! dipole moment of each Q 
  COMPLEX (DP), ALLOCATABLE :: &
       int1(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2(:,:,:,:,:),     &! nhm, nhm, 3,nat, nat),&
       int3(:,:,:,:,:),     &! nhm, nhm, npert, nat, nspin),&
       int3_paw(:,:,:,:,:), &! nhm, nhm, npert, nat, nspin),&
       int4(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
       int5(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nat),&
       int1_nc(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2_so(:,:,:,:,:,:),   &! nhm, nhm, 3, nat,nat,nspin),&
       int3_nc(:,:,:,:,:),     &! nhm, nhm, npert, nat, nspin),&
       int4_nc(:,:,:,:,:,:),   &! nhm, nhm, 3, 3, nat, nspin),&
       int5_so(:,:,:,:,:,:,:), &! nhm*(nhm+1)/2, 3, 3, nat, nat, nspin),&
!
!  These variables contains the five integrals defined in PRB 64, 35118 (2001)
!  int1 -> \int V_eff d/du (Q) d^3r
!  int2 -> \int d/du (V_loc) Q d^3r
!  int3 -> \int d\du (V_Hxc) Q d^3r
!  int4 -> \int V_eff d^2/dudu (Q) d^3r
!  int5 -> \int d/du (V_loc) d/du (Q) d^3r
!  int3_paw contains d/du (D^1-\tilde D^1) 
! 
       becsum_nc(:,:,:,:),     &! nhm*(nhm+1)/2,nat,npol,npol)
       becsumort(:,:,:,:),     &! nhm*(nhm+1)/2,nat,nspin,3*nat)
       alphasum_nc(:,:,:,:,:), &! nhm*(nhm+1)/2,3,nat,npol,npol)
       dpqq_so(:,:,:,:,:)       ! nhm, nhm, nspin, 3, ntyp 
!
!  becsum contains \sum_i <\psi_i | \beta_n><\beta_m| \psi_i > + (m-n)
!  besumort contains alphasum+\sum_i <\psi_i | \beta_n><\beta_m| \delta \psi_i >
!  dpqq_so dipole moment of each Q multiplied by the fcoef factors
!
  type (bec_type),  ALLOCATABLE, TARGET :: &
       becp1(:)              ! (nksq); (nkbtot, nbnd)
  !
  ! becp1 contains < beta_n | \psi_i > 
  !
  type (bec_type),  ALLOCATABLE, TARGET :: &
       alphap(:,:)           ! nkbtot, nbnd, 3, nksq)
  !
  ! alphap contains < d\du (\beta_n) | psi_i> 
  !
END MODULE gwus
!
!
MODULE partial
  USE kinds, ONLY :  DP
  !
  ! ... the variables needed for partial computation of dynamical matrix
  !
  SAVE
  !  
  INTEGER, ALLOCATABLE :: &
       comp_irr(:),           &! 3 * nat ),
       done_irr(:),           &! 3 * nat), &
       list(:),               &! 3 * nat),
       atomo(:)                ! nat)
  ! if 1 this representation has to be computed
  ! if 1 this representation has been done
  ! a list of representations (optionally read in input)
  ! list of the atoms that moves
  INTEGER :: nat_todo, nrapp
  ! number of atoms to compute
  ! The representation to do
  LOGICAL :: all_comp
  ! if .TRUE. all representation have been computed
  !
END MODULE partial
!
MODULE gamma_gamma
  INTEGER, ALLOCATABLE :: &
           has_equivalent(:),  &  ! 0 if the atom has to be calculated
           with_symmetry(:),   &  ! calculated by symmetry
           n_equiv_atoms(:),   &  ! number of equivalent atoms
           equiv_atoms(:,:)       ! which atoms are equivalent

  INTEGER :: n_diff_sites,    &   ! Number of different sites
             nasr                 ! atom calculated with asr
                                  !
  LOGICAL :: asr                  ! if true apply the asr

END MODULE gamma_gamma
!
MODULE control_gw
  USE kinds, ONLY :  DP
  USE parameters, ONLY: npk
  !
  ! ... the variables controlling the GW run
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter = 100
  ! maximum number of iterations
  INTEGER :: niter_gw, nmix_gw, nbnd_occ(npk), &
             start_irr, last_irr, current_iq, start_q, last_q
  ! maximum number of iterations (read from input)
  ! mixing type
  ! occupated bands in metals
  ! starting representation
  ! initial representation
  ! last representation of this run
  ! current q point
  ! initial q in the list, last_q in the list
  real(DP) :: tr2_gw, tr2_green
  !
  real(DP) :: eta
  ! threshold for gw calculation
  REAL (DP) :: alpha_mix(maxter), time_now, alpha_pv
  ! the mixing parameter
  ! CPU time up to now
  ! the alpha value for shifting the bands
  CHARACTER(LEN=10)  :: where_rec='no_recover'! where the gw run recovered
  CHARACTER(LEN=256) :: flmixdpot, tmp_dir_gw
  INTEGER :: rec_code, &   ! code for recover
             rec_code_read=-1000 ! code for recover. Not changed during the run

  INTEGER :: maxter_green

  LOGICAL :: lgamma,      &! if .TRUE. this is a q=0 computation
             lgamma_gamma,&! if .TRUE. this is a q=0 computation with k=0 only 
             convt,       &! if .TRUE. the GW has converged
             epsil,       &! if .TRUE. computes dielec. const and eff. charges
             done_epsil=.FALSE.,  &! .TRUE. when diel. constant is available
             trans,       &! if .TRUE. computes GW
             elgw,        &! if .TRUE. computes electron-gw interaction coeffs
             zue,         &! if .TRUE. computes eff. charges as induced polarization
             done_zue=.FALSE., &! .TRUE. when the eff. charges are available
             zeu,         &! if .TRUE. computes eff. charges as induced forces
             done_zeu=.FALSE., &! .TRUE. when the eff. charges are available
             recover,     &! if .TRUE. the run restarts
             ext_restart, &! if .TRUE. there is a restart file
             ext_recover, &! if .TRUE. there is a recover file
             lrpa,        &! if .TRUE. calculates the RPA dielectric constant
             lnoloc,      &! if .TRUE. calculates the dielectric constant
                           ! neglecting local field effects
             search_sym,  &! if .TRUE. search the mode symmetry
             lnscf,       &! if .TRUE. the run makes first a nscf calculation
             ldisp,       &! if .TRUE. the run calculates full GW dispersion
             reduce_io,   &! if .TRUE. reduces needed I/O
             done_bands,  &! if .TRUE. the bands have been calculated
             bands_computed=.FALSE., & ! if .TRUE. the bands were computed 
                                       ! in this run
             nogg,        &! if .TRUE. gamma_gamma tricks are disabled
             u_from_file=.FALSE.,  & ! if true the u are on file
             recover_read=.FALSE., & ! if true the recover data have been read
             all_done, &      ! if .TRUE. all representations have been done
             modielec, & ! if .TRUE. uses a model dielectric function to calculate W.
             do_coulomb, &
             do_sigma_c, &
             do_sigma_exx, &
             do_green, &
             do_sigma_matel,&
             do_q0_only,&
             godbyneeds,&
             padecont,&
             cohsex,&
             multishift,&
             do_sigma_extra

END MODULE control_gw
!
!
MODULE freq_gw
  !
  USE kinds,   ONLY : DP
  !
  SAVE
  ! ... the variables for computing frequency dependent dielectric constant
  LOGICAL :: fpol ! if .TRUE. dynamic dielectric constant is computed
  INTEGER, PARAMETER :: nfsmax=50  ! # of maximum frequencies
  INTEGER :: nfs                   ! # of frequencies
  REAL (KIND=DP) :: fiu(nfsmax)    ! values  of frequency
  !variables for convolution
  INTEGER :: nwcoul, nwgreen, nwalloc, nwsigma

  !The wsigmamin, wsigmamax, etc is currently being set in freqbins. 
  !I will change this so that it becomes a user defined option in the punchcard
  !with default values eventually. 

  REAL(DP) :: wsigmamin, wsigmamax, deltaw, wcoulmax, wgreenmin, wgreenmax
  REAL(DP), ALLOCATABLE :: wtmp(:), wcoul(:), wgreen(:), wsigma(:) 
  INTEGER, ALLOCATABLE :: ind_w0mw (:,:), ind_w0pw (:,:)
  REAL(DP) :: plasmon
  REAL(DP) :: greenzero

END MODULE freq_gw
!
!
MODULE units_gw
  !
  ! ... the units of the files and the record lengths
  !
  SAVE
  !
  INTEGER :: &
       iuwfc, lrwfc, iuvkb, iubar, lrbar, iuebar, lrebar, iudwf, iupsir, &
       lrdwf, iudrhous, lrdrhous, iudyn, iupdyn, iunrec, iudvscf, iudrho, &
       lrdrho, iucom, lrcom, iudvkb3, lrdvkb3, iuncoul, iungreen, iunsigma, &
       iudwfm, iudwfp, lrgrn, lrcoul, lrsigma, iuwfcna, iunsex, lrsex, &
       lrresid, lralphabeta, iunresid, iunalphabeta, iunsigext, lrsigext


  ! iunit with the wavefunctions
  ! the length of wavefunction record
  ! unit with vkb
  ! unit with the part DV_{bare}
  ! length of the DV_{bare}
  ! unit with D psi
  ! unit with evc in real space
  ! length of D psi record
  ! the unit with the products
  ! the length of the products
  ! the unit for the dynamical matrix
  ! the unit for the partial dynamical matrix
  ! the unit with the recover data
  ! the unit where the delta Vscf is written
  ! the unit where the delta rho is written
  ! the length of the deltarho files
  ! the unit of the bare commutator in US case
  ! the length  of the bare commutator in US case
  ! the screened coulomb potential

  logical, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_gw
!
!
MODULE output
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=256) :: fildyn, fildvscf, fildrho
  ! output file for the dynamical matrix
  ! output file for deltavscf
  ! output file for deltarho
  !
END MODULE output
!
!
MODULE disp
   !
   USE kinds, ONLY: DP
   !
   SAVE
   !
   INTEGER, PARAMETER :: nqmax = 1000
   !
   INTEGER :: nq1, nq2, nq3
    ! number of q-points in each direction
   INTEGER :: iq1, iq2, iq3
    ! specific q point from the regular grid
    ! (i.e., iq1/nq1,iq2/nq2,iq3/nq3)
   INTEGER :: nqs
    ! number of q points to be calculated 
   REAL (DP), ALLOCATABLE :: x_q(:,:)
    ! coordinates of the q points
   INTEGER, ALLOCATABLE :: done_iq(:)
    ! if 1 this q point has been already calculated
   INTEGER, ALLOCATABLE :: comp_iq(:)
    ! if 1 this q point has to be calculated
   INTEGER, ALLOCATABLE :: rep_iq(:)
    ! number of irreducible representation per q point
   INTEGER, ALLOCATABLE :: done_rep_iq(:,:)
    ! which representation have been already done in each q
   INTEGER, ALLOCATABLE :: nsymq_iq(:)
    ! dimension of the small group of q
   INTEGER, ALLOCATABLE :: comp_irr_iq(:,:)
    ! for each q, comp_irr. Used for image parallelization 
   INTEGER, ALLOCATABLE :: npert_iq(:,:)
    ! for each q, the number of perturbation of each irr
   REAL(DP) , ALLOCATABLE :: wq(:)

  !HL Variables required for shuffling the k/q grid. 
   INTEGER, ALLOCATABLE :: gmap(:,:)
   REAL(DP) :: g0vec(3,27)
   REAL(DP), ALLOCATABLE :: eval_occ(:,:) ! array of eigenvalues after folding.

 ! kpoints: default false.
 ! if true specifies the k-points we want to look 
 ! at in the brillouin zone specified by user in punch card.  
   LOGICAL  :: kpoints 
   REAL(DP) :: xk_kpoints(3,10)
   INTEGER  :: num_k_pts
   INTEGER  :: w_of_q_start
  

END MODULE disp

MODULE coulomb_app
CONTAINS
SUBROUTINE spheric_coulomb(ngmpol, nfs, diel_matrix, xq_coul)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE cell_base,     ONLY : alat, tpiba2, omega
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE disp,          ONLY : nqs, nq1, nq2, nq3
  IMPLICIT NONE
!SUBROUTINE SPHERIC CUT  COULOMB POTENTIAL Appliead to v_q.
!Spencer/Alavi truncation of the bare coulomb interaction
![PRB 77,193110 (2008)]
  LOGICAL     :: limq
  REAL(DP)    :: xq_coul(3)
  REAL(DP)    :: qg, rcut, spal
  INTEGER     :: ig, igp, iw, ngmpol, nfs
  INTEGER     :: unf_recl, recl, ios
  COMPLEX(DP), INTENT(INOUT) :: diel_matrix(ngmpol,ngmpol,nfs) 
  !COMPLEX(DP) :: diel_matrix(:,:,:) 
  REAL(DP)    :: qg2, qg2coul


rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))

DO iw = 1, nfs
   DO ig = 1, ngmpol
       qg2 = (g(1,ig) + xq_coul(1))**2 + (g(2,ig) + xq_coul(2))**2 + (g(3,ig) + xq_coul(3))**2
       !if(qg2.lt.eps8) limq =.true.
       limq = (qg2.lt.eps8) 
       IF(.not.limq) then
           DO igp = 1, ngmpol
              diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
           ENDDO
       ELSE 
           if(ig.ne.1) then
              DO igp = 1, ngmpol
                 diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
              ENDDO
           endif
       ENDIF

       qg = sqrt(qg2)
       spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)

!Normal case using truncated coulomb potential.
       if(.not.limq) then
          do igp = 1, ngmpol
              diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx(spal, 0.0d0)
          enddo
       else
!should only occur case iq->0, ig = 0 use vcut (q(0) = (4pi*e2*Rcut^{2})/2
             write(6,'("Taking Limit.")')
             write(6,*) (fpi*e2*(rcut**2))/2.0d0
             write(6,*) ig, iw
             write(6,*) g(:, ig)
             !for omega=0,q-->0, G=0 the real part of the head of the dielectric matrix should be real
             !we enforce that here:
         if(iw.eq.1) then
            diel_matrix(ig, igp, iw) = real(diel_matrix(ig,igp,iw))
         endif
         do igp = 1, ngmpol
            diel_matrix(ig, igp, iw) = diel_matrix(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
         enddo
       endif
   ENDDO !ig
ENDDO!iw
END SUBROUTINE
END MODULE coulomb_app

MODULE gwsigma
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : omega, alat
  USE qpoint,      ONLY : xq, igkq
  
  PUBLIC :: fft6_g2r
  
  SAVE

  COMPLEX (DP), ALLOCATABLE :: scrcoul(:,:,:,:)
  COMPLEX (DP), ALLOCATABLE :: green(:,:)

! HL self energy is a huge quantity!
  COMPLEX (DP), ALLOCATABLE :: sigma_ex(:,:)
  COMPLEX (DP), ALLOCATABLE :: sigma_g_ex(:,:)

  COMPLEX (DP), ALLOCATABLE :: sigma(:,:,:)
  COMPLEX (DP), ALLOCATABLE :: sigma_g(:,:,:)

  INTEGER, ALLOCATABLE ::         nlsig(:)
  INTEGER, TARGET, ALLOCATABLE :: nlsex(:)
  INTEGER, TARGET, ALLOCATABLE :: nlsco(:)
  
! Cutoff for the sigma + exchange/correlation.
  REAL(DP) :: ecutsig
  REAL(DP) :: ecutpol
  REAL(DP) :: ecutgrn
  REAL(DP) :: ecutsex
  REAL(DP) :: ecutsco

  REAL(DP) :: gcutmsig
  INTEGER  :: nbnd_sig
  INTEGER  :: ngmsig, ngmsco, ngmsex, ngmpol, ngmgrn

! Real space mesh for description of self-energy.
  INTEGER :: nr1sig, nr2sig, nr3sig, nrsig
  INTEGER :: nr1sco, nr2sco, nr3sco, nrsco
  INTEGER :: nr1sex, nr2sex, nr3sex, nrsex

  CONTAINS

  SUBROUTINE fft6_g2r(ngmtmp, nrtmp, nltmp, f_g, f_r, corx)
   USE qpoint,  ONLY :  npwq, igkq 

   IMPLICIT NONE
   INTEGER                          ::  ios
   INTEGER                          ::  ig, igp, ir, irp
   COMPLEX(DP)                      ::  czero

  ! corx: correlation or exchange grid
  ! gorc: green's function (transforms igkq), or coulomb

   INTEGER, INTENT(IN) ::  corx
   INTEGER             ::  nr1tmp, nr2tmp, nr3tmp, nrtmp, ngmtmp, ngmtmp1
   COMPLEX(DP)         :: f_g (ngmtmp, ngmtmp)
   COMPLEX(DP)         :: f_r (nrtmp,nrtmp)
   COMPLEX(DP)         :: aux (nrtmp)
   INTEGER             :: nltmp(:)

     if (corx.eq.1) then
       nr1tmp = nr1sco
       nr2tmp = nr2sco
       nr3tmp = nr3sco
     else if (corx.eq.2) then 
       nr1tmp = nr1sex
       nr2tmp = nr2sex
       nr3tmp = nr3sex
     else
        WRITE(6,'("error fft6_g2r_coulomb what to do?")')
        STOP
     endif

  ! if ngmsex exceeds npwq we are going beyond the highest G-vector describing the wavefunction. 
     czero = (0.0d0, 0.0d0)
     f_r(:,:) = czero
     do ig = 1, ngmtmp
        aux(:) = czero
        WRITE(6,'("FFT G-prime")')
        do igp = 1, ngmtmp
           aux(nltmp(igp)) = f_g(ig,igp)
        enddo
        WRITE(6,*)nr1tmp, nr2tmp, nr3tmp, nrtmp
        call cft3s (aux, nr1tmp, nr2tmp, nr3tmp, nr1tmp, nr2tmp, nr3tmp, +1)
        do irp = 1, nrtmp
           f_r(ig, irp) = aux(irp) / omega
        enddo
     enddo

  ! the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
  ! following the convention set in the paper
  ! [because the standard transform is sum_G f(G) exp(iGr) ]
     WRITE(6,'("FFT G")')
     do irp = 1, nrtmp
        aux = czero
        do ig = 1, ngmtmp
           aux(nltmp(ig)) = conjg( f_r(ig,irp) )
        enddo
        call cft3s (aux, nr1tmp, nr2tmp, nr3tmp, nr1tmp, nr2tmp, nr3tmp, +1)
        f_r(1:nrtmp,irp) = conjg ( aux )
     enddo

 END SUBROUTINE fft6_g2r
END MODULE gwsigma


MODULE gwsymm
       INTEGER :: ngmunique
       INTEGER, ALLOCATABLE :: ig_unique(:)
       INTEGER, ALLOCATABLE :: sym_ig(:)
       INTEGER, ALLOCATABLE :: sym_friend(:)
       LOGICAL   :: use_symm
END MODULE gwsymm


MODULE gwcom
  USE modes
  USE dynmat
  USE qpoint
  USE eqv
  USE efield_mod
  USE nlcc_gw
  USE gc_gw
  USE gwus
  USE partial
  USE control_gw
  USE freq_gw
  USE units_gw
  USE output
  USE gamma_gamma
  USE disp 
END MODULE gwcom

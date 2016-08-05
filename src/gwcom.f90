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
!
!... Common variables for the GW program
!  
MODULE eqv_gw
  USE kinds, ONLY :  DP
  USE eqv
  !
  ! ... The wavefunctions at point k+q 
  !
  SAVE
  !
  ! ... The variable describing the linear response problem 
  !
  COMPLEX (DP), ALLOCATABLE :: dpsip(:,:), dpsim(:,:)
  ! dpsip change of wavefunction for positive frequency
  ! dpsim change of wavefunction for negative frequency

  COMPLEX (DP), ALLOCATABLE :: dvbare(:)

  REAL (DP), ALLOCATABLE :: eprectot(:,:) ! needed for preconditioning
  !
END MODULE eqv_gw
!
!
MODULE efield_mod
  USE kinds, ONLY :  DP
  !
  ! ... the variables for the electric field perturbation
  !  
  SAVE
  !
  REAL (DP) :: epsilon (3,3)
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
  USE control_lr
  USE gw_type_mod, ONLY : output_type, name_length
  !
  ! ... the variables controlling the GW run
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter = 1000
  ! maximum number of iterations
  INTEGER :: niter_gw, nmix_gw, &
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
  REAL (DP) :: alpha_mix(maxter), time_now
  ! the mixing parameter
  ! CPU time up to now
  ! the alpha value for shifting the bands
  CHARACTER(LEN=10)  :: where_rec='no_recover'! where the gw run recovered
  CHARACTER(LEN=256) :: flmixdpot, tmp_dir_gw, tmp_dir_coul
  INTEGER :: rec_code, &   ! code for recover
             rec_code_read=-1000 ! code for recover. Not changed during the run

  INTEGER :: truncation ! method used to truncate in reciprocal space

  INTEGER :: maxter_coul, maxter_green, w_green_start

  LOGICAL :: lgamma_gamma,&! if .TRUE. this is a q=0 computation with k=0 only 
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
             do_sigma_exxG, &
             do_green, &
             do_sigma_matel,&
             do_q0_only,&
             prec_direct,&
             prec_shift,&
             godbyneeds,&
             padecont,&
             cohsex,&
             multishift,&
             do_sigma_extra,&
             solve_direct,&
             tinvert,&
             coul_multishift,&
             trunc_2d,&
             do_diag_w,&
             do_diag_g,&
             do_epsil,& !in case you want to set xq point
             do_imag,&    !from input file for epsilon
             do_pade_coul,&
             newgrid = .FALSE.,&
             loqua   = .FALSE.,&
             high_io, &  
             freq_gl, &
             just_corr,&
             double_grid

   TYPE(output_type) output

END MODULE control_gw
!
!
MODULE freq_gw
  !
  USE kinds,   ONLY : DP
  !
  SAVE
  ! ... the variables for computing frequency dependent dielectric constant
  INTEGER, PARAMETER :: nfsmax=350  ! # of maximum frequencies
  INTEGER :: nfs                   ! # of frequencies
  !REAL (KIND=DP) :: fiu(nfsmax)    ! values  of frequency
  COMPLEX (KIND=DP) :: fiu(nfsmax)    ! values  of frequency
  !variables for convolution
  INTEGER :: nwcoul, nwgreen, nwalloc, nwsigma, nwsigwin

  !The wsigmamin, wsigmamax, etc is currently being set in freqbins. 
  !I will change this so that it becomes a user defined option in the punchcard
  !with default values eventually. 

  REAL(DP) :: wsigmamin, wsigmamax, deltaw, wcoulmax, wgreenmin, wgreenmax 
  !Grid for the analytic continuation
  REAL(DP) :: wsig_wind_max, wsig_wind_min, deltaws
  REAL(DP), ALLOCATABLE :: wtmp(:), wcoul(:), wgreen(:), wsigma(:), wgtcoul(:) 
  INTEGER, ALLOCATABLE :: ind_w0mw (:,:), ind_w0pw (:,:)
  REAL,    ALLOCATABLE :: w0pmw (:,:)
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
END MODULE units_gw
!
!
MODULE output_mod
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=256) :: fildyn, fildvscf, fildrho, filsigc, filsigx, filcoul
  ! output file for the dynamical matrix
  ! output file for deltavscf
  ! output file for deltarho
  !
END MODULE output_mod
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
   REAL(DP) :: xk_kpoints(3, 2000)
   INTEGER  :: num_k_pts
   INTEGER  :: w_of_q_start, w_of_k_start, w_of_k_stop
  

END MODULE disp

MODULE gwsigma
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : omega, alat
  USE qpoint,      ONLY : xq, igkq
  USE fft_custom,  ONLY : fft_cus, set_custom_grid, ggent, gvec_init

  SAVE

  TYPE(fft_cus) sigma_x_st   ! Grid for \Sigma^{x} -> real space
  TYPE(fft_cus) sigma_c_st   ! Grid for real space -> restricted G space

  COMPLEX(DP), ALLOCATABLE :: sigma_band_exg(:)

! HL self energy is a huge quantity!
  INTEGER  :: nbnd_sig
! Cutoff for the sigma + exchange/correlation.
  REAL(DP) :: ecutsex
  REAL(DP) :: ecutsco
  REAL(DP) :: ecutprec
  INTEGER :: gexcut
! To easily test convergence at the end 
! of the calculation.
  REAL(DP)    :: corr_conv
  REAL(DP)    :: exch_conv
! Old FFT routines
! Real space mesh for description of self-energy.
  REAL(DP)    :: gcutmsig
  INTEGER     :: ngmsig, ngmsco, ngmsex, ngmpol, ngmgrn, gcutcorr
END MODULE gwsigma


MODULE gwsymm
       INTEGER :: ngmunique
       INTEGER, ALLOCATABLE :: ig_unique(:)
       INTEGER, ALLOCATABLE :: sym_ig(:)
       INTEGER, ALLOCATABLE :: sym_friend(:)
       LOGICAL   :: use_symm
END MODULE gwsymm


MODULE gwcom
  USE qpoint
  USE eqv_gw
  USE efield_mod
  USE nlcc_gw
  USE partial
  USE control_gw
  USE freq_gw
  USE units_gw
  USE output_mod
  USE gamma_gamma
  USE disp 
END MODULE gwcom

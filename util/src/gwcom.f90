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
!
!... Common variables for the GW program
!  
MODULE eqv_gw
  USE kinds, ONLY :  DP
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
  !
END MODULE eqv_gw
!
!
MODULE nlcc_gw
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for non-linear core correction
  !
  SAVE
  !
  LOGICAL :: nlcc_any
  ! .T. if any atom-type has nlcc
  !
END MODULE nlcc_gw
!
!
MODULE control_gw
  USE kinds, ONLY :  DP
  USE gw_type_mod, ONLY : output_type
  !
  ! ... the variables controlling the GW run
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter = 1000
  ! maximum number of iterations
  INTEGER :: niter_gw, nmix_gw
  ! maximum number of iterations (read from input)
  ! mixing type
  !
  REAL(DP) :: tr2_gw, tr2_green
  ! threshold for G and W
  !
  INTEGER lmax_gw, lmax_green
  ! lmax values for green and w
  !
  REAL(DP) :: eta
  ! small imaginary component
  !
  REAL(DP) :: alpha_mix(maxter)
  ! the mixing parameter
  !
  CHARACTER(LEN=256) :: tmp_dir_gw, tmp_dir_coul
  ! output directory
  !
  INTEGER :: truncation
  ! method used to truncate in reciprocal space
  !
  INTEGER :: model_coul
  ! screening model used in the calculation
  !
  INTEGER :: maxter_coul, maxter_green
  ! maximum iteration for G and W
  !
  LOGICAL :: convt,       &! if .TRUE. the GW has converged
             lnscf,       &! if .TRUE. the run makes first a nscf calculation
             ldisp,       &! if .TRUE. the run calculates full GW dispersion
             reduce_io,   &! if .TRUE. reduces needed I/O
             set_alpha_pv, & ! if .TRUE. use automatic method to determine alpha_pv
             do_coulomb, &
             do_sigma_c, &
             do_sigma_exx, &
             do_sigma_matel,&
             do_q0_only,&
             solve_direct,&
             coul_multishift,&
             do_epsil,& !in case you want to set xq point
             do_imag,&    !from input file for epsilon
             newgrid = .FALSE.,&
             double_grid,&
             plot_coul = .FALSE. ! plot the Coulomb potential

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
  INTEGER :: nfs                   ! # of frequencies
  COMPLEX(KIND = dp), ALLOCATABLE :: fiu(:)    ! values  of frequency
  !variables for convolution
  INTEGER :: nwcoul, nwsigma, nwsigwin

  !The wsigmamin, wsigmamax, etc is currently being set in freqbins. 
  !I will change this so that it becomes a user defined option in the punchcard
  !with default values eventually. 

  REAL(DP) :: wsigmamin, wsigmamax, wcoulmax
  !Grid for the analytic continuation
  REAL(DP) :: wsig_wind_max, wsig_wind_min

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
       iuwfc, lrwfc,      & ! wave functions
       iubar, lrbar,      & ! DV_{bare}
       iuncoul, lrcoul,   & ! screened Coulomb interaction
       iunsigma, lrsigma, & ! correlation self energy
       iunsex, lrsex        ! exchange self energy

END MODULE units_gw
!
!
MODULE output_mod
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=256) :: filsigc, filsigx, filcoul
  ! output file for correlation self energy
  ! output file for exchange self energy
  ! output file for screened Coulomb interaction
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
   REAL(DP) , ALLOCATABLE :: wq(:)
    ! weight of the q points
   INTEGER  :: num_k_pts
    ! number of k-points to be calculated
   REAL(DP) :: xk_kpoints(3, 2000)
    ! coordinates of the k points 
   INTEGER  :: w_of_q_start, w_of_k_start, w_of_k_stop
    ! impose restriction on interval of k and/or q points

END MODULE disp

MODULE gwsigma
  USE kinds,       ONLY : DP

  SAVE

  INTEGER  :: nbnd_sig
   ! bands for which expectation value is evaluated

  REAL(DP) :: ecutsex
  REAL(DP) :: ecutsco
   ! Cutoff for the sigma + exchange/correlation.
  REAL(DP)    :: corr_conv
  REAL(DP)    :: exch_conv
   ! test convergence of calculation

END MODULE gwsigma


MODULE gwsymm
       INTEGER :: ngmunique
       INTEGER, ALLOCATABLE :: ig_unique(:)
       INTEGER, ALLOCATABLE :: sym_ig(:)
       INTEGER, ALLOCATABLE :: sym_friend(:)
       LOGICAL   :: use_symm
END MODULE gwsymm


MODULE gwcom
  USE eqv_gw
  USE nlcc_gw
  USE control_gw
  USE freq_gw
  USE units_gw
  USE output_mod
  USE disp 
END MODULE gwcom

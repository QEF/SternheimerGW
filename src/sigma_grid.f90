!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2017
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
!> This module provides the subroutines to initialize the Fourier transform
!! grids used for exchange and correlation.
!!
!! Because exchange and correlation are significantly less dispersive than the
!! charge density, we can reduce the energy cutoff used to calculate these
!! quantities. This module provides the routines to initialize the Fourier
!! transform types, generate the associated Fourier transform grid and print
!! some information about these grids.
!!
MODULE sigma_grid_module

  USE fft_custom, ONLY: fft_cus

  IMPLICIT NONE

  PUBLIC sigma_grid, sigma_grid_type

  !> contains the FFT grids used for exchange and correlation
  TYPE sigma_grid_type

    !> grid used for exchange
    TYPE(fft_cus) exch

    !> grid used for correlation
    TYPE(fft_cus) corr

    !> grid used for correlation (parallelized over images)
    TYPE(fft_cus) corr_par

  END TYPE

CONTAINS

  !> Create a Fourier transform grid with a different energy cutoff.
  SUBROUTINE sigma_grid_create(comm, gamma_only, tpiba2, ecut_cust, fft_cust)

    USE kinds,      ONLY: dp
    USE fft_custom, ONLY: fft_cus, ggent, gvec_init

    !> communicator over which the routines are parallelized
    INTEGER,  INTENT(IN) :: comm

    !> Special case for \f$\Gamma\f$-only calculation
    LOGICAL,  INTENT(IN) :: gamma_only

    !> \f$(2 \pi / a)^2\f$ of the lattice
    REAL(dp), INTENT(IN) :: tpiba2

    !> The energy cutoff used for the custom type.
    REAL(dp), INTENT(IN) :: ecut_cust

    !> The custom FFT type initialized by this routine
    TYPE(fft_cus), INTENT(OUT) :: fft_cust

    !> number of g vectors on this process
    INTEGER num_g

    !!
    !! 1. set the energy cutoff of this type
    !!
    fft_cust%ecutt = ecut_cust

    !!
    !! 2. converts the energy cutoff to a cutoff for the G vectors
    !!
    fft_cust%gcutmt = fft_cust%ecutt / tpiba2

    !!
    !! 4. initialize the fft type
    !!
    CALL wrapper_fft_type_init(comm, gamma_only, fft_cust)

    !!
    !! 5. initialize the number of G vectors and allocate arrays
    !!
    num_g = fft_cust%dfftt%ngl(fft_cust%dfftt%mype + 1)
    IF (gamma_only) num_g = (num_g + 1) / 2
    CALL gvec_init(fft_cust, num_g, comm)

    !!
    !! 6. generate the FFT grid
    !!
    CALL ggent(fft_cust)
    fft_cust%initialized = .TRUE.

  END SUBROUTINE sigma_grid_create

  !> This routine provides a wrapper to interface the fft_type_init routine with
  !! the necessary quantities extracted from QE's global modules.
  SUBROUTINE wrapper_fft_type_init(comm, gamma_only, fft_cust)

    USE cell_base,  ONLY: at, bg
    USE fft_custom, ONLY: fft_cus
    USE fft_types,  ONLY: fft_type_init
    USE mp_bands,   ONLY: nyfft
    USE stick_base, ONLY: sticks_map

    !> The communicator over which we parallelize
    INTEGER,       INTENT(IN)    :: comm

    !> Special case for \f$\Gamma\f$-only calculation
    LOGICAL,       INTENT(IN)    :: gamma_only

    !> The custom FFT type initialized by this routine
    TYPE(fft_cus), INTENT(INOUT) :: fft_cust

    !> the sticks map created by fft type init
    TYPE(sticks_map) smap
 
#if defined(__MPI) && ! defined(__USE_3D_FFT)
    LOGICAL, PARAMETER :: lpara = .TRUE.
#else
    LOGICAL, PARAMETER :: lpara = .FALSE.
#endif

    CALL fft_type_init(fft_cust%dfftt, smap, "rho", gamma_only, lpara, comm, at, bg, fft_cust%gcutmt, nyfft=nyfft)

  END SUBROUTINE wrapper_fft_type_init
  
  !> Print info on local and global dimensions for real space grids
  SUBROUTINE sigma_grid_info(fft_cust, label)
  
    USE fft_custom, ONLY: fft_cus
    USE fft_types,  ONLY: fft_type_descriptor
    USE io_global,  ONLY: stdout
  
    IMPLICIT NONE
  
    !> The type we want to print info about.
    TYPE(fft_cus), INTENT(IN) :: fft_cust
  
    !> The label of the type.
    CHARACTER(*), INTENT(IN)  :: label
  
    !> The description added to the label
    CHARACTER(*), PARAMETER   :: descr = 'Grid'
  
    !> The length of the label
    INTEGER len_label
  
    !> counter on the number of processes
    INTEGER iproc
  
    !> abbreviation for fft_cust%dfftt
    TYPE(fft_type_descriptor) dfft

    dfft = fft_cust%dfftt
    !
    ! add 1 for the space between label and description
    len_label = LEN_TRIM(label) + LEN_TRIM(descr) + 1
    !
    WRITE(stdout,*)
    WRITE(stdout,'(5x,a,1x,a)') TRIM(label), TRIM(descr)
    WRITE(stdout,'(5x,a)') REPEAT('-', len_label)
    WRITE(stdout,'(5x,a)') 'E_cutoff(Ry) G_cutoff(crystal) num G vec'
    WRITE(stdout,'(5x,f8.2,7x,f8.2,7x,i6)') fft_cust%ecutt, fft_cust%gcutmt, fft_cust%ngmt
    WRITE(stdout,'(5x,a)') 'nr1t  nr2t  nr3t'
    WRITE(stdout,'(5x,3(1x,i5))') fft_cust%nr1t, fft_cust%nr2t, fft_cust%nr3t
    WRITE(stdout,'(5x,a)') 'Global Dimensions   Local  Dimensions   Processor Grid'
    WRITE(stdout,'(5x,a)') '.X.   .Y.   .Z.     .X.   .Y.   .Z.     .X.   .Y.   .Z.'
    WRITE(stdout,'(2x,3(1x,i5),2x,3(1x,i5),2x,3(1x,i5))') &
      dfft%nr1, dfft%nr2, dfft%nr3, dfft%nr1, dfft%my_nr2p, dfft%my_nr3p, 1, dfft%nproc2, dfft%nproc3
    WRITE(stdout,'(5x,a,3(1x,i5))') 'Array leading dimensions ( nr1x, nr2x, nr3x )   = ', &
      dfft%nr1x, dfft%nr2x, dfft%nr3x
    WRITE(stdout,'(5x,a,1x,i9)') 'Local number of cell to store the grid ( nrxx ) = ', &
      dfft%nnr
    WRITE(stdout,'(5x,a)') 'size of the "Z" section of each processor'
    ! print local part for 10 processors at once
    DO iproc = 1, dfft%nproc3, 10
      WRITE(stdout,'(5x,a,10i5)') 'nr3p = ', dfft%nr3p(iproc : MIN(dfft%nproc3, iproc + 10))
    END DO ! iproc
    !
  END SUBROUTINE sigma_grid_info

  !> This subroutine generates the FFT grids for Exchange and Correlation.
  !!
  !! At the moment the necessary quantities are taken from the global modules.
  !! We call a grid generation routine for exchange and correlation. It prints
  !! some information about the generated grids and estimates the necessary
  !! memory.
  !!
  SUBROUTINE sigma_grid(freq, ecut_x, ecut_c, grid)
    !
    USE cell_base,        ONLY : tpiba2
    USE control_flags,    ONLY : gamma_only
    USE freqbins_module,  ONLY : freqbins_type
    USE io_global,        ONLY : stdout, ionode
    USE kinds,            ONLY : dp
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp_images,        ONLY : inter_image_comm, nimage
    ! 
    IMPLICIT NONE
    !
    !> type that defines the frequencies used
    TYPE(freqbins_type),   INTENT(IN)  :: freq
    !
    !> the energy cutoff used for exchange
    REAL(dp),              INTENT(IN)  :: ecut_x
    !
    !> the energy cutoff used for correlation
    REAL(dp),              INTENT(IN)  :: ecut_c
    !
    !> the FFT grid used in the code
    TYPE(sigma_grid_type), INTENT(OUT) :: grid
    !
    !> number of bytes in a complex
    REAL(dp), PARAMETER :: complex_byte = 16.0_dp
    !> number of bytes in a mega byte
    REAL(dp), PARAMETER :: mega_byte = 10e6_dp
    !> convert a number of complex to a size in MB
    REAL(dp), PARAMETER :: in_MB = complex_byte / mega_byte
    !
    !> format to write the memory information
    CHARACTER(*), PARAMETER :: myfmt = '(5x,a,2x,f9.1,1x,a)'
    !
    !> helper variables for the number of G in exchange and correlation grid
    REAL(dp) num_g_x, num_g_c, num_g_ci
    !> helper variables for the number of r in exchange and correlation grid
    REAL(dp) num_r_x, num_r_c, num_r_ci
    !> helper product for the correlation grid
    REAL(dp) num_g_c_pr, num_r_c_pr

    !> number of frequencies used in the linear solver
    INTEGER num_freq

    !> number of frequencies used for the Green's function
    INTEGER num_green

    !> number of frequencies used for the self energy integration
    INTEGER num_sigma

    num_freq  = freq%num_freq()
    num_green = 2 * freq%num_coul()
    num_sigma = freq%num_sigma()

    !
    ! Generate the exchange grid
    !
    CALL sigma_grid_create(intra_bgrp_comm, gamma_only, tpiba2, ecut_x, grid%exch)
    IF (ionode) CALL sigma_grid_info(grid%exch, 'Exchange')
  
    !
    ! Generate the correlation grid
    !
    CALL sigma_grid_create(intra_bgrp_comm, gamma_only, tpiba2, ecut_c, grid%corr)
    IF (ionode) CALL sigma_grid_info(grid%corr, 'Correlation')
    !
    IF (nimage == 1) THEN
      ! reuse the same grid if only 1 image is used
      grid%corr_par = grid%corr
    ELSE
      ! create a grid parallelized over images
      CALL sigma_grid_create(inter_image_comm, gamma_only, tpiba2, ecut_c, grid%corr_par)
      IF (ionode) CALL sigma_grid_info(grid%corr_par, 'Correlation (images)')
    END IF

    !
    ! Print info about array size
    !
    num_g_x  = REAL(grid%exch%ngmt, KIND=dp)
    num_g_c  = REAL(grid%corr%ngmt, KIND=dp)
    num_g_ci = REAL(grid%corr_par%ngmt, KIND=dp)
    num_r_x  = REAL(grid%exch%dfftt%nnr, KIND=dp)
    num_r_c  = REAL(grid%corr%dfftt%nnr, KIND=dp)
    num_r_ci = REAL(grid%corr_par%dfftt%nnr, KIND=dp)
    !
    ! evaluate product of grid size
    num_g_c_pr = num_g_c * num_g_ci
    num_r_c_pr = num_r_c * num_r_ci
    !
    WRITE(stdout,*)
    WRITE(stdout,'(5x,a)') 'Memory Usage:'
    WRITE(stdout,*)
    WRITE(stdout, myfmt) "G(G, G'; w):    ", num_g_c_pr * REAL(num_green, KIND=dp) * in_MB, 'MB'
    WRITE(stdout, myfmt) "G(r, r'; w):    ", num_r_c_pr * REAL(num_green, KIND=dp) * in_MB, 'MB'
    WRITE(stdout, myfmt) "W(G, G'; w):    ", num_g_c**2 * REAL(num_freq,  KIND=dp) * in_MB, 'MB'
    WRITE(stdout, myfmt) "W(r, r'):       ", num_r_c_pr                            * in_MB, 'MB'
    WRITE(stdout, myfmt) "Sigma(G, G'; w):", num_g_c_pr * REAL(num_sigma, KIND=dp) * in_MB, 'MB'
    WRITE(stdout, myfmt) "Sigma(r, r'; w):", num_r_c_pr * REAL(num_sigma, KIND=dp) * in_MB, 'MB'
    WRITE(stdout, myfmt) "V(G, G'):       ", num_g_x**2                            * in_MB, 'MB'
    WRITE(stdout, myfmt) "V(r, r'):       ", num_r_x**2                            * in_MB, 'MB'
    WRITE(stdout,*)

  END SUBROUTINE sigma_grid

END MODULE sigma_grid_module

!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
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
!> This module provides the routines to evaluate the exchange self-energy.
!!
!! The exchange self energy is defined as
!! \f{equation}{
!!   \Sigma_k^{\text{x}}(G, G') = -\sum_{nq} w_{nk-q} \langle k + G, n k - q \lvert 
!!   V \rvert n k - q, k + G'\rangle~,
!! \f}
!! where \f$\lvert k + G\rangle\f$ is a plane wave, \f$\lvert n k\rangle\f$
!! are the eigenstates of the Hamiltonian
!! \f{equation}{
!!   H \lvert n k \rangle = \epsilon_{nk} \lvert n k \rangle~,
!! \f}
!! \f$w_{nk}\f$ is weight times occupation of state \f$\lvert n k\rangle\f$, and \f$V\f$
!! is the (truncated) Coulomb potential.
!!
!! To evaluate this equation, we expand the eigenstates in the plane wave basis
!! \f{equation}{
!!   \lvert n k \rangle = \sum_{G} c_{nk}(G) \lvert k + G \rangle~.
!! \f}
!! This leads to the following expression for the self-energy
!! \f{equation}{
!!   \Sigma_k^{\text{x}}(G, G') = -\sum_{nq} w_{nk-q} \sum_{G_1 G_2} c^\ast_{nk-q}(G_1) c_{nk-q}(G_2)
!!   \langle k + G, k - q + G_1 \lvert V \rvert k - q + G_2, k + G'\rangle~.
!! \f}
!! We can evaluate the exchange integral in the basis of plane waves
!! \f{equation}{
!!   \langle k + G, k - q + G_1 \lvert V \rvert k - q + G_2, k + G'\rangle
!!   = \frac{1}{\Omega} V_q(G - G_2) \delta(G_1 - [G_2 + G' - G])~.
!! \f}
!! Here, \f$V_q(G)\f$ is the Fourier transform of the (truncated) Coulomb potential.
!! Introducing \f$G'' = G - G_2\f$ and resolving the summation over \f$G_1\f$ with the
!! \f$\delta\f$ function yields
!! \f{equation}{
!!   \Sigma_k^{\text{x}}(G, G') = -\frac{1}{\Omega} \sum_{nq} w_{nk-q} \sum_{G''}
!!     c^\ast_{nk-q}(G' - G'') c_{nk-q}(G - G'') V_q(G'')~.
!! \f}
MODULE exchange_module

  IMPLICIT NONE

  PRIVATE

  PUBLIC exchange_wrapper

  !> a number indicating that the map is out of the boundary of the array
  INTEGER, PARAMETER :: out_of_bound = 0

CONTAINS

  !> Evaluate the convolution over \f$G''\f$.
  !!
  !! The individual parts of the exchange self-energy are given as
  !! \f{equation}{
  !!   \Sigma_{k,q}^{\text{x}} = w_{nk-q} \sum_{G''}V_q(G'') 
  !!   c_{nk-q}^\ast(G' - G'') c_{nk-q}(G - G'')~.
  !! \f}
  !! This routine adds the current element to the previously calculated ones,
  !! so that the array contains the sum over all \w$(nq)\f$ in the end.
  SUBROUTINE exchange_convolution(occupation, coulomb, evec, map, sigma)

    USE kinds, ONLY: dp

    !> The occupation \f$w_{nk-q}\f$ of the state \f$\lvert nk-q\rangle\f$.
    REAL(dp),    INTENT(IN)    :: occupation

    !> The coulomb potential \f$V_q(G'')\f$.
    REAL(dp),    INTENT(IN)    :: coulomb(:)

    !> The eigenvector \f$c_{nk-q}(G)\f$.
    COMPLEX(dp), INTENT(IN)    :: evec(:)

    !> A map of two G indices onto the index of their difference, i.e.,
    !! if \f$G_i - G_j = G_k\f$ then \f$\text{map}(i,j) = k\f$. If an
    !! element is not a valid index of evec, set it to out_of_bound.
    INTEGER,     INTENT(IN)    :: map(:,:)

    !> The self-energy containing the sum of all \f$(nq)\f$ parts.
    COMPLEX(dp), INTENT(INOUT) :: sigma(:,:)

    !> the number of G points
    INTEGER num_g

    !> counter on the \f$G,~G'~\text{and}~G''\f$
    INTEGER ig, igp, igpp

    !> pointer to \f$G - G''\f$ and \f$G' - G''\f$
    INTEGER g_gpp, gp_gpp

    !
    ! sanity test
    !
    ! num_g contains number of G
    num_g = SIZE(coulomb)
    IF (ANY(map > SIZE(evec))) &
      CALL errore(__FILE__, "map points out of the bounds of the array", 1)
    IF (SIZE(sigma, 1) /= num_g) &
      CALL errore(__FILE__, "1st dim of self-energy not consistent with eigenvector", 1)
    IF (SIZE(sigma, 2) /= num_g) &
      CALL errore(__FILE__, "2nd dim of self-energy not consistent with eigenvector", 1)

    !
    ! summation over G''
    !
    DO igpp = 1, num_g
      !
      ! loop over G'
      DO igp = 1, num_g
        !
        ! skip element if G' - G'' is out of bounds
        gp_gpp = map(igp, igpp)
        IF (gp_gpp == out_of_bound) CYCLE
        !
        ! loop over G
        DO ig = 1, num_g
          !
          ! skip element if G - G'' is out of bounds
          g_gpp = map(ig, igpp)
          IF (g_gpp == out_of_bound) CYCLE
          !
          ! sigma -= w_{nk-q} V_q(G'') c*_nk-q(G' - G'') c_nk-q(G - G'')
          sigma(ig, igp) = sigma(ig, igp) - occupation * coulomb(igpp) &
                                          * CONJG(evec(gp_gpp)) * evec(g_gpp)
          !
        END DO ! ig
      END DO ! igp
    END DO ! igpp

  END SUBROUTINE exchange_convolution

  !> Construct a map from G and G' to G - G'.
  SUBROUTINE exchange_map(unit_cell, grid, mill, index_g, map)

    USE kinds,             ONLY: dp
    USE sigma_grid_module, ONLY: sigma_grid_type

    !> The unit cell of the crystal
    REAL(dp),             INTENT(IN)  :: unit_cell(3,3)

    !> The grid used to calculate the exchange
    TYPE(sigma_grid_type),INTENT(IN)  :: grid

    !> The global list of G vectors
    INTEGER,              INTENT(IN)  :: mill(:,:)

    !> The index of a G vector in the wave function
    INTEGER,              INTENT(IN)  :: index_g(:)

    !> The map from G and G' to G - G'
    INTEGER, ALLOCATABLE, INTENT(OUT) :: map(:,:)

    !> Flag that indicates that we transform to crystal units.
    INTEGER,    PARAMETER :: to_crystal = -1

    !> The G vectors in crystal coordinates (float)
    REAL(dp), ALLOCATABLE :: gvec_r(:,:)

    !> The G vectors in crystal coordinates (integer)
    INTEGER,  ALLOCATABLE :: gvec_i(:,:)

    !> the index of a certain G vector
    INTEGER,  ALLOCATABLE :: g_index(:,:,:)

    !> the maximum value in the gvec array
    INTEGER max_g(3)

    !> the difference of G and G'
    INTEGER delta_g(3)

    !> counter on the G and G' vector
    INTEGER ig, igp

    !
    ! transform the grid to crystal coordinates
    !
    ALLOCATE(gvec_r(3, grid%exch_fft%ngm))
    gvec_r = grid%exch_gvec
    CALL cryst_to_cart(grid%exch_fft%ngm, gvec_r, unit_cell, to_crystal)
    !
    ALLOCATE(gvec_i(3, grid%exch_fft%ngm))
    gvec_i = NINT(gvec_r)
    DEALLOCATE(gvec_r)

    !
    ! create the map for G and G' to G - G'
    !
    ! determine size of helper array
    DO ig = 1, 3
      max_g(ig) = MAXVAL(ABS(mill(ig,:)))
    END DO ! ig
    !
    ! create helper array
    ALLOCATE(g_index(-max_g(1):max_g(1), -max_g(2):max_g(2), -max_g(3):max_g(3)))
    g_index = out_of_bound
    !
    ! create map from vector to index
    DO ig = 1, SIZE(index_g)
      !
      ! skip last elements not present on current k point
      IF (index_g(ig) == out_of_bound) EXIT
      !
      ! distance of the origin
      delta_g = mill(:, index_g(ig))
      !
      ! set pointer to wave function index
      g_index(delta_g(1), delta_g(2), delta_g(3)) = ig
      !
    END DO ! ig
    !
    ! now create the output map
    ALLOCATE(map(grid%exch_fft%ngm, grid%exch_fft%ngm))
    map = out_of_bound
    !
    DO igp = 1, grid%exch_fft%ngm
      !
      DO ig = 1, grid%exch_fft%ngm
        !
        ! evaluate vector G - G' and test if it is within array boundaries
        delta_g = gvec_i(:,ig) - gvec_i(:,igp)
        IF (ANY(ABS(delta_g) > max_g)) CYCLE
        !
        ! find the index corresponding to G - G'
        map(ig, igp) = g_index(delta_g(1), delta_g(2), delta_g(3))
        !
      END DO ! ig
      !
    END DO ! igp

  END SUBROUTINE exchange_map

  !> Evaluate the Coulomb potential for all vectors within the exchange grid.
  SUBROUTINE exchange_coulomb(tpiba, method, vcut, grid, qvec, coulomb)

    USE kinds,               ONLY: dp
    USE sigma_grid_module,   ONLY: sigma_grid_type
    USE truncation_module,   ONLY: truncate, vcut_type

    !> \f$2 \pi / a\f$ where \f$a\f$ is the first dimension of the lattice
    REAL(dp), INTENT(IN) :: tpiba

    !> The ID of the truncation method used
    INTEGER,  INTENT(IN) :: method

    !> The truncated Coulomb potential
    TYPE(vcut_type), INTENT(IN) :: vcut

    !> The grid used to calculate the exchange
    TYPE(sigma_grid_type), INTENT(IN) :: grid

    !> The q-point at which the Coulomb potential is evaluated
    REAL(dp), INTENT(IN) :: qvec(3)

    !> The truncated Coulomb potential for all G vectors in the grid
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: coulomb(:)

    !> the vector q + G
    REAL(dp) q_G(3)

    !> counter on the G vectors
    INTEGER ig

    !
    ! allocate arrays of the appropriate size
    !
    ALLOCATE(coulomb(grid%exch_fft%ngm))

    !
    ! generate the Coulomb potential for all q + G
    !
    DO ig = 1, grid%exch_fft%ngm
      !
      q_G = (qvec + grid%exch_gvec(:,ig)) * tpiba
      !
      ! evaluate the truncated coulomb potential
      coulomb(ig) = truncate(method, vcut, q_G)
      !
    END DO ! ig

  END SUBROUTINE exchange_coulomb

  !> Extract the necessary quantities and evaluate the exchange according
  !! to the following algorithm.
  SUBROUTINE exchange_wrapper(ikpt, grid, vcut)

    USE buffers,            ONLY: get_buffer
    USE cell_base,          ONLY: tpiba, at, omega
    USE constants,          ONLY: degspin
    USE control_gw,         ONLY: output, truncation
    USE control_lr,         ONLY: nbnd_occ 
    USE disp,               ONLY: xk_kpoints
    USE eqv,                ONLY: evq
    USE gvect,              ONLY: mill
    USE io_global,          ONLY: meta_ionode
    USE kinds,              ONLY: dp
    USE klist,              ONLY: xk, igk_k, ngk
    USE mp_images,          ONLY: inter_image_comm, root_image
    USE mp_pools,           ONLY: inter_pool_comm, root_pool
    USE parallel_module,    ONLY: parallel_task, mp_root_sum
    USE qpoint,             ONLY: nksq, ikqs, npwq
    USE sigma_grid_module,  ONLY: sigma_grid_type
    USE sigma_io_module,    ONLY: sigma_io_write_x
    USE units_gw,           ONLY: iunsex, lrsex, iuwfc, lrwfc
    USE wvfct,              ONLY: wg
    USE timing_module,      ONLY: time_sigma_x
    USE truncation_module,  ONLY: vcut_type

    !> The index of the k-point for which the exchange is evaluated
    INTEGER, INTENT(IN) :: ikpt

    !> the FFT grid for exchange
    TYPE(sigma_grid_type), INTENT(IN) :: grid

    !> The truncated Coulomb potential
    TYPE(vcut_type), INTENT(IN) :: vcut

    !> temporary array to distribute the work
    INTEGER, ALLOCATABLE :: num_task(:)

    !> counter on the q points
    INTEGER iq

    !> the first and last band on this process
    INTEGER iband_start, iband_stop

    !> counter on the bands
    INTEGER iband

    !> index of the wave function k - q
    INTEGER ikq

    !> the q point at which the Coulomb potential is evaluated
    REAL(dp) qvec(3)

    !> a map from two G indices on the index of their difference
    INTEGER,     ALLOCATABLE :: map(:,:)

    !> the truncated Coulomb potential
    REAL(dp),    ALLOCATABLE :: coulomb(:)

    !> the exchange self-energy
    COMPLEX(dp), ALLOCATABLE :: sigma(:,:)

    !> status of allocation
    INTEGER ierr

    !> complex constant of 0
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    CALL start_clock(time_sigma_x)

    ! allocate array for self energy and initialize to 0
    ALLOCATE(sigma(grid%exch_fft%ngm, grid%exch_fft%ngm), STAT=ierr)
    CALL errore(__FILE__, "error allocating array for exchange self energy", ierr)
    sigma = zero

    !!
    !! 1. Extract the wave function of every q-point
    !!
    DO iq = 1, nksq
      !
      ikq  = ikqs(iq)
      npwq = ngk(ikq)
      !
      CALL get_buffer(evq, lrwfc, iuwfc, ikq)
      evq(npwq + 1:, :) = zero
      !!
      !! 2. Distribute the work over the process grid
      !!
      CALL parallel_task(inter_image_comm, nbnd_occ(ikq), iband_start, iband_stop, num_task)
      !!
      !! 3. construct the map from G and G' to G - G'
      !!
      CALL exchange_map(at, grid, mill, igk_k(:,ikq), map)
      !!
      !! 4. construct the Coulomb potential
      !!
      ! q = k - (k - q)
      qvec = xk_kpoints(:, ikpt) - xk(:, ikq)
      CALL exchange_coulomb(tpiba, truncation, vcut, grid, qvec, coulomb)
      !!
      !! 5. every process evaluates his contribution to sigma
      !!
      DO iband = iband_start, iband_stop
        !
        CALL exchange_convolution(wg(iband, ikq) / degspin, coulomb, evq(:,iband), map, sigma)
        !
      END DO ! iband
      !
    END DO ! iq

    !!
    !! 6. collect the self-energy on the root process
    !!
    CALL mp_root_sum(inter_image_comm, root_image, sigma)
    CALL mp_root_sum(inter_pool_comm,  root_pool,  sigma)

    !!
    !! 7. write the self-energy to file
    !!
    IF (meta_ionode) THEN
      !
      sigma = sigma / omega
      ! write unformatted
      CALL davcio(sigma, lrsex, iunsex, ikpt, 1)
      ! write formatted
      CALL sigma_io_write_x(output%unit_sigma, ikpt, sigma)
      !
    END IF ! meta_ionode

    CALL stop_clock(time_sigma_x)

  END SUBROUTINE exchange_wrapper

END MODULE exchange_module

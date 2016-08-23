!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
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
!> This module provides the routines to evaluate the exchange self-energy.
!!
!! In real space the exchange energy is defined as
!! \f{equation}{
!!   \Sigma_{k}^{\text{x}}(r, r') = -\sum_{nq} f_{nk-q} 
!!   \frac{\psi_{nk-q}^\ast(r) \psi_{nk-q}(r')}{|r - r'|}~.
!! \f}
!! Here, the \f$f_{nk}\f$ and \f$\psi_{nk}\f$ are the occupation and the
!! eigenfunction for band \f$n\f$ and wave vector \f$k\f$, respectively.
!! We are interested in the quantity in reciprocal space, which is related to
!! the real space self-energy in the following way
!! \f{equation}{
!!   \Sigma_{k}^{\text{x}}(G, G') = \Omega \sum_{r,r'} e^{i (k + G) r}
!!   \Sigma_{k}^{\text{x}}(r, r') e^{-i (k + G') r'}~.
!! \f}
!! To evaluate this quantity, we conveniently evaluate the Coulomb potential
!! in reciprocal space
!! \f{equation}{
!!   V_k(G) = \sum_{r} \frac{e^{i (k + G) r}}{r} = \frac{4\pi}{|k + G|^2}~.
!! \f}
!! For practical calculations, we need to truncate the integration, which is
!! handled in the truncation_module.
!!
!! We remind of the fact that the wave function contains a Bloch phase and
!! a lattice periodic part
!! \f{equation}{
!!   \psi_{nk}(r) = e^{ikr} u_{nk}(r)
!! \f}
!! Replacing the Coulomb potential by the Fourier transform resolves the
!! convolution of \f$r\f$ and \f$r'\f$. The product of exponential and lattice
!! periodic part of the wave function corresponds to the Fourier transform of
!! the lattice periodic part. We obtain the following single \f$G''\f$ summation
!! \f{equation}{
!!   \Sigma_{k}^{\text{x}}(G, G') = -\sum_{nq} f_{nk-q} \sum_{G''}
!!   V_q(G'') u_{nk-q}^\ast(G - G'') u_{nk-q}(G'' - G')~.
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
  !!   \Sigma_{k,q}^{\text{x}} = f_{nk-q} \sum_{G''}V_q(G'') 
  !!   u_{nk-q}^\ast(G - G'') u_{nk-q}(G'' - G')~.
  !! \f}
  !! This routine adds the current element to the previously calculated ones,
  !! so that the array contains the sum over all \f$(nq)\f$ in the end.
  SUBROUTINE exchange_convolution(occupation, coulomb, evec, map, sigma)

    USE kinds, ONLY: dp

    !> The occupation \f$f_{nk-q}\f$ of the state \f$(nk-q)\f$.
    REAL(dp),    INTENT(IN)    :: occupation

    !> The coulomb potential \f$V_q(G'')\f$.
    REAL(dp),    INTENT(IN)    :: coulomb(:)

    !> The eigenvector \f$u_{nk-q}(G)\f$.
    COMPLEX(dp), INTENT(IN)    :: evec(:)

    !> A map of two G indices onto the index of their difference, i.e.,
    !! if \f$G_i - G_j = G_k\f$ then \f$\text{map}(i,j) = k\f$. If an
    !! element is not a valid index of evec, set it to out_of_bound.
    INTEGER,     INTENT(IN)    :: map(:,:)

    !> The self-energy containing the sum of all \f$(nq)\f$ parts.
    COMPLEX(dp), INTENT(INOUT) :: sigma(:,:)

    !> the number of G points
    INTEGER num_g

    !> counter on the G, G' and G''
    INTEGER ig, igp, igpp

    !> pointer to G - G'' and G'' - G'
    INTEGER g_gpp, gpp_gp

    !
    ! sanity test
    !
    ! num_g contains number of G
    num_g = SIZE(evec)
    IF (ANY(map > num_g)) &
      CALL errore(__FILE__, "map points out of the bounds of the array", 1)
    IF (SIZE(sigma, 1) /= num_g) &
      CALL errore(__FILE__, "1st dim of self-energy not consistent with eigenvector", 1)
    IF (SIZE(sigma, 2) /= num_g) &
      CALL errore(__FILE__, "2nd dim of self-energy not consistent with eigenvector", 1)

    !
    ! summation over G''
    !
    DO igpp = 1, SIZE(coulomb)
      !
      ! loop over G'
      DO igp = 1, num_g
        !
        ! skip element if G'' - G' is out of bounds
        gpp_gp = map(igpp, igp)
        IF (gpp_gp == out_of_bound) CYCLE
        !
        ! loop over G
        DO ig = 1, num_g
          !
          ! skip element if G - G'' is out of bounds
          g_gpp = map(ig, igpp)
          IF (g_gpp == out_of_bound) CYCLE
          !
          ! sigma += f_{nk-q} V_q(G'') u*_nk-q(G - G'') u_nk-q(G'' - G')
          sigma(ig, igp) = sigma(ig, igp) + occupation * coulomb(igpp) &
                                          * CONJG(evec(g_gpp)) * evec(gpp_gp)
          !
        END DO ! ig
      END DO ! igp
    END DO ! igpp

  END SUBROUTINE exchange_convolution

  !> Construct a map from G and G' to G - G'.
  SUBROUTINE exchange_map(gvec, index_g_evec, map)

    !> The list of G vectors used
    INTEGER,  INTENT(IN) :: gvec(:,:)

    !> The indices of the G vectors used for the wave function
    INTEGER,  INTENT(IN) :: index_g_evec(:)

    !> The map from G and G' to G - G'
    INTEGER,  ALLOCATABLE, INTENT(OUT) :: map(:,:)

    ! the maximum value in the gvec array
    INTEGER max_g(3)

    ! the index of a certain G vector
    INTEGER, ALLOCATABLE :: g_index(:,:,:)

    ! the number of G vectors
    INTEGER num_g

    ! counter on the G and G' vector
    INTEGER ig, igp

    ! the difference of G and G'
    INTEGER delta_g(3)

    !
    ! create the map for G and G' to G - G'
    !
    ! determine size of helper array
    DO ig = 1, 3
      max_g(ig) = MAXVAL(ABS(gvec(ig,:)))
    END DO ! ig
    !
    ! create helper array
    ALLOCATE(g_index(-max_g(1):max_g(1), -max_g(2):max_g(2), -max_g(3):max_g(3)))
    g_index = out_of_bound
    !
    ! create map from vector to index
    DO ig = 1, SIZE(gvec, 2)
      g_index(gvec(:,1), gvec(:,2), gvec(:,3)) = ig
    END DO ! ig
    !
    ! now create the output map
    num_g = SIZE(index_g_evec)
    ALLOCATE(map(num_g, num_g))
    map = out_of_bound
    !
    DO igp = 1, num_g
      !
      ! skip elements that are out of bounds
      IF (index_g_evec(igp) == out_of_bound) CYCLE
      !
      DO ig = 1, num_g
        !
        ! skip elements that are out of bounds
        IF (index_g_evec(ig) == out_of_bound) CYCLE
        !
        ! evaluate vector G - G' and test if it is within array boundaries
        delta_g = gvec(:,ig) - gvec(:,igp)
        IF (ANY(ABS(delta_g) > max_g)) CYCLE
        !
        ! find the index corresponding to G - G'
        map(ig, igp) = g_index(delta_g(1), delta_g(2), delta_g(3))
        !
      END DO ! ig
    END DO ! igp

  END SUBROUTINE exchange_map

  !> Evaluate the Coulomb potential for all vectors within the exchange grid.
  SUBROUTINE exchange_coulomb(tpiba, method, exchange_grid, qvec, gvec, coulomb, index_g_coul)

    USE fft_custom,        ONLY: fft_cus
    USE kinds,             ONLY: dp
    USE truncation_module, ONLY: truncate

    !> \f$2 \pi / a\f$ where \f$a\f$ is the first dimension of the lattice
    REAL(dp), INTENT(IN) :: tpiba

    !> The ID of the truncation method used
    INTEGER,  INTENT(IN) :: method

    !> The grid used to calculate the exchange defining the cutoff \f$q_{\text{cut}}\f$
    TYPE(fft_cus), INTENT(IN) :: exchange_grid

    !> The q-point at which the Coulomb potential is evaluated
    REAL(dp), INTENT(IN) :: qvec(3)

    !> The list of G vectors used
    REAL(dp), INTENT(IN) :: gvec(:,:)

    !> The truncated Coulomb potential for all \f$|q + G| < q_{\text{cut}}\f$
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: coulomb(:)

    !> The indices of the G vectors used for the Coulomb potential
    INTEGER,  ALLOCATABLE, INTENT(OUT) :: index_g_coul(:)

    !> the vector q + G
    REAL(dp) q_G(3)

    !> counter on the G vectors
    INTEGER ig

    !> counter on the found G vectors
    INTEGER ifound

    !
    ! allocate arrays of the appropriate size
    !
    ALLOCATE(coulomb(exchange_grid%ngmt))
    ALLOCATE(index_g_coul(exchange_grid%ngmt))
    index_g_coul = out_of_bound

    !
    ! find all |q + G| < q_cut and determine the Coulomb potential for it
    !
    ifound = 0
    DO ig = 1, SIZE(gvec, 2)
      !
      q_G = qvec + gvec(:,ig)
      !
      ! skip points outside of q_cut
      IF (SQRT(SUM(q_G**2)) >= exchange_grid%gcutmt) CYCLE
      !
      ! store the new point
      ifound = ifound + 1
      index_g_coul(ifound) = ig
      !
      ! evaluate the truncated coulomb potential
      coulomb(ifound) = truncate(method, q_G * tpiba)
      !
    END DO ! ig

  END SUBROUTINE exchange_coulomb

  !> Extract the necessary quantities and evaluate the exchange according
  !! to the following algorithm.
  SUBROUTINE exchange_wrapper(ikpt)

    USE cell_base,          ONLY: tpiba
    USE control_gw,         ONLY: output, nbnd_occ, truncation
    USE eqv_gw,             ONLY: evq
    USE gvect,              ONLY: mill, g
    USE gwsigma,            ONLY: sigma_x_st
    USE io_global,          ONLY: meta_ionode
    USE kinds,              ONLY: dp
    USE klist,              ONLY: xk, igk_k
    USE mp_images,          ONLY: inter_image_comm, root_image
    USE mp_pools,           ONLY: inter_pool_comm, root_pool
    USE qpoint,             ONLY: nksq, ikks, ikqs
    USE units_gw,           ONLY: iunsex, lrsex, iuwfc, lrwfc
    USE timing_module,      ONLY: time_sigma_x
    USE wvfct,              ONLY: wg

    !> The index of the k-point for which the exchange is evaluated
    INTEGER, INTENT(IN) :: ikpt

    !> temporary array to distribute the work
    INTEGER, ALLOCATABLE :: num_task(:)

    !> the first and last q-point on this process
    INTEGER iq_start, iq_stop

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

    !> the indices of the G vectors at q
    INTEGER,     ALLOCATABLE :: index_g_coul(:)

    !> the truncated Coulomb potential
    REAL(dp),    ALLOCATABLE :: coulomb(:)

    !> the exchange self-energy
    COMPLEX(dp), ALLOCATABLE :: sigma(:,:)

    CALL start_clock(time_sigma_x)

    !!
    !! 1. Distribute the work over the process grid
    !!
    CALL parallel_task(inter_pool_comm, nksq, iq_start, iq_stop, num_task)
    CALL parallel_task(inter_image_comm, nbnd_occ(ikpt), iband_start, iband_stop, num_task)

    !!
    !! 2. Extract the wave function of every q-point
    !!
    DO iq = iq_start, iq_stop
      !
      ! sanity check - k-point of input and in module are the same
      IF (ikks(iq) /= ikpt) &
        CALL errore(__FILE__, "unexpected mismatch of ikks and ikpt", 1)
      !
      ikq  = ikqs(iq)
      !
      CALL get_buffer(evq, lrwfc, iuwfc, ikq)
      !
      !!
      !! 3. construct the map from G and G' to G - G'
      !!
      CALL exchange_map(mill, igk_k(:, ikq), map)
      !!
      !! 4. construct the Coulomb potential
      !!
      ! q = k - (k - q)
      qvec = xk(:, ikpt) - xk(:, ikq)
      CALL exchange_coulomb(tpiba, truncation, sigma_x_st, qvec, g, coulomb, index_g_coul)
      !!
      !! 5. every process evaluates his contribution to sigma
      !!
      DO iband = iband_start, iband_stop
        !
        CALL exchange_convolution(wg(ikq, iband), coulomb, evq(:,iband), map, sigma)
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
      ! write unformatted
      CALL davcio(sigma, lrsex, iunsex, ikpt, 1)
      ! write formatted
      CALL sigma_io_write_x(output%unit_sigma, ikpt, sigma)
      !
    END IF ! meta_ionode

    CALL stop_clock(time_sigma_x)

  END SUBROUTINE exchange_wrapper

END MODULE exchange_module

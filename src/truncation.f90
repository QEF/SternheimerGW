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
!> Provides the routines to truncate a quantity to improve the k-point convergence.
MODULE truncation_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  !
  ! definition of different truncation methods
  !
  !> length of the truncation method
  INTEGER, PARAMETER :: trunc_length = 80

  !> no truncation - use bare Coulomb potential
  INTEGER, PARAMETER :: NO_TRUNCATION = 0
  !> \cond
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_1 = 'none'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_2 = 'off'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_3 = 'false'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_4 = 'no'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_5 = 'no truncation'
  !> \endcond

  !> spherical truncation - truncate potential at certain distance
  INTEGER, PARAMETER :: SPHERICAL_TRUNCATION = 1
  !> \cond
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_1 = 'on'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_2 = 'true'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_3 = 'yes'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_4 = 'spherical'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_5 = 'spherical truncation'
  !> \endcond

  !> film geometry truncation (expects film in x-y plane) -
  !! truncate potential at certain height
  INTEGER, PARAMETER :: FILM_TRUNCATION = 2
  !> \cond
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_1 = 'film'
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_2 = 'film truncation'
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_3 = '2d'
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_4 = '2d truncation'
  !> \endcond

  PRIVATE truncate_bare, truncate_spherical, truncate_film
 
CONTAINS

  !> Evaluate how the quantity associated with a reciprocal vector is truncated.
  !!
  !! Calculate a factor to scale a quantity defined in reciprocal space. There are
  !! different methods implemented to truncate the quantity, which are selected by
  !! the first parameter.
  !!
  !! \par No truncation 
  !! The bare Coulomb potential is used and only the divergent terms are removed.
  !! To activate this truncation scheme, set truncation to 'none', 'off', 'false',
  !! 'no', or 'no truncation' in the input file.
  !!
  !! \par Spherical truncation
  !! We truncate in real space at a certain radius R. In reciprocal space this leads
  !! to the prefactor \f$[1 - \cos(k R)]\f$. This prefactor expands to \f$(k R)^2 / 2\f$
  !! for small k cancelling the divergence of the Coulomb potential there. To activate
  !! this truncation scheme, set truncation to 'on', 'true', 'yes', 'spherical', or
  !! 'spherical truncation' in the input file.
  !!
  !! \par Film truncation
  !! We truncate at a certain height Z, which eliminates the divergence in reciprocal
  !! space. For details refer to Ismail-Beigi, Phys. Rev. B 73, 233103 (2006). To
  !! activate this truncation scheme, set truncation to 'film', '2d', 'film truncation',
  !! or '2d truncation' in the input file.
  !!
  !! \param[in] method Truncation method used; must be one of the integer constants
  !!            defined in this module. 
  !! \param[in] kpt Reciprocal lattice vector for which the quantity is truncated.
  !! \return Coulomb potential in reciprocal space scaled according to the specified
  !!         truncation scheme.
  REAL(dp) FUNCTION truncate(method, kpt) RESULT (factor)

    USE cell_base, ONLY: at, alat, omega
    USE constants, ONLY: fpi
    USE disp,      ONLY: nq1, nq2, nq3

    INTEGER,  INTENT(IN) :: method
    REAL(dp), INTENT(IN) :: kpt(3)

    REAL(dp) length_cut

    SELECT CASE (method)

    CASE (NO_TRUNCATION)
      factor = truncate_bare(kpt)

    CASE (SPHERICAL_TRUNCATION)

      ! cutoff radius
      length_cut = (3 * omega * nq1 * nq2 * nq3 / fpi)**(1.0 / 3.0)
      factor = truncate_spherical(kpt, length_cut)

    CASE (FILM_TRUNCATION)

      ! cutoff height
      length_cut = 0.5 * SQRT(SUM(at(:,3)**2)) * alat * nq3
      factor = truncate_film(kpt, length_cut)

    END SELECT ! method

  END FUNCTION truncate

  !> Implements the bare Coulomb potential.
  !! \param[in] kpt The point in reciprocal space.
  !! \see truncate for the details.
  REAL(dp) FUNCTION truncate_bare(kpt) RESULT (factor)

    USE constants, ONLY : eps8, fpi, e2

    REAL(dp), INTENT(IN) :: kpt(3)

    ! |k| and |k|^2
    REAL(dp) length_k, length_k2

    length_k2 = SUM(kpt**2)
    length_k = SQRT(length_k2)

    ! for large k vector
    IF (length_k > eps8) THEN

      ! bare Coulomb potential 4 pi e^2 / k^2
      factor = fpi * e2 / length_k2

    ! small k are removed
    ELSE

      factor = 0

    END IF

  END FUNCTION truncate_bare

  !> Implements the spherical truncation.
  !! \param[in] kpt The point in reciprocal space.
  !! \param[in] rcut The distance at which the potential is cut in real space.
  !! \see truncate
  REAL(dp) FUNCTION truncate_spherical(kpt, rcut) RESULT (factor)

    USE constants, ONLY: eps8, tpi, fpi, e2

    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: rcut

    ! |k| and |k|^2
    REAL(dp) length_k, length_k2

    length_k2 = SUM(kpt**2)
    length_k = SQRT(length_k2)

    ! for large k vector
    IF (length_k > eps8) THEN

      ! Coulomb potential 4 pi e^2 / k^2 is scaled by (1 - cos(k r))
      factor = fpi * e2 / length_k2 * (1 - COS(rcut * length_k))

    ! limit of small values
    ELSE

      ! (1 - cos(k r)) ~ (k r)^2 / 2
      ! with prefactor 4 pi e^2 / k^2, this yields 2 pi e^2 r^2
      factor = tpi * e2 * rcut**2

    END IF

  END FUNCTION truncate_spherical

  !> Implements the film truncation.
  !! \param[in] kpt The point in reciprocal space.
  !! \param[in] zcut The height at which the potential is cut in real space.
  !! \see truncate
  REAL(dp) FUNCTION truncate_film(kpt, zcut) RESULT (factor)

    USE constants, ONLY: eps6, tpi, fpi, e2

    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: zcut

    ! |k| and |k|^2
    REAL(dp) length_k, length_k2

    ! length of vector in z direction and in xy plane
    REAL(dp) length_kz, length_kxy

    length_k2 = SUM(kpt**2)
    length_k = SQRT(length_k2)
    length_kxy = SQRT(SUM(kpt(1:2)**2))
    length_kz = kpt(3)

    ! general case - large vector
    IF (length_k > eps6) THEN

      ! Coulomb potential 4 pi e^2 / k^2 is scaled by (1 - exp(-kxy z) * cos(kz z))
      factor = fpi * e2 / length_k2 * (1 - EXP(-length_kxy * zcut) * COS(length_kz * zcut))

    ! special case - small vector
    ELSE

      ! 1 - cos(kz z) ~ 1 - (kz z)^2 / 2 = (kz z)^2 / 2
      ! with prefactor 4 pi e^2 / k^2, this yields 2 pi e^2 z^2
      factor = tpi * e2 * zcut**2

    END IF

  END FUNCTION truncate_film

END MODULE truncation_module

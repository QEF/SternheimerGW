!> Provides the routines to truncate a quantity to improve the k-point convergence.
MODULE truncation_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  !
  ! definition of different truncation methods
  !
  !> length of the truncation method
  INTEGER, PARAMETER :: trunc_length = 80

  !> no truncation
  INTEGER, PARAMETER :: NO_TRUNCATION = 0
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_1 = 'none'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_2 = 'off'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_3 = 'false'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_4 = 'no'
  CHARACTER(LEN=trunc_length), PARAMETER :: NO_TRUNCATION_5 = 'no truncation'

  !> spherical truncation
  INTEGER, PARAMETER :: SPHERICAL_TRUNCATION = 1
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_1 = 'on'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_2 = 'true'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_3 = 'yes'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_4 = 'spherical'
  CHARACTER(LEN=trunc_length), PARAMETER :: SPHERICAL_TRUNCATION_5 = 'spherical truncation'

  !> film geometry truncation (expects film in x-y plane)
  INTEGER, PARAMETER :: FILM_TRUNCATION = 2
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_1 = 'film'
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_2 = 'film truncation'
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_3 = '2d'
  CHARACTER(LEN=trunc_length), PARAMETER :: FILM_TRUNCATION_4 = '2d truncation'

  PRIVATE truncate_spherical, truncate_film
 
CONTAINS

  !> Evaluate how the quantity associated with a reciprocal vector is truncated.
  !!
  !! Calculate a factor to scale a quantity defined in reciprocal space. There are
  !! different methods implemented to truncate the quantity, which are selected by
  !! the first parameter.
  !!
  !! \par[No truncation]
  !! The quantity remains unchanged, i.e., the function returns the value 1.
  !!
  !! \par[Spherical truncation]
  !!
  !! \par[Film truncation]
  !!
  !! \param[in] method Truncation method used; must be one of the integer constants defined in this module. 
  !! \param[in] kpt Reciprocal lattice vector for which the quantity is truncated.
  !! \return Factor with which the quantity should be scaled.
  REAL(dp) FUNCTION truncate(method, kpt) RESULT (factor)

    INTEGER,  INTENT(IN) :: method
    REAL(dp), INTENT(IN) :: kpt(3)

    SELECT CASE (method)

    CASE (NO_TRUNCATION)
      factor = 1.0

    CASE (SPHERICAL_TRUNCATION)
      factor = truncate_spherical()

    CASE (FILM_TRUNCATION)
      factor = truncate_film()

    END SELECT ! method

  END FUNCTION truncate

  !> Implements the spherical truncation.
  REAL(dp) FUNCTION truncate_spherical() RESULT (factor)
  END FUNCTION truncate_spherical

  !> Implements the film truncation.
  REAL(dp) FUNCTION truncate_film() RESULT (factor)
  END FUNCTION truncate_film

END MODULE truncation_module

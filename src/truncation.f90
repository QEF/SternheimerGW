!> Provides the routines to truncate a quantity to improve the k-point convergence.
MODULE truncation_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  !
  ! definition of different truncation methods
  !
  !> no truncation
  INTEGER, PARAMETER :: NO_TRUNCATION = 0
  !> spherical truncation
  INTEGER, PARAMETER :: SPHERICAL_TRUNCATION = 1
  !> film geometry truncation (expects film in x-y plane)
  INTEGER, PARAMETER :: FILM_TRUNCATION = 2

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

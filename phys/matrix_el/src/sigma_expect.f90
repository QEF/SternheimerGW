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
!> Evaluate expectation value of \f$\Sigma\f$.
!!
!! This module contains several routines to assist with the calculation of the
!! expectation value of a wavefunction \f$\phi\f$ with the self-energy
!! \f$\Sigma\f$. If the self-energy is given in imaginary frequency, we can use
!! a Pade approximation to perform the continuation to the real axis.
MODULE sigma_expect_mod

  USE kinds, ONLY: dp

  IMPLICIT NONE

  INTERFACE sigma_expect
    MODULE PROCEDURE sigma_expect_2d, sigma_expect_3d
  END INTERFACE sigma_expect

CONTAINS

  !> evaluate matrix elements of \f$\Sigma\f$ stored on file
  !!
  !! For a given set of wave functions, read the corresponding \f$\Sigma\f$ from the
  !! specified file and evaluate the matrix elements with it
  !! \param iunit unit pointing to the file in which \f$\Sigma\f$ is stored
  !! \param irec record which contains \f$\Sigma\f$
  !! \param wavef multiple wave functions for which the expectation value is computed
  !! \param ngm maximum G allowed
  !! \param igk map from G's to local k-point
  !! \param matel output the resulting matrix elements
  !! \param ndim third array dimension (optional: default = 1)
  SUBROUTINE sigma_expect_file(iunit,irec,wavef,ngm,igk,matel,ndim)

    USE reorder_mod, ONLY : reorder, create_map

    INTEGER,           INTENT(IN)  :: iunit
    INTEGER,           INTENT(IN)  :: irec
    COMPLEX(dp),       INTENT(IN)  :: wavef(:,:)
    INTEGER,           INTENT(IN)  :: ngm
    INTEGER,           INTENT(IN)  :: igk(:)
    COMPLEX(dp),       INTENT(OUT) :: matel(:,:,:)
    INTEGER, OPTIONAL, INTENT(IN)  :: ndim

    INTEGER ndim_loc
    INTEGER iloc
    INTEGER irec_loc

    INTEGER,     ALLOCATABLE :: map(:)
    COMPLEX(dp), ALLOCATABLE :: wavef_ordered(:,:), sigma(:,:,:)

    ! use ndim or default to 1
    IF ( PRESENT(ndim) ) THEN
      ndim_loc = ndim
    ELSE
      ndim_loc = 1
    END IF

    ! test array sizes
    CALL errore("sigma_expect_mod->sigma_expect_file", "array size mismatch", &
                size(matel) /= size(wavef,2)**2 * ndim_loc)

    !
    ! read sigma from file
    !

    ! allocate array to contain sigma
    ALLOCATE( sigma(ngm,ngm,ndim_loc) )

    ! read from file (factor 2 for complex)
    DO iloc = 1, ndim_loc
      irec_loc = (irec - 1) * ndim_loc + iloc
      CALL davcio( sigma(:,:,iloc), 2*ngm*ngm, iunit, irec_loc, -1 )
    END DO ! iloc

    !
    ! reorder wave function
    !

    ! create copy of wavef
    ALLOCATE( wavef_ordered(SIZE(wavef,1),SIZE(wavef,2)) )
    wavef_ordered = wavef

    ! create map to order wavef
    ALLOCATE( map(size(igk)) )
    map = create_map(igk,ngm)

    ! reorder wavef so that it is compatible with igk
    CALL reorder(wavef_ordered,map)
    
    !
    ! evaluate the expectation value
    !
    matel = sigma_expect( sigma, wavef_ordered(:ngm,:) )

    ! deallocate arrays
    DEALLOCATE( sigma )
    DEALLOCATE( wavef_ordered )
    DEALLOCATE( map )

  END SUBROUTINE sigma_expect_file

  !> Evaluate expectation value of \f$\Sigma\f$ for single wave function.
  !!
  !! \f{equation}{
  !!   \bigl\langle \phi_\text{l} \bigl\lvert \Sigma \bigr\rvert \phi_\text{r} \bigr\rangle
  !! \f}
  !! \param left_wavef wave function \f$\phi_\text{l}\f$
  !! \param sigma self-energy \f$\Sigma\f$
  !! \param right_wavef wave function \f$\phi_\text{r}\f$
  !! \return matrix element \f$\langle \phi_\text{l} \lvert \Sigma \rvert \phi_\text{r} \rangle\f$
  FUNCTION expectation(left_wavef,sigma,right_wavef) RESULT(energy)

    COMPLEX(dp), INTENT(IN)  :: sigma(:,:)
    COMPLEX(dp), INTENT(IN)  :: left_wavef(:)
    COMPLEX(dp), INTENT(IN)  :: right_wavef(:)
    COMPLEX(dp)              :: energy

    ! sanity check of the input
    CALL errore("sigma_expect_mod->expectation", "array size mismatch", &
                size(left_wavef) /= size(sigma,1))
    CALL errore("sigma_expect_mod->expectation", "array size mismatch", &
                size(right_wavef) /= size(sigma,2))

    ! evaluate < phi_l | Sigma | phi_r >
    energy = dot_product( left_wavef, matmul( sigma, right_wavef ) )

  END FUNCTION expectation

  !> Evaluate expectation value of \f$\Sigma\f$ for multiple wave functions.
  !!
  !! \f{equation}{
  !!   \bigl\langle \phi_n \bigl\lvert \Sigma \bigr\rvert \phi_m \bigr\rangle
  !! \f}
  !! \param sigma self-energy \f$\Sigma\f$
  !! \param wavef set of wave functions \f$\phi_n\f$
  !! \return matrix element \f$\langle \phi_n \lvert \Sigma \rvert \phi_m \rangle\f$
  FUNCTION sigma_expect_2d(sigma,wavef) RESULT (energy)

    COMPLEX(dp), INTENT(IN) :: sigma(:,:)
    COMPLEX(dp), INTENT(IN) :: wavef(:,:)

    COMPLEX(dp)             :: energy(size(wavef,2),size(wavef,2))

    INTEGER iband, jband

    ! loop over all bands
    DO jband = 1, size(wavef,2)
      DO iband = 1, size(wavef,2)
        energy(iband,jband) = expectation(wavef(:,jband),sigma,wavef(:,iband))
      END DO ! iband
    END DO ! jband

  END FUNCTION sigma_expect_2d

  !> Evaluate expectation value of \f$\Sigma\f$ at multiple frequencies and wave functions.
  !!
  !! \f{equation}{
  !!   \bigl\langle \phi_n \bigl\lvert \Sigma(\omega) \bigr\rvert \phi_m \bigr\rangle
  !! \f}
  !! \param sigma self-energy \f$\Sigma(\omega)\f$
  !! \param wavef set of wave functions \f$\phi_n\f$
  !! \return matrix element \f$\langle \phi_n \lvert \Sigma(\omega) \rvert \phi_m \rangle\f$
  FUNCTION sigma_expect_3d(sigma,wavef) RESULT (energy)

    COMPLEX(dp), INTENT(IN) :: sigma(:,:,:)
    COMPLEX(dp), INTENT(IN) :: wavef(:,:)

    COMPLEX(dp)             :: energy(size(wavef,2),size(wavef,2),size(sigma,3))

    INTEGER iband, jband, ifreq

    ! loop over all bands
    DO ifreq = 1, size(sigma,3)
      DO jband = 1, size(wavef,2)
        DO iband = 1, size(wavef,2)
          energy(iband,jband,ifreq) = expectation(wavef(:,jband),sigma(:,:,ifreq),wavef(:,iband))
        END DO ! iband
      END DO ! jband
    END DO ! nband

  END FUNCTION sigma_expect_3d

END MODULE sigma_expect_mod

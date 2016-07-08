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
!> Evaluate expectation value of \f$\Sigma\f$.
!!
!! This module contains several routines to assist with the calculation of the
!! expectation value of a wavefunction \f$\phi\f$ with the self-energy
!! \f$\Sigma\f$. If the self-energy is given in imaginary frequency, we can use
!! a Pade approximation to perform the continuation to the real axis.
MODULE sigma_expect_mod

  USE kinds, ONLY: dp

  IMPLICIT NONE

CONTAINS

  !> Evaluate expectation value of \f$\Sigma\f$.
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

END MODULE sigma_expect_mod

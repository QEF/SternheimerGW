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
!> Evaluate the self-energy \f$\Sigma = G W\f$
!!
!! There are two aspects to this product. One is the spacial component and the
!! other is the frequency one. In real space, no convolution is needed and
!! \f$\Sigma\f$ is just the product of \f$G\f$ and \f$W\f$. Hence, we take the
!! Fourier transform to get from reciprocal space where \f$G\f$ and \f$W\f$ are
!! defined to real space. Then we multiply both and Fourier transform \f$\Sigma\f$
!! back to reciprocal space.
!!
!! For the frequency an integration is necessary. There are two different
!! implementations to evaluate this integration to obtain the self-energy
!! \f$\Sigma(\omega)\f$. The essential difference is whether the frequencies
!! \f$\omega\f$ are on the real or the imaginary axis. Common to both cases is
!! that \f$W(\omega)\f$ is given on a few points on the imaginary axis. To
!! get \f$W\f$ for an arbitary frequency in the complex plane, we use the
!! Pade approximation.
!!
!! <h4>The real axis integration</h4>
!! If the frequencies of interest are on the real axis, we perform the following
!! integral [1]
!! \f{equation}{ \label{sig:real}
!!   \Sigma_k(\omega) = \frac{i}{N} \sum_{q} \int \frac{d\omega'}{2\pi} 
!!   G_{k-q}(\omega + \omega') W_{q}(\omega') e^{i \delta \omega'}~.
!! \f}
!! Here, \f$\delta\f$ is a small positive number that ensures the correct time
!! order. \f$N\f$ is the number of \f$q\f$ points in the calculation. Note, that
!! we suppressed the spacial coordiantes for clarity.
!!
!! <h4>The imaginary axis integration</h4>
!! The Fourier transform between frequency and time domain are defined as
!! \f{equation}{
!!   f(t) = \int \frac{d\omega}{2\pi} f(\omega) e^{i\omega t}
!! \f}
!! and
!! \f{equation}{
!!   f(\omega) = \int dt f(t) e^{-i\omega t}~.
!! \f}
!! One can extend this concept to time and frequency on the imaginary axis [2]
!! \f{equation}{
!!   f(it) = i \int \frac{d\omega}{2\pi} f(i\omega) e^{i\omega t}
!! \f}
!! and
!! \f{equation}{
!!   f(i\omega) = -i \int dt f(it) e^{-i\omega t}~.
!! \f}
!! Notice the occurence of the extra factors \f$\pm i\f$.
!!
!! Taking the Fourier transform of \eqref{sig:real}, we obtain due to the Fourier
!! convolution theorem
!! \f{equation}{
!!   \Sigma_k(t) = \frac{i}{N} \sum_q G_{k-q}(t) W_q(\delta - t)
!! \f}
!! the analytic continuation of this equation to the complex plane yields
!! \f{equation}{
!!   \Sigma_k(it) = \frac{i}{N} \sum_q G_{k-q}(it) W_q(\delta - it)~.
!! \f}
!! The equivalent of the Fourier convolution theorem for the complex Fourier
!! transform brings us back to frequency space
!! \f{equation}{
!!   \Sigma_k(i\omega) = -\frac{1}{N} \sum_q \int \frac{d \omega'}{2\pi}
!!   G_{k-q}(i\omega + i\omega') W_{q}(i\omega') e^{-\delta \omega'}.
!! \f}
!! 
!! <h4>References</h4>
!! [1] <a href="http://link.aps.org/doi/10.1103/PhysRevB.81.115105">
!!       Giustino, Cohen, Louie, Phys. Rev. B **81**, 115105 (2010)
!!     </a>
!!
!! [2] <a href="http://www.sciencedirect.com/science/article/pii/S001046559800174X">
!!       Rieger, *et al.*, Comput. Phys. Comm. **117**, 211 (1999)
!!     </a>
MODULE sigma_module

  IMPLICIT NONE

  PRIVATE

  PUBLIC sigma_real, sigma_imag

CONTAINS

  !> Implements the product on the real axis.
  SUBROUTINE sigma_real
  END SUBROUTINE sigma_real

  !> Implements the product on the imaginary axis.
  SUBROUTINE sigma_imag
  END SUBROUTINE sigma_imag

END MODULE sigma_module

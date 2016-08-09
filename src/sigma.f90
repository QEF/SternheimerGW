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

  PUBLIC sigma_real, sigma_imag, sigma_correlation

CONTAINS

  !> Implements the product on the real axis.
  SUBROUTINE sigma_real
  END SUBROUTINE sigma_real

  !> Implements the product on the imaginary axis.
  SUBROUTINE sigma_imag
  END SUBROUTINE sigma_imag

  !> Evaluate the product \f$\Sigma = \alpha G W\f$.
  !!
  !! The product is evaluated in real space. Because the Green's function can
  !! be used for several W values, we expect that the Green's function is
  !! already transformed to real space by the calling routine.
  SUBROUTINE sigma_prod(omega, fft_cust, alpha, gmapsym, green, array)

    USE fft_custom,     ONLY: fft_cus
    USE fft6_module,    ONLY: invfft6, fwfft6
    USE kinds,          ONLY: dp
    USE timing_module,  ONLY: time_GW_product

    !> volume of the unit cell
    REAL(dp),      INTENT(IN)    :: omega

    !> type that defines the custom Fourier transform for Sigma
    TYPE(fft_cus), INTENT(IN)    :: fft_cust

    !> The prefactor with which the self-energy is multiplied
    COMPLEX(dp),   INTENT(IN)    :: alpha

    !> the symmetry mapping from the reduced to the full G mesh
    INTEGER,       INTENT(IN)    :: gmapsym(:)

    !> the Green's function in real space \f$G(r, r')\f$
    COMPLEX(dp),   INTENT(IN)    :: green(:,:)

    !> *on input* the screened Coulomb interaction with both indices in reciprocal space <br>
    !! *on output* the self-energy with both indices in reciprocal space
    COMPLEX(dp),   INTENT(INOUT) :: array(:,:)

    !> the number of points in real space
    INTEGER num_r

    !> the number of G vectors defining the reciprocal space
    INTEGER num_g

    !> counter on G vectors
    INTEGER ig

    CALL start_clock(time_GW_product)

    ! determine the helper variables
    num_r = fft_cust%dfftt%nnr
    num_g = fft_cust%ngmt

    !
    ! sanity check of the input
    !
    IF (SIZE(green, 1) /= num_r) &
      CALL errore(__FILE__, "size of G inconsistent with FFT definition", 1)
    IF (SIZE(green, 2) /= num_r) &
      CALL errore(__FILE__, "Green's function not a square matrix", 1)
    IF (SIZE(array, 1) /= num_r) &
      CALL errore(__FILE__, "size of array inconsistent with FFT definition", 1)
    IF (SIZE(array, 2) /= num_r) &
      CALL errore(__FILE__, "array is not a square matrix", 1)
    IF (SIZE(gmapsym) /= num_g) &
      CALL errore(__FILE__, "gmapsym is inconsistent with FFT definition", 1)

    !!
    !! 1. We Fourier transform \f$W(G, G')\f$ to real space.
    !!
    ! array contains W(r, r')
    CALL invfft6('Custom', array, fft_cust%dfftt, fft_cust%nlt(gmapsym), omega)

    !!
    !! 2. We evaluate the product in real space 
    !!    \f$\Sigma(r, r') / \alpha = G(r, r') W(r, r') \f$
    !!
    ! array contains Sigma(r, r') / alpha
    array = green * array

    !!
    !! 3. The resulting is transformed back to reciprocal space \f$\Sigma(G, G')\f$.
    !!
    !. array contains Sigma(G, G') / alpha
    CALL fwfft6('Custom', array, fft_cust%dfftt, fft_cust%nlt, omega)

    !! 4. We multiply with the prefactor \f$\alpha\f$.
    DO ig = 1, num_g
      array(:num_g, ig) = alpha * array(:num_g, ig)
    END DO ! ig

    CALL stop_clock(time_GW_product)

  END SUBROUTINE sigma_prod

  !> Evaluate the correlation self-energy for a given frequency range.
  !!
  !! The algorithm consists of the following steps:
  SUBROUTINE sigma_correlation(omega, fft_cust, multishift, lmax, threshold,  &
                               mu, alpha, kpt, freq, first_sigma, &!eval, evec, &
                               gmapsym, coulomb, sigma)

    USE fft_custom,      ONLY: fft_cus
    USE fft6_module,     ONLY: invfft6
    USE freqbins_module, ONLY: freqbins_type
    USE green_module,    ONLY: green_prepare, green_function
    USE kinds,           ONLY: dp
    USE mp_images,       ONLY: inter_image_comm

    !> Volume of the unit cell
    REAL(dp),      INTENT(IN)  :: omega

    !> Type that defines the custom Fourier transform for Sigma
    TYPE(fft_cus), INTENT(IN)  :: fft_cust

    !> If this flag is set, we use the multishift solver
    LOGICAL,       INTENT(IN)  :: multishift

    !> The l value of the BiCGStab(l) algorithm.
    INTEGER,       INTENT(IN)  :: lmax

    !> The convergence threshold for the solver.
    REAL(dp),      INTENT(IN)  :: threshold

    !> The chemical potential of the system.
    REAL(dp),      INTENT(IN)  :: mu

    !> The prefactor with which the self-energy is multiplied
    COMPLEX(dp),   INTENT(IN)  :: alpha

    !> The k point at which the Green's function is evaluated
    REAL(dp),      INTENT(IN)  :: kpt(:)

    !> This type defines the frequencies used for the integration
    TYPE(freqbins_type), INTENT(IN) :: freq

    !> the first frequency done on this process
    INTEGER,       INTENT(IN)  :: first_sigma

!    !> The eigenvalues \f$\epsilon\f$ for the occupied bands.
!    COMPLEX(dp),   INTENT(IN)  :: eval(:)
!
!    !> The eigenvectors \f$u_{\text v}(G)\f$ of the occupied bands.
!    COMPLEX(dp),   INTENT(IN)  :: evec(:,:)

    !> the symmetry mapping from the reduced to the full G mesh
    INTEGER,       INTENT(IN)  :: gmapsym(:)

    !> The screened Coulomb interaction in reciprocal space
    COMPLEX(dp),   INTENT(IN)  :: coulomb(:,:,:)

    !> The self-energy \f$\Sigma\f$ at the specified frequency points
    COMPLEX(dp),   INTENT(OUT) :: sigma(:,:,:)

    !> the number of real space points used for the correlation
    INTEGER num_r_corr

    !> The number of plane-waves in the NSCF calculation.
    INTEGER num_g

    !> the number of G vectors used for the correlation
    INTEGER num_g_corr

    !> the number of tasks done by this process
    INTEGER num_task

    !> counter on the number of tasks
    INTEGER itask

    !> the number of frequencies for the Green's function
    INTEGER num_green

    !> counter on the Green's function frequencies
    INTEGER igreen

    !> counter on the Coulomb frequencies
    INTEGER icoul

    !> Counter for the frequencies of the self energy.
    INTEGER isigma

    !> The map from G-vectors at current k to global array.
    INTEGER,     ALLOCATABLE :: map(:)

    !> complex value of the chemical potential
    COMPLEX(dp) mu_

    !> The scaling prefactor in the integration, this is the product of alpha
    !! and the weight of the frequency point.
    COMPLEX(dp) alpha_weight

    !> The current frequency used for the Coulomb interaction
    COMPLEX(dp) freq_coul

    !> The frequencies used for the Green's function
    COMPLEX(dp), ALLOCATABLE :: freq_green(:)

    !> The frequencies used for the self energy.
    COMPLEX(dp), ALLOCATABLE :: freq_sigma(:)

    !> The Green's function in real or reciprocal space
    COMPLEX(dp), ALLOCATABLE :: green(:,:,:)

    !> work array; contains either W or \f$\Sigma\f$
    COMPLEX(dp), ALLOCATABLE :: work(:,:)

    !> complex constant of 0
    COMPLEX(dp), PARAMETER   :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    !
    ! sanity check of the input
    !
    num_r_corr = fft_cust%dfftt%nnr
    num_g_corr = fft_cust%ngmt
    num_task = SIZE(sigma, 3)
    IF (SIZE(coulomb, 1) /= num_g_corr) &
      CALL errore(__FILE__, "screened Coulomb and FFT type inconsistent", 1)
    IF (SIZE(coulomb, 2) /= num_g_corr) &
      CALL errore(__FILE__, "screened Coulomb must be a square matrix", 1)
    IF (SIZE(sigma, 1) /= num_g_corr) &
      CALL errore(__FILE__, "self energy and FFT type inconsistent", 1)
    IF (SIZE(sigma, 2) /= num_g_corr) &
      CALL errore(__FILE__, "self energy must be a square matrix", 1)
    IF (first_sigma + num_task - 1 > freq%num_sigma()) &
      CALL errore(__FILE__, "frequency index is out of bounds of the frequency array", 1)
    IF (SIZE(gmapsym) /= num_g_corr) &
      CALL errore(__FILE__, "gmapsym and FFT type are inconsistent", 1)

    !!
    !! 1. shift the frequency grid so that the origin is at the Fermi energy
    !!
    mu_ = CMPLX(mu, 0.0_dp)
    num_green = 2 * freq%num_coul()
    ALLOCATE(freq_green(num_green))
    ALLOCATE(freq_sigma(freq%num_sigma()))
    freq_green = freq%green(mu_)
    freq_sigma = mu_ + freq%sigma

    !!
    !! 2. prepare the QE module so that we can evaluate the Green's function
    !!
    ! this will allocate map
    CALL green_prepare(kpt, num_g_corr, map, num_g)

    !!
    !! 3. evaluate the Green's function of the system
    !!
    ! allocate an array that can contain the Fourier transformed quanity
    ALLOCATE(green(num_r_corr, num_r_corr, num_green))
    ! after this call, we obtained G(G, G', w)
    CALL green_function(inter_image_comm, multishift, lmax, threshold, map, &
                        num_g, freq_green, green)

    !!
    !! 4. we add the nonanalytic part if on the real axis
    !!
!    IF (.NOT. freq%imag_sigma) THEN
!      CALL green_nonanalytic(map, freq_green, eval, evec, green)
!    END IF

    !!
    !! 5. Fourier transform Green's function to real space
    !!
    ! the result is G(r, r', w)
    DO igreen = 1, num_green
      CALL invfft6('Custom', green(:,:,igreen), fft_cust%dfftt, fft_cust%nlt, omega)
    END DO ! igreen

    ! create work array
    ALLOCATE(work(num_r_corr, num_r_corr))

    !!
    !! 6. distribute the frequencies of \f$\Sigma\f$ across the image
    !!
    DO itask = 1, num_task
      !
      isigma = first_sigma + itask - 1
      !
      DO igreen = 1, num_green
        !!
        !! 7. construct W for the frequency \f$\omega^{\Sigma} - \omega^{\text{green}}\f$.
        !!
        work = zero
        freq_coul = freq_sigma(isigma) - freq_green(igreen)
        ! work will contain W(G, G', wS - wG)
        CALL construct_w(coulomb, work(1:num_g_corr, 1:num_g_corr), ABS(freq_coul))
        !!
        !! 8. convolute G and W
        !!
        icoul = MOD(igreen - 1, freq%num_coul()) + 1
        alpha_weight = alpha * freq%weight(icoul)
        ! work will contain Sigma(G, G', wS)
        CALL sigma_prod(omega, fft_cust, alpha_weight, gmapsym, green(:,:,igreen), work)
        !
        !!
        !! 9. add the result to \f$\Sigma\f$
        !!
        sigma(:,:,itask) = sigma(:,:,itask) + work(1:num_g_corr, 1:num_g_corr)
        !
      END DO ! igreen
      !
    END DO ! isigma

  END SUBROUTINE sigma_correlation

END MODULE sigma_module

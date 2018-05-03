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
!! integral [1,2]
!! \f{equation}{
!!   \Sigma_k(\omega) = \frac{i}{2 \pi N} \sum_{q} \int d\omega'
!!   G_{k-q}(\omega + \omega') W_{q}(\omega') e^{-i \delta \omega'}~.
!! \f}
!! Here, \f$\delta\f$ is a small positive number that ensures the correct time
!! order. \f$N\f$ is the number of \f$q\f$ points in the calculation. Note, that
!! we suppressed the spacial coordiantes for clarity. Because \f$W\f$ is evaluated
!! on the imaginary axis and continued to the real axis, it is convenient to
!! shift the integration boundaries
!! \f{equation}{ \label{sig:real}
!!   \Sigma_k(\omega) = \frac{i}{2 \pi N} \sum_{q} \int d\omega'
!!   G_{k-q}(\omega') W_{q}(\omega - \omega') e^{-i \delta (\omega - \omega')}~.
!! \f}
!! In this way, we can evaluate the Green's function on the integration grid once
!! and only the Coulomb potential is reevaluated for every frequency.
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
!! One can extend this concept to time and frequency on the imaginary axis [3]
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
!!   G_{k-q}(i\omega') W_{q}(i\omega - i\omega') e^{-\delta \omega'}.
!! \f}
!! Because the Green's function is evaluated on the imaginary axis (where it
!! is smooth), we can take the limit \f$\delta \rightarrow 0\f$ without
!! numerical problems.
!! 
!! <h4>References</h4>
!! [1] <a href="http://link.aps.org/doi/10.1103/PhysRevB.81.115105">
!!       Giustino, Cohen, Louie, Phys. Rev. B **81**, 115105 (2010)
!!     </a>
!!
!! [2] <a href="http://link.aps.org/doi/10.1103/PhysRevB.88.075117">
!!       Lambert, Giustino, Phys. Rev. B **88**, 075117 (2013)
!!     </a>
!!
!! [3] <a href="http://www.sciencedirect.com/science/article/pii/S001046559800174X">
!!       Rieger, *et al.*, Comput. Phys. Comm. **117**, 211 (1999)
!!     </a>
MODULE sigma_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC sigma_wrapper

CONTAINS

  !> This function provides a wrapper that extracts the necessary information
  !! from the global modules to evaluate the self energy.
  SUBROUTINE sigma_wrapper(ikpt, grid, config_green, freq, vcut, config, debug)

    USE analytic_module,      ONLY: analytic_coeff
    USE cell_base,            ONLY: omega
    USE constants,            ONLY: tpi
    USE control_gw,           ONLY: output, tmp_dir_coul, model_coul, tr2_gw
    USE coulpade_module,      ONLY: coulpade
    USE debug_module,         ONLY: debug_type, debug_set, test_nan
    USE disp,                 ONLY: x_q
    USE ener,                 ONLY: ef
    USE fft6_module,          ONLY: fft_map_generate
    USE freqbins_module,      ONLY: freqbins_type
    USE gvect,                ONLY: mill
    USE io_files,             ONLY: prefix
    USE io_global,            ONLY: meta_ionode, ionode_id, stdout
    USE kinds,                ONLY: dp
    USE klist,                ONLY: xk, lgauss
    USE mp,                   ONLY: mp_bcast
    USE mp_images,            ONLY: inter_image_comm, root_image
    USE mp_pools,             ONLY: inter_pool_comm, root_pool, me_pool
    USE output_mod,           ONLY: filcoul
    USE parallel_module,      ONLY: mp_root_sum
    USE select_solver_module, ONLY: select_solver_type
    USE setup_nscf_module,    ONLY: sigma_config_type
    USE sigma_grid_module,    ONLY: sigma_grid_type
    USE sigma_io_module,      ONLY: sigma_io_write_c
    USE symm_base,            ONLY: nsym, s, invs, ftau, nrot
    USE timing_module,        ONLY: time_sigma_c, time_sigma_setup, &
                                    time_sigma_io, time_sigma_comm
    USE truncation_module,    ONLY: vcut_type
    USE units_gw,             ONLY: iuncoul, lrcoul, iunsigma, lrsigma

    !> index of the k-point for which the self energy is evaluated
    INTEGER, INTENT(IN) :: ikpt

    !> the FFT grids used for the Fourier transformation
    TYPE(sigma_grid_type), INTENT(IN) :: grid

    !> the configuration of the linear solver for the Green's function
    TYPE(select_solver_type), INTENT(IN) :: config_green

    !> type containing the information about the frequencies used for the integration
    TYPE(freqbins_type), INTENT(IN) :: freq

    !> the truncated Coulomb potential
    TYPE(vcut_type), INTENT(IN) :: vcut

    !> evaluate the self-energy for these configurations
    TYPE(sigma_config_type), INTENT(IN) :: config(:)

    !> the debug configuration of the calculation
    TYPE(debug_type), INTENT(IN) :: debug

    !> complex constant of zero
    COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    !> real constant of 0.5
    REAL(dp),    PARAMETER :: half = 0.5_dp

    !> number of G and G' vectors in correlation grid
    INTEGER num_g_corr, num_gp_corr

    !> counter on the configurations
    INTEGER icon

    !> store the index of the q-point (to avoid rereading the Coulomb matrix)
    INTEGER iq

    !> complex prefactor <br>
    !! for real frequency integration it is \f$ i / 2 \pi N\f$ <br>
    !! for imaginary frequencies it is \f$-1 / 2 \pi N\f$
    COMPLEX(dp) prefactor

    !> the prefactor including the q dependent parts
    COMPLEX(dp) alpha

    !> the symmetry map from reduced to full G mesh
    INTEGER, ALLOCATABLE :: gmapsym(:,:)

    !> the phase associated with the symmetry
    COMPLEX(dp), ALLOCATABLE :: eigv(:,:)

    !> the highest and lowest occupied state (in Ry)
    REAL(dp) ehomo, elumo

    !> the chemical potential of the system
    REAL(dp) mu

    !> the map from local to global G vectors
    INTEGER, ALLOCATABLE :: fft_map(:)

    !> the screened Coulomb interaction
    COMPLEX(dp), ALLOCATABLE :: coulomb(:,:,:)

    !> the self-energy at the current k-point
    COMPLEX(dp), ALLOCATABLE :: sigma(:,:,:)

    !> the self-energy at the current k-point collected on the root process
    COMPLEX(dp), ALLOCATABLE :: sigma_root(:,:,:)

    !> number of bytes in a real
    INTEGER, PARAMETER   :: byte_real = 8

    !> name of file in which the Coulomb interaction is store
    CHARACTER(:), ALLOCATABLE :: filename

    !> record in which the result are written
    INTEGER irec

    !> counter on the frequencies
    INTEGER ifreq

    !> error flag from file opening
    INTEGER ierr

    !> flag indicating whether Coulomb file is already opened
    LOGICAL opend

    !> debug correlation part of self energy
    LOGICAL debug_sigma

    CALL start_clock(time_sigma_c)
    CALL start_clock(time_sigma_setup)

    WRITE(stdout, '(a)')
    WRITE(stdout, '(5x, a, 3f8.4, a)') 'evaluate self energy for k = (', xk(:, ikpt), ' )'

    ! set helper variable
    num_g_corr  = grid%corr_fft%ngm
    num_gp_corr = grid%corr_par_fft%ngm
    debug_sigma = debug_set .AND. debug%sigma_corr
    CALL fft_map_generate(grid%corr_par_fft, mill, fft_map)

    !
    ! set the prefactor depending on whether we integrate along the real or
    ! the imaginary frequency axis
    !
    IF (freq%imag_sigma) THEN
      prefactor = CMPLX(-1.0_dp / (tpi * REAL(nsym, KIND=dp)), 0.0_dp, KIND=dp)
    ELSE
      prefactor = CMPLX(0.0_dp, 1.0_dp / (tpi * REAL(nsym, KIND=dp)), KIND=dp)
    END IF

    !
    ! determine symmetry
    !
    ALLOCATE(gmapsym(num_g_corr, nrot))
    ALLOCATE(eigv(num_g_corr, nrot))
    CALL gmap_sym(num_g_corr, nsym, s, ftau, gmapsym, eigv, invs)
    DEALLOCATE(eigv)

    !
    ! define the chemical potential
    !
    IF (.NOT.lgauss) THEN
      !
      ! for semiconductors choose the middle of the gap
      CALL get_homo_lumo(ehomo, elumo)
      mu = half * (ehomo + elumo)
      !
    ELSE
      !
      ! for metals set it to the Fermi energy
      mu = ef
      !
    END IF
    CALL mp_bcast(mu, ionode_id, inter_pool_comm)

    CALL stop_clock(time_sigma_setup)

    !
    ! initialize self energy
    !
    ALLOCATE(sigma(num_g_corr, num_gp_corr, freq%num_sigma()))
    sigma = zero

    !
    ! open the file containing the coulomb interaction
    ! TODO a wrapper routine for all I/O
    !
    INQUIRE(UNIT = iuncoul, OPENED = opend)
    IF (.NOT. opend) THEN
      !
      filename = TRIM(tmp_dir_coul) // TRIM(prefix) // "." // TRIM(filcoul) // "1"
      OPEN(UNIT = iuncoul, FILE = filename, IOSTAT = ierr, &
           ACCESS = 'direct', STATUS = 'old', RECL = byte_real * lrcoul)
      CALL errore(__FILE__, "error opening " // filename, ierr)
      !
    END IF ! open file
    !
    ! initialize index of q (set to 0 so that it is always read on first call)
    iq = 0
    !
    ! allocate array for the coulomb matrix
    ALLOCATE(coulomb(num_g_corr, num_g_corr, freq%num_freq()))

    !
    ! sum over all q-points
    !
    DO icon = 1, SIZE(config)
      !
      WRITE(stdout, '(5x,a,3f8.4,a)', ADVANCE='NO') 'k + q = (', xk(:,config(icon)%index_kq), ' )'
      !
      ! evaluate the coefficients for the analytic continuation of W
      !
      IF (config(icon)%index_q /= iq) THEN
        !
        iq = config(icon)%index_q
        CALL davcio(coulomb, lrcoul, iuncoul, iq, -1)
        CALL coulpade(x_q(:,iq), vcut, coulomb)
        CALL analytic_coeff(model_coul, tr2_gw, freq, coulomb)
        !
        ! check if any NaN occured in coulpade
        IF (debug_sigma) THEN
          IF (ANY(test_nan(coulomb))) THEN
            CALL errore(__FILE__, 'Found a NaN in Coulomb after analytic continuation', iq)
          END IF
        END IF
        !
      END IF
      !
      ! determine the prefactor
      !
      alpha = config(icon)%weight * prefactor
      !
      ! evaluate Sigma
      !
      CALL sigma_correlation(omega, grid, config_green,                 &
                             mu, alpha, config(icon)%index_kq, freq,    &
                             gmapsym(:, config(icon)%sym_op), &
                             coulomb, sigma, debug)
      !
    END DO ! icon

    DEALLOCATE(gmapsym)
    DEALLOCATE(coulomb)

    !
    ! collect sigma on a single process
    !
    ! TODO replace this by parallel I/O
    CALL start_clock(time_sigma_comm)

    ! first sum sigma across the pool
    CALL mp_root_sum(inter_pool_comm, root_pool, sigma)

    ! unpack sigma in the large array
    IF (me_pool == root_pool) THEN
      !
      ALLOCATE(sigma_root(num_g_corr, num_g_corr, freq%num_sigma()), STAT = ierr)
      IF (ierr /= 0) THEN
        CALL errore(__FILE__, "error allocating array to collect sigma", ierr)
        RETURN
      END IF
      !
      sigma_root = zero
      !
      DO ifreq = 1, freq%num_sigma()
        sigma_root(:, fft_map, ifreq) = sigma(:,:,ifreq)
      END DO
      !
      CALL mp_root_sum(inter_image_comm, root_image, sigma_root)
      !
    END IF

    DEALLOCATE(sigma)

    CALL stop_clock(time_sigma_comm)

    !
    ! the root process writes sigma to file
    !
    CALL start_clock(time_sigma_io)
    IF (meta_ionode .AND. ALLOCATED(sigma_root)) THEN
      !
      DO ifreq = 1, freq%num_sigma()
        irec = (ikpt - 1) * freq%num_sigma() + ifreq
        CALL davcio(sigma_root(:,:,ifreq), lrsigma, iunsigma, irec, 1)
      END DO ! ifreq
      !
      CALL sigma_io_write_c(output%unit_sigma, ikpt, sigma_root)
      !
    END IF ! ionode
    CALL stop_clock(time_sigma_io)

    CALL stop_clock(time_sigma_c)

  END SUBROUTINE sigma_wrapper

  !> Evaluate the product \f$\Sigma = \alpha G W\f$.
  !!
  !! The product is evaluated in real space. Because the Green's function can
  !! be used for several W values, we expect that the Green's function is
  !! already transformed to real space by the calling routine.
  SUBROUTINE sigma_prod(omega, grid, alpha, green, array, debug)

    USE debug_module,      ONLY: debug_type, debug_set, test_nan
    USE fft6_module,       ONLY: invfft6, fwfft6
    USE kinds,             ONLY: dp
    USE sigma_grid_module, ONLY: sigma_grid_type
    USE timing_module,     ONLY: time_GW_product

    !> volume of the unit cell
    REAL(dp),      INTENT(IN)    :: omega

    !> type that defines the custom Fourier transform for Sigma
    TYPE(sigma_grid_type), INTENT(IN) :: grid

    !> The prefactor with which the self-energy is multiplied
    COMPLEX(dp),   INTENT(IN)    :: alpha

    !> the Green's function in real space \f$G(r, r')\f$
    COMPLEX(dp),   INTENT(IN)    :: green(:,:)

    !> *on input* the screened Coulomb interaction with both indices in reciprocal space <br>
    !! *on output* the self-energy with both indices in reciprocal space
    COMPLEX(dp),   INTENT(INOUT) :: array(:,:)

    !> the debug configuration of the calculation
    TYPE(debug_type), INTENT(IN) :: debug

    !> the number of points in real space
    INTEGER num_r, num_rp

    !> the number of G vectors defining the reciprocal space
    INTEGER num_g, num_gp

    !> counter on G vectors
    INTEGER ig

    !> debug the convolution of G and W
    LOGICAL debug_sigma

    CALL start_clock(time_GW_product)

    ! determine the helper variables
    num_r  = grid%corr_fft%nnr
    num_rp = grid%corr_par_fft%nnr
    num_g  = grid%corr_fft%ngm
    num_gp = grid%corr_par_fft%ngm
    debug_sigma = debug_set .AND. debug%sigma_corr

    !
    ! sanity check of the input
    !
    IF (SIZE(green, 1) /= num_r) &
      CALL errore(__FILE__, "size of G inconsistent with G-vector FFT definition", 1)
    IF (SIZE(green, 2) /= num_rp) &
      CALL errore(__FILE__, "size of G inconsistent with G'-vector FFT definition", 1)
    IF (SIZE(array, 1) /= num_r) &
      CALL errore(__FILE__, "size of array inconsistent with G-vector FFT definition", 1)
    IF (SIZE(array, 2) /= num_rp) &
      CALL errore(__FILE__, "size of array inconsistent with G'-vector FFT definition", 1)
    IF (test_nan(alpha)) &
      CALL errore(__FILE__, "prefactor of the convolution is NaN", 1)

    !!
    !! 1. We Fourier transform \f$W(G, G')\f$ to real space.
    !!
    ! array contains W(r, r')
    CALL invfft6('Rho', array, grid%corr_fft, grid%corr_par_fft, omega)
    ! check for NaN after FFT
    IF (debug_sigma) THEN
      IF (ANY(test_nan(array))) THEN
        CALL errore(__FILE__, "FFT of W resulted in NaN", 1)
      END IF
    END IF

    !!
    !! 2. We evaluate the product in real space 
    !!    \f$\Sigma(r, r') / \alpha = G(r, r') W(r, r') \f$
    !!
    ! array contains Sigma(r, r') / alpha
    array = green * array
    ! check for NaN after multiplication
    IF (debug_sigma) THEN
      IF (ANY(test_nan(array))) THEN
        CALL errore(__FILE__, "G * W produced a NaN", 1)
      END IF
    END IF

    !!
    !! 3. The resulting is transformed back to reciprocal space \f$\Sigma(G, G')\f$.
    !!
    !. array contains Sigma(G, G') / alpha
    CALL fwfft6('Rho', array, grid%corr_fft, grid%corr_par_fft, omega)
    ! check for NaN after the FFT
    IF (debug_sigma) THEN
      IF (ANY(test_nan(array(:num_g, :num_gp)))) THEN
        CALL errore(__FILE__, "FFT of Sigma resulted in NaN", 1)
      END IF
    END IF

    !! 4. We multiply with the prefactor \f$\alpha\f$.
    DO ig = 1, num_gp
      array(:num_g, ig) = alpha * array(:num_g, ig)
    END DO ! ig

    CALL stop_clock(time_GW_product)

  END SUBROUTINE sigma_prod

  !> Evaluate the correlation self-energy for a given frequency range.
  !!
  !! The algorithm consists of the following steps:
  SUBROUTINE sigma_correlation(omega, grid, config, mu, alpha, ikq, freq, &
                               gmapsym, coulomb, sigma, debug)

    USE analytic_module,      ONLY: analytic_eval
    USE debug_module,         ONLY: debug_type, debug_set, test_nan
    USE fft6_module,          ONLY: invfft6
    USE freqbins_module,      ONLY: freqbins_type
    USE green_module,         ONLY: green_prepare, green_function
    USE io_global,            ONLY: stdout
    USE kinds,                ONLY: dp
    USE select_solver_module, ONLY: select_solver_type
    USE sigma_grid_module,    ONLY: sigma_grid_type
    USE timing_module,        ONLY: time_sigma_c

    !> Volume of the unit cell
    REAL(dp),      INTENT(IN)  :: omega

    !> Type that defines the custom Fourier transform for Sigma
    TYPE(sigma_grid_type), INTENT(IN) :: grid

    !> The configuration of the linear solver
    TYPE(select_solver_type), INTENT(IN) :: config

    !> The chemical potential of the system.
    REAL(dp),      INTENT(IN)  :: mu

    !> The prefactor with which the self-energy is multiplied
    COMPLEX(dp),   INTENT(IN)  :: alpha

    !> The index of the point k - q at which the Green's function is evaluated
    INTEGER,       INTENT(IN)  :: ikq

    !> This type defines the frequencies used for the integration
    TYPE(freqbins_type), INTENT(IN) :: freq

    !> the symmetry mapping from the reduced to the full G mesh
    INTEGER,       INTENT(IN)  :: gmapsym(:)

    !> The screened Coulomb interaction in reciprocal space
    COMPLEX(dp),   INTENT(IN)  :: coulomb(:,:,:)

    !> The self-energy \f$\Sigma\f$ at the specified frequency points
    COMPLEX(dp),   INTENT(INOUT) :: sigma(:,:,:)

    !> the debug configuration of the calculation
    TYPE(debug_type), INTENT(IN) :: debug

    !> debug correlation part of self energy
    LOGICAL debug_sigma

    !> the number of real space points used for the correlation
    INTEGER num_r_corr, num_rp_corr

    !> The number of plane-waves in the NSCF calculation.
    INTEGER num_g

    !> the number of G vectors used for the correlation
    INTEGER num_g_corr, num_gp_corr

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

    !> measure the time since start of the routine
    REAL(dp) get_clock

    !> time when routine started
    REAL(dp) start_time

    start_time = get_clock(time_sigma_c)

    !
    ! sanity check of the input
    !
    debug_sigma = debug_set .AND. debug%sigma_corr
    num_r_corr  = grid%corr_fft%nnr
    num_rp_corr = grid%corr_par_fft%nnr
    num_g_corr  = grid%corr_fft%ngm
    num_gp_corr = grid%corr_par_fft%ngm
    IF (SIZE(coulomb, 1) /= num_g_corr) &
      CALL errore(__FILE__, "screened Coulomb and G-vector FFT type inconsistent", 1)
    IF (SIZE(coulomb, 2) /= num_g_corr) &
      CALL errore(__FILE__, "screened Coulomb not a square matrix", 1)
    IF (SIZE(sigma, 1) /= num_g_corr) &
      CALL errore(__FILE__, "self energy and G-vector FFT type inconsistent", 1)
    IF (SIZE(sigma, 2) /= num_gp_corr) &
      CALL errore(__FILE__, "self energy and G'-vector FFT type inconsistent", 1)
    IF (SIZE(sigma, 3) /= freq%num_sigma()) &
      CALL errore(__FILE__, "frequency dimension of self energy not correct size", 1)
    IF (SIZE(gmapsym) /= num_g_corr) &
      CALL errore(__FILE__, "gmapsym and FFT type are inconsistent", 1)

    !!
    !! 1. shift the frequency grid so that the origin is at the Fermi energy
    !!
    mu_ = CMPLX(mu, 0.0_dp, KIND=dp)
    num_green = 2 * freq%num_coul()
    ALLOCATE(freq_green(num_green))
    ALLOCATE(freq_sigma(freq%num_sigma()))
    freq_green = freq%green(mu_)
    freq_sigma = mu_ + freq%sigma

    !!
    !! 2. prepare the QE module so that we can evaluate the Green's function
    !!
    ! this will allocate map
    CALL green_prepare(ikq, num_g_corr, map, num_g)

    !!
    !! 3. evaluate the Green's function of the system
    !!
    ! after this call, we obtained G(G, G', w)
    CALL green_function(grid, config, map, num_g, freq_green, green, debug)

    ! check for NaN in Green's function
    IF (debug_sigma) THEN
      IF (ANY(test_nan(green))) THEN
        CALL errore(__FILE__, "Green's function contains NaN in reciprocal space", ikq)
      END IF
    END IF


    WRITE(stdout,'(2x,a,f9.2,a)', ADVANCE='NO') 'G: ', get_clock(time_sigma_c) - start_time, 's'
    start_time = get_clock(time_sigma_c)

    !!
    !! 4. Fourier transform Green's function to real space
    !!
    ! the result is G(r, r', w)
    DO igreen = 1, num_green
      CALL invfft6('Rho', green(:,:,igreen), grid%corr_fft, grid%corr_par_fft, omega)
    END DO ! igreen


    ! check for NaN in Green's function
    IF (debug_sigma) THEN
      IF (ANY(test_nan(green))) THEN
        CALL errore(__FILE__, "Green's function contains NaN in real space", ikq)
      END IF
    END IF

    !!
    !! 5. distribute the frequencies of \f$\Sigma\f$ across the image
    !!
    DO isigma = 1, freq%num_sigma()
      !
      DO igreen = 1, num_green
        !!
        !! 6. construct W for the frequency \f$\omega^{\Sigma} - \omega^{\text{green}}\f$.
        !!
        freq_coul = freq_sigma(isigma) - freq_green(igreen)
        ! ensure frequency is in the correct quadrant
        IF (REAL(freq_coul) * AIMAG(freq_coul) < 0.0_dp) freq_coul = CONJG(freq_coul)
        !
        ! work will be allocated and contain W(G, G', wS - wG)
        CALL analytic_eval(gmapsym, grid, freq, coulomb, freq_coul, work)
        !
        ! check for NaN in screened Coulomb
        IF (debug_sigma) THEN
          IF (ANY(test_nan(work))) THEN
            CALL errore(__FILE__, "screened Coulomb interaction contains NaN", igreen)
          END IF
        END IF
        !!
        !! 7. convolute G and W
        !!
        icoul = MOD(igreen - 1, freq%num_coul()) + 1
        alpha_weight = alpha * freq%weight(icoul)
        ! work will contain Sigma(G, G', wS)
        CALL sigma_prod(omega, grid, alpha_weight, green(:,:,igreen), work, debug)
        !
        ! check for NaN after convolution
        IF (debug_sigma) THEN
          IF (ANY(test_nan(work(:num_g_corr,:num_gp_corr)))) THEN
            CALL errore(__FILE__, "convolution of G and W introduced NaN", igreen)
          END IF
        END IF
        !!
        !! 8. add the result to \f$\Sigma\f$
        !!
        sigma(:,:,isigma) = sigma(:,:,isigma) + work(1:num_g_corr, 1:num_gp_corr)
        !
      END DO ! igreen
      !
    END DO ! isigma

    WRITE(stdout,'(2x,a,f9.2,a)') 'G*W:', get_clock(time_sigma_c) - start_time, 's'
    FLUSH(stdout)

  END SUBROUTINE sigma_correlation

END MODULE sigma_module

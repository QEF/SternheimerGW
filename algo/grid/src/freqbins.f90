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
!> Contains routine and type to define frequencies
MODULE freqbins_module

  USE kinds, ONLY: dp

  IMPLICIT NONE

  !> do not use symmetry
  INTEGER, PARAMETER :: no_symmetry = 0

  !> double the number of frequency points by using \f$W(\omega) = W(-\omega)\f$
  INTEGER, PARAMETER :: even_symmetry = 1

  !> fit a function of \f$\omega^2\f$ so that W is even by construction
  INTEGER, PARAMETER :: square_symmetry = 2

  PRIVATE

  PUBLIC freqbins_type, freqbins, freqbins_symm, no_symmetry

  !> Contains the information about the frequencies used for the convolution of
  !! G and W to obtain the self-energy \f$\Sigma\f$.
  TYPE freqbins_type

    !> flag indicates whether we integrate on the real or imaginary axis
    LOGICAL imag_sigma

    !> indicates what symmetrization is used for frequencies of screened Coulomb
    INTEGER freq_symm_coul

    !> number of frequencies in the symmetrized mesh
    INTEGER num_freq_sym

    !> The coarse frequency grid used to determine W.
    COMPLEX(dp), ALLOCATABLE :: solver(:)

    !> The dense frequency grid onto which W is analytically continued.
    COMPLEX(dp), ALLOCATABLE :: coul(:)

    !> The dense frequency grid used for the self-energy.
    COMPLEX(dp), ALLOCATABLE :: sigma(:)

    !> The weights for the frequency integration.
    REAL(dp),    ALLOCATABLE :: weight(:)

    !> The grid used for the output of the self-energy.
    REAL(dp),    ALLOCATABLE :: window(:)

    !> The small shift from the real axis into the complex plane
    REAL(dp) eta

  CONTAINS

    PROCEDURE :: init       => freqbins_init
    PROCEDURE :: num_freq   => freqbins_num_freq
    PROCEDURE :: num_coul   => freqbins_num_coul
    PROCEDURE :: num_sigma  => freqbins_num_sigma
    PROCEDURE :: num_window => freqbins_num_window
    PROCEDURE :: green      => freqbins_green
    PROCEDURE :: symmetrize => freqbins_symmetrize

  END TYPE freqbins_type

CONTAINS

  !> Generate frequency bins
  !!
  !! We assume that the self-energy \f$\Sigma\f$ is needed for frequencies
  !! \f$\omega^\Sigma_\text{min} \le \omega^\Sigma \le \omega^\Sigma_\text{max}\f$.
  !! The convolution requires W for positive frequencies[1]
  !! \f$0 \le \omega^\text{coul} \le \omega^\text{coul}_\text{max}\f$.
  !! The Green's function is needed for
  !! \f{equation}{
  !!   \omega^\Sigma_\text{min} - \omega^\text{coul}_\text{max} \le
  !!   \omega^\text{green} \le \omega^\Sigma_\text{max} - \omega^\text{coul}_\text{max}~.
  !! \f}
  !! Because of the multishift algorithm, the frequency dependence of the
  !! Green's function is inexpensive, so we use the same spacing.
  !!
  !! @note We require \f$\omega^\text{coul}_\text{max} > 0 \f$, 
  !! \f$\omega^\Sigma_\text{min} \le 0\f$, and \f$\omega^\Sigma_\text{max} > 0\f$.
  !! The zero of energy is set to the Fermi level.
  !!
  !! All input frequencies should be set in Ry and only the nontrivial real or
  !! imaginary part should be passed. The code will automatically add the trivial
  !! zero part depending on the imag_sigma flag.
  !!
  !! The algorithm performs the following steps:
  !!
  ! TODO should set two frequency windows one fine grid for the range around the fermi level
  !      say  ef +/- 60 eV  down to the lowest pseudo state included! and a second 
  !      course window for everything outside this range.
  SUBROUTINE freqbins(imag_sigma, min_sigma, max_sigma, num_sigma, max_coul, num_coul, &
                      min_window, max_window, num_window, freq)
 
    IMPLICIT NONE 

    !> evaluate the self-energy on the real or imaginary axis
    LOGICAL,  INTENT(IN) :: imag_sigma

    !> minimal frequency for the self-energy (in Ry)
    REAL(dp), INTENT(IN) :: min_sigma

    !> maximal frequency for the self-energy (in Ry)
    REAL(dp), INTENT(IN) :: max_sigma

    !> number of frequency points for self-energy
    INTEGER,  INTENT(IN) :: num_sigma

    !> maximal frequency for the screened Coulomb potential (in Ry)
    REAL(dp), INTENT(IN) :: max_coul

    !> number of frequency points for the screened Coulomb potential
    INTEGER,  INTENT(IN) :: num_coul

    !> minimal frequency for the output window (in Ry)
    REAL(dp), INTENT(IN) :: min_window

    !> maximal frequency for the output window (in Ry)
    REAL(dp), INTENT(IN) :: max_window

    !> number of frequency points in output window
    INTEGER,  INTENT(IN) :: num_window

    !> contains the frequency information
    TYPE(freqbins_type), INTENT(INOUT) :: freq

    !> real constant of zero
    REAL(dp), PARAMETER :: zero = 0.0_dp

    !> temporary array for the frequency grid
    REAL(dp), ALLOCATABLE :: grid(:)

    !!
    !! 1. We create equidistant frequencies for \f$\Sigma\f$
    !!
    CALL freqbins_equidist_grid(min_sigma, max_sigma, num_sigma, grid)
    !!
    !! 2. Depending on the imag_sigma flag, we initialize either the
    !!    real or imaginary part.
    !!
    ALLOCATE(freq%sigma(num_sigma))
    !
    IF (imag_sigma) THEN 
      freq%sigma = CMPLX(zero, grid, KIND = dp)
    ELSE
      freq%sigma = CMPLX(grid, zero, KIND = dp)
    END IF
    !
    DEALLOCATE(grid)
    !!
    !! 3. We construct a frequency grid for W.
    !!
    freq%imag_sigma = imag_sigma
    IF (.NOT.imag_sigma) THEN
      !!
      !! - if imag_sigma is cleared, we use an equidistant grid on the real axis
      !!
      CALL freqbins_equidist_grid(zero, max_coul, num_coul, grid)
      !
      ALLOCATE(freq%coul(num_coul))
      ALLOCATE(freq%weight(num_coul))
      !
      freq%coul = CMPLX(grid, freq%eta, KIND = dp)
      ! first point has half the weight because we will use it twice
      freq%weight = 2.0 * max_coul / REAL(2 * num_coul - 1, KIND = dp)
      freq%weight(1) = 0.5 * freq%weight(1)
      !
      DEALLOCATE(grid)
      !
    ELSE
      !!
      !! - otherwise we use a Gauss-Legendre grid on the imaginary axis
      !!
      ALLOCATE(grid(num_coul))
      ALLOCATE(freq%weight(num_coul))
      !
      CALL gauleg_grid(zero, max_coul, grid, freq%weight, num_coul)
      !
      ALLOCATE(freq%coul(num_coul))
      !
      freq%coul = CMPLX(zero, grid, KIND = dp)
      !
    END IF
    !!
    !! 4. create grid for output
    !!
    CALL freqbins_equidist_grid(min_window, max_window, num_window, freq%window)

    ! initialize the type
    CALL freq%init()

    !! References
    !! [1] <a href="http://link.aps.org/doi/10.1103/PhysRevB.74.035101">
    !!     Shishkin, Kresse, Phys. Rev. B, **74**, 035101 (2006)
    !!     </a>

  END SUBROUTINE freqbins

  !> create an equidistant grid
  SUBROUTINE freqbins_equidist_grid(min_freq, max_freq, num_freq, grid)

    !> the minimal value of the grid
    REAL(dp), INTENT(IN) :: min_freq

    !> the maximal value of the grid
    REAL(dp), INTENT(IN) :: max_freq

    !> the number of grid points
    INTEGER,  INTENT(IN) :: num_freq

    !> the generated grid
    REAL(dp), INTENT(OUT), ALLOCATABLE :: grid(:)

    !> the spacing between the grid points
    REAL(dp) delta_freq

    !> counter for the array
    INTEGER ifreq

    ! determine the spacing of the grid points
    ! note: subtract 1, so that the boundaries are included
    delta_freq = (max_freq - min_freq) / (num_freq - 1)

    ! allocate the grid
    ALLOCATE(grid(num_freq))

    ! generate the grid
    DO ifreq = 1, num_freq
      grid(ifreq) = min_freq + delta_freq * (ifreq - 1)
    END DO ! igrid

  END SUBROUTINE freqbins_equidist_grid

  !> extract the number of frequencies in coarse mesh for W
  INTEGER FUNCTION freqbins_num_freq(this) RESULT (num_freq)

    !> The frequency type of which the number of frequencies are extracted
    CLASS(freqbins_type), INTENT(IN) :: this

    ! the init routine of the type needs to be called before this routine
    IF (this%num_freq_sym < 0) &
       CALL errore(__FILE__, "freqbins_type used before it's init routine was called", 1)

    num_freq = this%num_freq_sym

  END FUNCTION freqbins_num_freq

  !> extract the number of frequencies in the dense mesh of W
  INTEGER FUNCTION freqbins_num_coul(this) RESULT (num_coul)

    !> The frequency type of which the number of frequencies are extracted
    CLASS(freqbins_type), INTENT(IN) :: this

    num_coul = SIZE(this%coul)

  END FUNCTION freqbins_num_coul

  !> extract the number of frequencies in the dense mesh of Sigma
  INTEGER FUNCTION freqbins_num_sigma(this) RESULT (num_sigma)

    !> The frequency type of which the number of frequencies are extracted
    CLASS(freqbins_type), INTENT(IN) :: this

    num_sigma = SIZE(this%sigma)

  END FUNCTION freqbins_num_sigma

  !> extract the number of frequencies in the output window
  INTEGER FUNCTION freqbins_num_window(this) RESULT (num_window)

    !> The frequency type of which the number of frequencies are extracted
    CLASS(freqbins_type), INTENT(IN) :: this

    num_window = SIZE(this%window)

  END FUNCTION freqbins_num_window

  !> extract the frequencies used for the Green's function
  FUNCTION freqbins_green(this, freq_sigma) RESULT (freq_green)

    !> The frequency type for which, we want to evaluate the Green's function frequencies
    CLASS(freqbins_type), INTENT(IN) :: this

    !> frequency of self-energy relative to which we determine the frequencies
    COMPLEX(dp), INTENT(IN)  :: freq_sigma

    !> the frequencies used for the Green's function
    COMPLEX(dp), ALLOCATABLE :: freq_green(:)

    ! helper for number of frequencies for W
    INTEGER num_coul

    num_coul = this%num_coul()

    ALLOCATE(freq_green(2 * num_coul))
    freq_green(:num_coul) = freq_sigma + this%coul
    freq_green(num_coul + 1:) = freq_sigma - this%coul

  END FUNCTION freqbins_green

  !> Expand an array from the symmetry reduced to the full frequency mesh.
  !!
  !! We assume an array with the following property
  !! \f{equation}{
  !!   A(G, G', \omega) = A(G, G', -\omega)~.
  !! \f}
  !! Then, we can extend an given input array for which only the first half
  !! of the frequencies was calculated to the full mesh by adding the complex
  !! conjugate to the end.
  SUBROUTINE freqbins_symm(freq_in, freq_out, array)

    USE constants, ONLY: eps14
    USE kinds,     ONLY: dp

    !> the definition of the input frequency mesh
    TYPE(freqbins_type), INTENT(IN)       :: freq_in

    !> the full frequency array
    COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: freq_out(:)

    !> *on input*: The first half of the array \f$A\f$ <br>
    !! *on output*: The full array \f$A\f$ (extended by its complex conjugate).
    COMPLEX(dp), INTENT(INOUT), OPTIONAL  :: array(:,:,:)

    !> original number of frequencies
    INTEGER num_freq

    !> number of frequencies that are vanishingly small
    INTEGER num_zero

    !> new number of frequencies
    INTEGER num_freq_sym

    !> counter on initial frequency mesh
    INTEGER ifreq

    !> counter on final frequency mesh
    INTEGER ifreq_sym

    ! set helper variable
    num_freq = SIZE(freq_in%solver)

    ! trivial case - symmetry not used => return same frequency used for solver
    IF (freq_in%freq_symm_coul == no_symmetry) THEN
      ALLOCATE(freq_out(num_freq))
      freq_out = freq_in%solver
      RETURN
    END IF

    ! fit to w^2 so that even symmetry is fulfilled by construction
    IF (freq_in%freq_symm_coul == square_symmetry) THEN
      ALLOCATE(freq_out(num_freq))
      freq_out = freq_in%solver**2
      RETURN
    END IF

    ! else -> even_symmetry
    !!
    !! For the symmetrized frequency mesh, we duplicate all nonzero frequencies.
    !! We only for at most one frequency that is 0.
    !!
    num_zero = COUNT(ABS(freq_in%solver) < eps14)
    IF (num_zero > 1) &
      CALL errore(__FILE__, "only a single frequency may be smaller than 1e-14", num_zero)

    ! determine number of frequencies
    num_freq_sym = 2 * num_freq - num_zero
    ALLOCATE(freq_out(num_freq_sym))

    !
    ! sanity test
    !
    IF (PRESENT(array)) THEN

      ! frequency and array should habe same dimension
      IF (SIZE(array, 3) /= num_freq_sym) THEN
        CALL errore(__FILE__, "array and frequency mesh inconsistent", 1)
      END IF

    END IF

    !
    ! copy array elements corresponding to nonzero frequency
    !
    ifreq_sym = num_freq
    DO ifreq = 1, num_freq
      !
      ! copy frequency
      freq_out(ifreq) = freq_in%solver(ifreq)
      !
      IF (ABS(freq_out(ifreq)) >= eps14) THEN
        !
        ! duplicate nonzero elements
        ifreq_sym = ifreq_sym + 1
        freq_out(ifreq_sym) = -freq_in%solver(ifreq)
        IF (PRESENT(array)) array(:, :, ifreq_sym) = array(:, :, ifreq)
        !
      END IF
      !
    END DO ! ifreq

  END SUBROUTINE freqbins_symm

  !> initialize the frequency bins type
  SUBROUTINE freqbins_init(this)

    USE kinds, ONLY: dp

    !> the frequency bins type which is initialized
    CLASS(freqbins_type), INTENT(INOUT) :: this

    !> temporary array for the symmetrized mesh
    COMPLEX(dp), ALLOCATABLE :: freq(:)

    !!
    !! 1. determines the symmetrized frequency grid (if applicable)
    !!
    ! evaluate the size of the symmetrized mesh and store it
    CALL freqbins_symm(this, freq)
    this%num_freq_sym = SIZE(freq)

  END SUBROUTINE freqbins_init

  !> symmetrize the frequency
  !!
  !! if no symmetry or even symmetry is used, return the input frequency
  !! if square symmetry is used, return the input frequency squared
  FUNCTION freqbins_symmetrize(this, freq) RESULT (res)

    USE kinds, ONLY: dp

    !> the frequency bins type defining the symmetry
    CLASS(freqbins_type), INTENT(IN) :: this

    !> the frequency to be symmetrized
    COMPLEX(dp), INTENT(IN) :: freq

    !> resulting symmetrized frequency
    COMPLEX(dp) :: res

    ! for square symmetry square the frequency
    IF (this%freq_symm_coul == square_symmetry) THEN
      res = freq**2

    ! otherwise return input frequency
    ELSE
      res = freq

    END IF

  END FUNCTION freqbins_symmetrize

END MODULE freqbins_module

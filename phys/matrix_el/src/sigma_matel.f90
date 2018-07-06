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
SUBROUTINE sigma_matel(ik0, grid, freq)

  USE buffers,              ONLY : get_buffer, close_buffer
  USE buiol,                ONLY : buiol_check_unit
  USE cell_base,            ONLY : tpiba2
  USE constants,            ONLY : RYTOEV, eps14
  USE control_gw,           ONLY : do_imag
  USE control_lr,           ONLY : lgamma
  USE disp,                 ONLY : xk_kpoints
  USE ener,                 ONLY : ef
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE freq_gw,              ONLY : nwsigma, nwsigwin
  USE freqbins_module,      ONLY : freqbins_type
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecw,                ONLY : ecutwfc
  USE gwsigma,              ONLY : nbnd_sig, corr_conv, exch_conv, ecutsco, ecutsex
  USE io_files,             ONLY : diropn 
  USE io_global,            ONLY : stdout, meta_ionode
  USE kinds_gw,             ONLY : i8b
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, lgauss
  USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_images,            ONLY : my_image_id
  USE mp_world,             ONLY : mpime
  USE output_mod,           ONLY : filsigx, filsigc
  USE qpoint,               ONLY : npwq
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE sigma_expect_mod,     ONLY : sigma_expect, sigma_expect_file
  USE sigma_grid_module,    ONLY : sigma_grid_type
  USE timing_module,        ONLY : time_matel
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma, lrsex, iunsex
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin

IMPLICIT NONE

  !> the FFT grid used
  TYPE(sigma_grid_type), INTENT(IN) :: grid

  !> the frequency grid used
  TYPE(freqbins_type),   INTENT(IN) :: freq

  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)
  COMPLEX(DP)               :: psic(dffts%nnr), vpsi(ngm)
  COMPLEX(DP)               :: ZdoTC, sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                               sigma_band_x(nbnd_sig, nbnd_sig, 1), vxc(nbnd_sig,nbnd_sig)
  REAL(DP)                  :: vtxc, etxc
  INTEGER                   :: igk(npwx), ikq
  INTEGER                   :: ig, ibnd, jbnd, ipol, ik0, ir
  INTEGER                   :: ng
  INTEGER                   :: sigma_c_ngm, sigma_x_ngm
  LOGICAL                   :: exst, opnd

  !> used to check nonzero file size
  INTEGER(i8b) file_size

  !> energy of the highest occupied state
  REAL(dp) ehomo

  !> energy of the lowest unoccupied state
  REAL(dp) elumo

  !> the chemical potential
  REAL(dp) mu

  !> real constant of 0.5
  REAL(dp),    PARAMETER :: half = 0.5_dp

  !> complex constant of 0
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

  CALL start_clock(time_matel)

  IF ( .NOT. meta_ionode ) RETURN

  nbnd   = nbnd_sig 
  lgamma = .true.

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

  ! check for gamma point
  IF(ALL(ABS(xk_kpoints(:,ik0)) <= eps14)) THEN
     ikq = 1
  ELSE
     ikq = 2
  END IF

  WRITE(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ik0, (xk_kpoints(ipol,ik0) , ipol = 1, 3)
  WRITE(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)

  ! set matrix elements to 0
  vxc          = 0
  sigma_band_x = 0
  sigma_band_c = 0

  WRITE(1000+mpime,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)

  ! create map to G ordering at current k-point
  CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                npw, igk, g2kin )
  npwq = npw
  ! read wave functions at current k-point
  CALL get_buffer (evc, lrwfc, iuwfc, ikq)
  evc(npw + 1:, :) = zero

  !
  ! expectation value of V_xc
  !

  ! generate v_xc(r) in real space:
  v%of_r = 0
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )

  WRITE(1000+mpime, '("Taking Matels.")')
  WRITE(1000+mpime, '("Taking NPWQ.", i4)') npwq
  WRITE(1000+mpime, '("my_image_id", i4)') my_image_id

  ! loop over all bands
  DO jbnd = 1, nbnd_sig

    ! extract right wave function
    psic = 0
    DO ig = 1, npwq
      ! map G-vectors according to smooth igk
      psic(dffts%nl(igk(ig))) = evc(ig, jbnd)
    END DO

    ! FFT wave function to real space
    CALL invfft ('Wave', psic(:), dffts)

    ! product of wave function and potential
    DO ir = 1, dfftp%nnr
      psic (ir) = psic(ir) * v%of_r (ir,1)
    END DO

    ! FFT back to reciprocal space
    CALL fwfft ('Wave', psic(:), dffts)

    ! reorder back to original order
    DO ig = 1, npwq
      vpsi(ig) = psic(dffts%nl(igk(ig)))
    END DO
 
    ! multiply with left wave function 
    DO ibnd = 1, nbnd_sig
      vxc(ibnd,jbnd) = ZdoTC (npwq, evc (1, ibnd), 1, vpsi, 1)
    END DO

  END DO ! bands

  WRITE(1000+mpime, '(4x,"VXC (eV)")')
  WRITE(1000+mpime, '(8(1x,f7.3))') real(vxc(:,:))*RYTOEV
  WRITE(1000+mpime, '("Max number Plane Waves WFC ", i4)') npwx
  WRITE(1000+mpime, '("Sigma_Ex Matrix Element")') 

  !
  ! expectation value of Sigma_x
  !
  ! open file containing exchange part of sigma
  INQUIRE( UNIT=iunsex, OPENED=opnd )
  IF (.NOT. opnd) CALL diropn( iunsex, filsigx, lrsex, exst )

  ! sanity check
  IF (ABS(exch_conv - ecutsex) < eps14 .OR. &
      ABS(exch_conv) < eps14) THEN
    sigma_x_ngm = grid%exch_fft%ngm
  ELSE IF((exch_conv < ecutsex) .AND. (exch_conv > 0.0)) THEN
    DO ng = 1, ngm
       IF ( gl( igtongl (ng) ) <= (exch_conv/tpiba2)) sigma_x_ngm = ng
    END DO
  ELSE
    CALL errore("sigma_matel", "Exch Conv must be greater than zero and less than ecutsex", 1)
  END IF

  ! check file size is not zero
  INQUIRE(UNIT = iunsex, SIZE = file_size)

  ! evaluate matrix elements for exchange
  IF (file_size /= 0) CALL sigma_expect_file(iunsex,ik0,evc,sigma_x_ngm,igk,sigma_band_x)
  WRITE(1000+mpime,*) 
  WRITE(1000+mpime,'(4x,"sigma_x (eV)")')
  WRITE(1000+mpime,'(8(1x,f7.3))') real(sigma_band_x)*RYTOEV
  WRITE(1000+mpime,*) 
  WRITE(1000+mpime,'(8(1x,f7.3))') aimag(sigma_band_x)*RYTOEV

  !
  ! expectation value of Sigma_c:
  !
  WRITE(1000+mpime,*) 
  WRITE(1000+mpime, '("sigma_c matrix element")') 

  ! open file containing correlation part of sigma
  INQUIRE( UNIT=iunsigma, OPENED=opnd )
  IF (.NOT. opnd) CALL diropn( iunsigma, filsigc, lrsigma, exst )

  ! For convergence tests corr_conv can be set at input lower than ecutsco.
  ! This allows you to calculate the correlation energy at lower energy cutoffs
  IF (ABS(corr_conv - ecutsco) < eps14) THEN
    sigma_c_ngm = grid%corr_fft%ngm
  ELSE IF(corr_conv < ecutsco .AND. corr_conv > 0.0) THEN
    DO ng = 1, ngm
      IF (gl( igtongl (ng) ) <= (corr_conv/tpiba2)) sigma_c_ngm = ng
    END DO
  ELSE
    CALL errore("sigma_matel", "Corr Conv must be greater than zero and less than ecutsco", 1)
  END IF

  WRITE(1000+mpime, *)
  WRITE(1000+mpime, '(5x, "G-Vects CORR_CONV:")')
  WRITE(1000+mpime, '(5x, f6.2, i5)') corr_conv, sigma_c_ngm
  WRITE(1000+mpime, *)

  ! check file size is not zero
  INQUIRE(UNIT = iunsigma, SIZE = file_size)

  ! evaluate expectation value of wave function
  IF (file_size /= 0) CALL sigma_expect_file(iunsigma,ik0,evc,sigma_c_ngm,igk,sigma_band_c,nwsigma)
  WRITE (1000+mpime,'("Finished Sigma_c")')

  !
  ! analytic continuation from imaginary frequencies to real ones
  ! 
  IF (do_imag) THEN
    ! We can set arbitrary \Sigma(\omega) energy windows with analytic continuation:
    ALLOCATE (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))
    ! print selfenergy on the imaginary axis.
    CALL print_matel_im(ikq, vxc(1,1), sigma_band_x(1,1,1), sigma_band_c(1,1,1), AIMAG(freq%sigma), nwsigma)
    ! do analytic continuation and print selfenergy on the real axis.
    sigma_band_con(:,:,:) = zero
    CALL sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), mu, AIMAG(freq%sigma), freq%window, nwsigwin)
    CALL print_matel(ikq, vxc(1,1), sigma_band_x(1,1,1), sigma_band_con(1,1,1), freq%window, nwsigwin)
   deallocate(sigma_band_con)
  ELSE
    ! print sigma on real axis
    CALL print_matel(ikq, vxc(1,1), sigma_band_x(1,1,1), sigma_band_c(1,1,1), REAL(freq%sigma) + mu, nwsigma)
  END IF

  CALL stop_clock(time_matel)

END SUBROUTINE sigma_matel

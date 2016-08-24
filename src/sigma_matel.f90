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
SUBROUTINE sigma_matel (ik0, freq)

  USE buffers,              ONLY : get_buffer, close_buffer
  USE buiol,                ONLY : buiol_check_unit
  USE cell_base,            ONLY : tpiba2
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE control_gw,           ONLY : lgamma, do_imag, output
  USE disp,                 ONLY : xk_kpoints
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE fft_custom,           ONLY : fft_cus, set_custom_grid, ggent, gvec_init
  USE freq_gw,              ONLY : nwsigma, nwsigwin
  USE freqbins_module,      ONLY : freqbins_type
  USE gvecs,                ONLY : nls
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecw,                ONLY : ecutwfc
  USE gwsigma,              ONLY : sigma_x_st, sigma_c_st, nbnd_sig, corr_conv, exch_conv, &
                                   sigma_band_exg, gcutcorr
  USE io_files,             ONLY : diropn 
  USE io_global,            ONLY : stdout, meta_ionode
  USE kinds,                ONLY : DP
  USE kinds_gw,             ONLY : i8b
  USE klist,                ONLY : xk
  USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_images,            ONLY : my_image_id
  USE mp_world,             ONLY : mpime
  USE output_mod,           ONLY : filsigx, filsigc
  USE qpoint,               ONLY : npwq
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE sigma_expect_mod,     ONLY : sigma_expect, sigma_expect_file
  USE timing_module,        ONLY : time_matel
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma, lrsex, iunsex
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin

IMPLICIT NONE

  TYPE(freqbins_type), INTENT(IN) :: freq

  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:)
  COMPLEX(DP)               :: psic(dffts%nnr), vpsi(ngm)
  COMPLEX(DP)               :: ZdoTC, sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                               sigma_band_x(nbnd_sig, nbnd_sig, 1), vxc(nbnd_sig,nbnd_sig)
  REAL(DP)                  :: vtxc, etxc
  INTEGER                   :: igk(npwx), ikq
  INTEGER                   :: ig, iw, ibnd, jbnd, ipol, ik0, ir
  INTEGER                   :: ng
  INTEGER                   :: sigma_c_ngm, sigma_x_ngm
  LOGICAL                   :: exst, opnd

  CALL start_clock(time_matel)

  IF ( .NOT. meta_ionode ) RETURN

  nbnd   = nbnd_sig 
  lgamma = .true.

  ! check for gamma point
  IF( ALL(xk_kpoints(:,ik0) == 0.0) ) THEN
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
      psic ( nls (igk(ig)) ) = evc(ig, jbnd)
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
      vpsi(ig) = psic(nls(igk(ig)))
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
  IF ((exch_conv == sigma_x_st%ecutt) .OR. (exch_conv == 0.0)) THEN
      sigma_x_ngm = sigma_x_st%ngmt
  ELSE IF((exch_conv < sigma_x_st%ecutt) .AND. (exch_conv > 0.0)) THEN
    DO ng = 1, ngm
       IF ( gl( igtongl (ng) ) <= (exch_conv/tpiba2)) sigma_x_ngm = ng
    END DO
  ELSE
    CALL errore("sigma_matel", "Exch Conv must be greater than zero and less than ecut_sco", 1)
  END IF

  ! evaluate matrix elements for exchange
  CALL sigma_expect_file(iunsex,ik0,evc,sigma_x_ngm,igk,sigma_band_x)
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
  IF (corr_conv == sigma_c_st%ecutt) THEN
    sigma_c_ngm = gcutcorr
  ELSE IF(corr_conv < sigma_c_st%ecutt .AND. corr_conv > 0.0) THEN
    DO ng = 1, ngm
      IF (gl( igtongl (ng) ) <= (corr_conv/tpiba2)) sigma_c_ngm = ng
    END DO
  ELSE
    CALL errore("sigma_matel", "Corr Conv must be greater than zero and less than ecut_sco", 1)
  END IF

  WRITE(1000+mpime, *)
  WRITE(1000+mpime, '(5x, "G-Vects CORR_CONV:")')
  WRITE(1000+mpime, '(5x, f6.2, i5)') corr_conv, sigma_c_ngm
  WRITE(1000+mpime, *)

  ! evaluate expectation value of wave function
  CALL sigma_expect_file(iunsigma,ik0,evc,sigma_c_ngm,igk,sigma_band_c,nwsigma)
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
    sigma_band_con(:,:,:) = dcmplx(0.0d0, 0.d0)
    CALL sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), AIMAG(freq%sigma), freq%window, nwsigwin)
    CALL print_matel(ikq, vxc(1,1), sigma_band_x(1,1,1), sigma_band_con(1,1,1), freq%window, nwsigwin)
   deallocate(sigma_band_con)
  ELSE
    ! print sigma on real axis
    CALL print_matel(ikq, vxc(1,1), sigma_band_x(1,1,1), sigma_band_c(1,1,1), REAL(freq%sigma), nwsigma)
  END IF
  IF (allocated(sigma_band_exg)) DEALLOCATE(sigma_band_exg)

  CALL stop_clock(time_matel)

END SUBROUTINE sigma_matel

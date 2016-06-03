!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
!
!Modified by Henry Lambert to use v^{2D}(r'-r), i.e.
!2D truncation of bare Coulomb interaction.
!
!-----------------------------------------------------------------------
subroutine dv_of_drho (mode, dvscf, flag)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the self consistent potential
  !     due to the perturbation.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2, fpi, eps8
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : nl, ngm, g,nlm, gstart
  USE cell_base, ONLY : alat, tpiba2, at, tpiba
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,     ONLY : dft_is_gradient
  USE scf,       ONLY : rho, rho_core
  USE eqv,       ONLY : dmuxc
  USE nlcc_gw,   ONLY : nlcc_any
  USE qpoint,    ONLY : xq
  USE gc_lr,     ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE control_gw, ONLY : lrpa, truncation
  USE control_flags, only : gamma_only
  USE disp,          ONLY : nq1, nq2, nq3, iq1, iq2, iq3
  USE truncation_module, ONLY : truncate, no_truncation, spherical_truncation

  implicit none

  integer :: mode
  ! input: the mode to do

  complex(DP), intent(inout):: dvscf (dfftp%nnr, nspin_mag)
  ! input: the change of the charge,
  ! output: change of the potential

  logical :: flag
  ! input: if true add core charge

  integer :: ir, is, is1, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors

  real(DP) :: qg2, fac, qxy, qz
  real(DP) :: zcut, spal
  ! the modulus of (q+G)^2
  ! the structure factor
  real(DP) :: q_G(3)
  ! the vector q + G

  complex(DP), allocatable :: dvaux (:,:), drhoc (:)
  !  the change of the core charge
  complex(DP), allocatable :: dvhart (:,:) !required in gamma_only

  call start_clock ('dv_of_drho')
  allocate (dvaux( dfftp%nnr,  nspin_mag))
  dvaux (:,:) = (0.d0, 0.d0)
  if (flag) allocate (drhoc( dfftp%nnr))    
  !
  ! the exchange-correlation contribution is computed in real space
  !
  if (lrpa) goto 111
  fac = 1.d0 / DBLE (nspin_lsda)
  if (nlcc_any.and.flag) then
  !extra contribution from frozen core charge:
     do is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
       enddo
  endif

  do is = 1, nspin_mag
     do is1 = 1, nspin_mag
        do ir = 1, dfftp%nnr
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        enddo
     enddo
  enddo
  !
  ! add gradient correction to xc, NB: if nlcc is true we need to add here
  ! its contribution. grho contains already the core charge
  !
  if ( dft_is_gradient() ) call dgradcorr &
       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf, dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat, dvaux)
  !Here they subtract out the core charge.
  if (nlcc_any.and.flag) then
     do is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
     enddo
  endif


111 continue
  !
  ! copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  if (nspin_mag == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2) 
  end if
  !
  CALL fwfft ('Dense', dvscf(:,1), dfftp)
  !
  ! hartree contribution is computed in reciprocal space
  !
  IF (gamma_only) THEN
    ALLOCATE(dvhart(dfftp%nnr,nspin_mag))
    dvhart(:,:) = (0.d0,0.d0)
    DO is = 1, nspin_lsda
      DO ig = 1, ngm
        q_G = (xq + g(:,ig)) * tpiba
        IF (truncation == spherical_truncation) THEN
          dvhart(nl(ig),is) = truncate(no_truncation, q_G) * dvscf(nl(ig),1)
        ELSE
          dvhart(nl(ig),is) = truncate(truncation, q_G) * dvscf(nl(ig),1)
        END IF
        dvhart(nlm(ig),is) = conjg(dvhart(nl(ig),is))
      END DO
      !
      !  and transformed back to real space
      !
      CALL invfft('Dense', dvhart(:,is), dfftp)
    END DO
    !
    ! at the end the two contributes are added
    dvscf  = dvaux  + dvhart
    DEALLOCATE(dvhart)
  ELSE
    DO is = 1, nspin_lsda
      CALL fwfft ('Dense', dvaux (:, is), dfftp)
      DO ig = 1, ngm
        q_G = (xq + g(:,ig)) * tpiba
        IF (truncation == spherical_truncation) THEN
          dvaux(nl(ig),is) = dvaux(nl(ig),is) + truncate(no_truncation, q_G) * dvscf(nl(ig),1)
        ELSE
          dvaux(nl(ig),is) = dvaux(nl(ig),is) + truncate(truncation, q_G) * dvscf(nl(ig),1)
        END IF
      END DO
      CALL invfft ('Dense', dvaux (:, is), dfftp)
    END DO
    ! at the end the two contributes are added
    dvscf (:,:) = dvaux (:,:)
  END IF
  !
  if (flag) deallocate (drhoc)
  deallocate (dvaux)
  call stop_clock ('dv_of_drho')
  return
end subroutine dv_of_drho

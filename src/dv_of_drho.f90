!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dv_of_drho (mode, dvscf, flag)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the self consistent potential
  !     due to the perturbation.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2, fpi
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
  USE gc_gw,     ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE control_gw, ONLY : lrpa, trunc_2d
  USE control_flags, only : gamma_only
  !OBM: gamma_only is disregarded for phonon calculations, TDDFPT purposes only

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
  if (gamma_only) then
    allocate(dvhart(dfftp%nnr,nspin_mag))
    dvhart(:,:) = (0.d0,0.d0)
    do is = 1, nspin_lsda
      do ig = 1, ngm
         qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
         if (qg2 > 1.d-8) then
            dvhart(nl(ig),is) = e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
            dvhart(nlm(ig),is) = conjg(dvhart(nl(ig),is))
         endif
      enddo
      !
      !  and transformed back to real space
      !
      CALL invfft('Dense', dvhart(:,is), dfftp)
    enddo
    !
    ! at the end the two contributes are added
    dvscf  = dvaux  + dvhart
    deallocate(dvhart)
  else
    do is = 1, nspin_lsda
       CALL fwfft ('Dense', dvaux (:, is), dfftp)
       IF(.not.trunc_2d) then
         do ig = 1, ngm
            qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
            if (qg2 > 1.d-8) then
                dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                                 e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
            endif
         enddo
       ELSE
           zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
           DO ig = 1, ngm
                qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
             IF (qg2 > 1.d-8) then
                qxy  = sqrt((g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2)
                qz   = sqrt((g(3,ig)+xq(3))**2)
                spal = 1.0d0 - EXP(-tpiba*qxy*zcut)*cos(tpiba*qz*zcut)
                dvaux(nl(ig), is) = dvaux(nl(ig), is) + dvscf(nl(ig), 1)*dcmplx(fpi*e2/(tpiba2*qg2)*spal, 0.0d0)
             ENDIF 
           ENDDO
          !DO ig = 1, ngm
          !   qg2 = (g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2 + (g(3,ig)+xqloc(3))**2
          !   qxy  = sqrt((g(1,ig) + xqloc(1))**2 + (g(2,ig) + xqloc(2))**2)
          !   qz   = sqrt((g(3,ig) + xqloc(3))**2)
          !   IF(qxy.gt.eps8) then
          !      spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
          !      DO igp = 1, sigma_c_st%ngmt
          !         dvaux(nl(ig)) = dvaux(nl(ig)) + dvscf(nl(ig),1)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
          !      ENDDO
          !   ELSE IF(qxy.lt.eps8.and.qz.gt.eps8) then
          !      spal = 1.0d0 + EXP(-tpiba*qxy*zcut)*((qz/qxy)*sin(tpiba*qz*zcut) - cos(tpiba*qz*zcut))
          !      DO igp = 1, sigma_c_st%ngmt
          !         dvaux(nl(ig)) = dvaux(nl(ig)) + dvscf(nl(ig),1)*dcmplx((e2*fpi/(tpiba2*qg2))*spal, 0.0d0)
          !      ENDDO
          !   ELSE
          !      dvaux(nl(ig)) = 0.0d0
          !   ENDIF
          !ENDDO
       ENDIF
       CALL invfft ('Dense', dvaux (:, is), dfftp)
    enddo
    ! at the end the two contributes are added
    dvscf (:,:) = dvaux (:,:)
  endif
  !
  if (flag) deallocate (drhoc)
  deallocate (dvaux)
  call stop_clock ('dv_of_drho')
  return
end subroutine dv_of_drho

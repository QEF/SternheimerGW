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
SUBROUTINE sigma_matel (ik0)
  USE io_global,            ONLY : stdout, meta_ionode
  USE io_files,             ONLY : prefix
  USE buffers,              ONLY : get_buffer, close_buffer
  USE kinds,                ONLY : DP
  USE kinds_gw,             ONLY : i8b
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecs,                ONLY : nls
  USE constants,            ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,              ONLY : nwsigma, wsigma, wsig_wind_min, wsig_wind_max, deltaws, nwsigwin
  USE klist,                ONLY : xk, nkstot
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin
  USE gvecw,                ONLY : ecutwfc
  USE qpoint,               ONLY : npwq
  USE units_gw,             ONLY : iunsigma, iuwfc, lrwfc, lrsigma, lrsex, iunsex
  USE control_gw,           ONLY : lgamma, do_imag, do_serial, do_sigma_exxG, tmp_dir_coul
  USE wavefunctions_module, ONLY : evc
  USE gwsigma,              ONLY : sigma_x_st, sigma_c_st, nbnd_sig, corr_conv, exch_conv, &
                                   sigma_band_exg, gcutcorr
  USE disp,                 ONLY : xk_kpoints
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE cell_base,            ONLY : tpiba2, at, alat
  USE buiol,                ONLY : buiol_check_unit
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE fft_custom,           ONLY : fft_cus, set_custom_grid, ggent, gvec_init
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_world,             ONLY : mpime
  USE mp_images,            ONLY : my_image_id, inter_image_comm
  USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
  USE reorder_mod,          ONLY : reorder, create_map
  USE sigma_expect_mod,     ONLY : sigma_expect

IMPLICIT NONE
  complex(DP), allocatable  :: sigma_band_con(:,:,:)
  complex(DP), allocatable  :: sigma_g_ex(:,:)
  complex(DP)               ::   czero
  complex(DP)               ::   psic(dffts%nnr), vpsi(ngm)
  complex(DP)               ::   ZdoTC, sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                                 sigma_band_ex(nbnd_sig, nbnd_sig), vxc(nbnd_sig,nbnd_sig)
  complex(DP), allocatable  ::   sigma(:,:,:)
  real(DP), allocatable     ::   wsigwin(:)
  real(DP)                  ::   w_ryd(nwsigma)
  real(DP)                  ::   one, zcut
  real(DP)    :: vtxc, etxc
  real(DP)    :: zero(3)
  integer, allocatable      ::   map(:) 
  integer, allocatable      ::   igkq_ig(:) 
  integer, allocatable      ::   igkq_tmp(:) 
  integer                   ::   ikq, ikq_head
  integer                   ::   ig, iw, ibnd, jbnd, ios, ipol, ik0, ir
  integer                   ::   ng
  integer     :: sigma_c_ngm, sigma_x_ngm
  integer     :: kpoolid(nkstot), iqrec1(nkstot)
  integer(i8b) :: unf_recl
  logical, external :: eqvect
  logical :: found_k
  character(len=256) :: tempfile, filename
  real(DP), parameter :: eps=1.e-5_dp

#define DIRECT_IO_FACTOR 8

  ALLOCATE( map(npwx) )
  allocate (igkq_tmp(npwx))
  allocate (igkq_ig(npwx))

  one   = 1.0d0 
  czero = (0.0d0, 0.0d0)
  w_ryd = wsigma(:nwsigma)/RYTOEV
  nbnd = nbnd_sig 
  zero(:) = 0.d0
  lgamma=.true.
  ikq = 1
  found_k = .false.

  if((xk_kpoints(1,ik0).eq.0.0).and.(xk_kpoints(2,ik0).eq.0.0).and.(xk_kpoints(3,ik0).eq.0.0))then 
     ikq_head = 1
  else
     ikq_head = 2
  endif

  ikq = ikq_head
  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ik0, (xk_kpoints(ipol,ik0) , ipol = 1, 3)
  write(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)
  kpoolid = 0
  iqrec1  = 0

!All pools need access to sigma file now:
  filename = trim(prefix)//"."//"sigma1"
  tempfile = trim(tmp_dir_coul) // trim(filename)
  unf_recl = DIRECT_IO_FACTOR * int(lrsigma, kind=kind(unf_recl))
  open(iunsigma, file = trim(adjustl(tempfile)), iostat = ios, &
  form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
  write(1000+mpime,*) tempfile, ios

  if(.not. do_sigma_exxG) then
     filename = trim(prefix)//"."//"sigma_ex1"
     tempfile = trim(tmp_dir_coul) // trim(filename)
     unf_recl = DIRECT_IO_FACTOR * int(lrsex, kind=kind(unf_recl))
     open(iunsex, file = trim(adjustl(tempfile)), iostat = ios, &
     form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
     write(1000+mpime,*) tempfile, ios
  endif

!ONLY THE POOL WITH THIS KPOINT CALCULATES THE CORRECT MATRIX ELEMENT.
  vxc(:,:) = czero
  sigma_band_ex (:, :) = czero
  sigma_band_c (:,:,:) = czero

  if (meta_ionode) THEN
      write(1000+mpime,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, (xk(ipol,ikq) , ipol = 1, 3)
      CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                    npw, igk, g2kin )
      npwq = npw
      call get_buffer (evc, lrwfc, iuwfc, ikq)
      zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
! generate v_xc(r) in real space:
      v%of_r(:,:) = (0.0d0)
      CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
      write(1000+mpime, '("Taking Matels.")')
      write(1000+mpime, '("Taking NPWQ.", i4)') npwq
      write(1000+mpime, '("my_image_id", i4)') my_image_id
      do jbnd = 1, nbnd_sig
         psic = czero
         do ig = 1, npwq
            psic ( nls (igk(ig)) ) = evc(ig, jbnd)
         enddo
!Need to do this fft according to igkq arrays and switching between serial/parallel routines. 
         CALL invfft ('Wave', psic(:), dffts)
         do ir = 1, dfftp%nnr
            psic (ir) = psic(ir) * v%of_r (ir,1)
         enddo
         CALL fwfft ('Wave', psic(:), dffts)
         do ig = 1, npwq
            vpsi(ig) = psic(nls(igk(ig)))
         enddo
         do ibnd = 1, nbnd_sig
            vxc(ibnd,jbnd) = ZdoTC (npwq, evc (1, ibnd), 1, vpsi, 1)
         enddo
      enddo
      write(1000+mpime, '(4x,"VXC (eV)")')
      write(1000+mpime, '(8(1x,f7.3))') real(vxc(:,:))*RYTOEV
      write(1000+mpime, '("Max number Plane Waves WFC ", i4)') npwx
      write(1000+mpime, '("Sigma_Ex Matrix Element")') 

    IF( .NOT. do_sigma_exxG) THEN

      ALLOCATE (sigma_g_ex (sigma_x_st%ngmt, sigma_x_st%ngmt))

      IF ((exch_conv == sigma_x_st%ecutt) .OR. (exch_conv == 0.0)) THEN
          sigma_x_ngm = sigma_x_st%ngmt
      ELSE IF((exch_conv < sigma_x_st%ecutt) .AND. (exch_conv > 0.0)) THEN
        DO ng = 1, ngm
           IF ( gl( igtongl (ng) ) <= (exch_conv/tpiba2)) sigma_x_ngm = ng
        END DO
      ELSE
        CALL errore("sigma_matel", "Exch Conv must be greater than zero and less than ecut_sco", 1)
      END IF

      ! read sigma_x from file
      ios = 0 
      READ( UNIT = iunsex, REC = ik0, IOSTAT = ios ) sigma_g_ex

      ! if reading failed
      IF (ios /= 0) THEN
        WRITE(1000+mpime, '(5x, "Could not read Sigma_X file. Have you calculated it?")') 
        sigma_band_ex = 0.0

      ! evaluate expectation value of wave function
      ELSE

        ! create map for reordering
        map = create_map(igk,sigma_x_st%ngmt)

        ! reorder evc array so that it is compatible with current igk
        CALL reorder(evc,map)

        ! evaluate the expectation value
        sigma_band_ex = sigma_expect( sigma_g_ex, evc(:sigma_x_st%ngmt,:) )

        WRITE(1000+mpime,*) 
        WRITE(1000+mpime,'(4x,"sigma_ex (eV)")')
        WRITE(1000+mpime,'(8(1x,f7.3))') real(sigma_band_ex(:,:))*RYTOEV
        WRITE(1000+mpime,*) 
        WRITE(1000+mpime,'(8(1x,f7.3))') aimag(sigma_band_ex(:,:))*RYTOEV

      END IF

      DEALLOCATE(sigma_g_ex)

    ELSE ! do_sigma_exxG = .true.

      DO ibnd = 1, nbnd_sig
        sigma_band_ex(ibnd,ibnd) = sigma_band_exg(ibnd)
      END DO

    END IF

!MATRIX ELEMENTS OF SIGMA_C:
    WRITE(1000+mpime,*) 
    WRITE(1000+mpime, '("sigma_c matrix element")') 
    ALLOCATE (sigma(gcutcorr, gcutcorr,nwsigma)) 

!For convergence tests corr_conv can be set at input lower than ecutsco.
!This allows you to calculate the correlation energy at lower energy cutoffs
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

    ! read sigma from file
    sigma = 0.0

    IF(do_serial) THEN
      DO iw = 1, nwsigma
        CALL davcio (sigma(:,:, iw), lrsigma, iunsigma, iw, -1)
      END DO
    ELSE
      ! let's us avoid crash if we haven't calculated one of these things yet:
      READ( UNIT = iunsigma, REC = ik0, IOSTAT = ios ) sigma
      WRITE(1000+mpime, *) ios

      ! if reading sigma failed
      IF (ios /= 0) THEN
        WRITE(1000+mpime, '("Could not read Sigma_C file. Have you calculated it?")')
        sigma_band_c (:,:,:) = czero

      ! evaluate expectation value of wave function
      ELSE

        ! create map for reordering
        map = create_map(igk,sigma_c_ngm)

        ! read evc
        CALL get_buffer (evc, lrwfc, iuwfc, ikq)

        ! reorder evc array so that it is compatible with current igk
        CALL reorder(evc,map)

        ! evaluate the expectation value
        sigma_band_c = sigma_expect( sigma, evc(:sigma_c_ngm,:) )

        DEALLOCATE (sigma)
        write (1000+mpime,'("Finished Sigma_c")')

      END IF

    END IF
!Need to broadcast from the current pool to all the nodes

  END If !on pool with K-point
  CALL mp_barrier(inter_pool_comm)

!Now first pool should always have
!the kpoint we are looking for.
  if(meta_ionode) THEN
     if(do_imag) then 
!We can set arbitrary \Sigma(\omega) energy windows with analytic continuation:
        allocate (wsigwin(nwsigwin))
        do iw = 1, nwsigwin
            wsigwin(iw) = wsig_wind_min + (wsig_wind_max-wsig_wind_min)/float(nwsigwin-1)*float(iw-1)
        enddo
        allocate (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))
!print selfenergy on the imaginary axis.
        call print_matel_im(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
        !call print_matel_im(2, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
!do analytic continuation and print selfenergy on the real axis.
        sigma_band_con(:,:,:) = dcmplx(0.0d0, 0.d0)
        call sigma_pade(sigma_band_c(1,1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
        call print_matel(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_con(1,1,1), wsigwin(1), nwsigwin)
     else
        call print_matel(ikq_head, vxc(1,1), sigma_band_ex(1,1), sigma_band_c(1,1,1), wsigma(1), nwsigma)
     endif
  endif
  if(allocated(sigma_band_con)) deallocate(sigma_band_con)
  if(allocated(igkq_tmp)) deallocate(igkq_tmp)
  if(allocated(igkq_ig))  deallocate(igkq_ig)
  if(allocated(sigma_band_exg)) deallocate(sigma_band_exg)
  IF(ALLOCATED(map)) DEALLOCATE(map)

call mp_barrier(inter_pool_comm)
call mp_barrier(inter_image_comm)
return
end SUBROUTINE sigma_matel

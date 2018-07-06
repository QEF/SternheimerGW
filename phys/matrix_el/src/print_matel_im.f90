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
subroutine print_matel_im(ikq, vxc, sigma_band_ex, sigma_band_c, w_ryd, nwsigma) 

  USE constants,            ONLY : RYTOEV, pi
  USE gwcom,                ONLY : output
  USE gwsigma,              ONLY : nbnd_sig
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk
  USE pp_output_mod,        ONLY : pp_output, pp_output_xml
  USE wvfct,                ONLY : nbnd, et

implicit none

INTEGER                   ::   nwsigma
REAL(DP)                  ::   wsigma(nwsigma) 
INTEGER                   ::   iw, ibnd, jbnd
REAL(DP)                  ::   w_ryd(nwsigma)
REAL(DP)                  ::   one
COMPLEX(DP)               ::   sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                               sigma_band_ex(nbnd_sig, nbnd_sig), vxc(nbnd_sig,nbnd_sig)
REAL(DP)                  ::   resig_diag(nwsigma,nbnd_sig), imsig_diag(nwsigma,nbnd_sig),&
                               et_qp(nbnd_sig), a_diag(nwsigma,nbnd_sig)
REAL(DP)                  ::   dresig_diag(nwsigma,nbnd_sig), vxc_tr, vxc_diag(nbnd_sig),&
                               sigma_ex_tr, sigma_ex_diag(nbnd_sig)

REAL(DP)                  ::   resig_diag_tr(nwsigma), imsig_diag_tr(nwsigma), a_diag_tr(nwsigma),&
                               et_qp_tr, z_tr, z(nbnd_sig)

COMPLEX(DP)               ::   czero

INTEGER                   ::   iman, nman, ndeg(nbnd_sig), ideg, ikq
LOGICAL                   ::   single_line

     one   = 1.0d0 
     czero = (0.0d0, 0.0d0)
     nbnd  = nbnd_sig 
     wsigma(:) = w_ryd(:) * RYTOEV

     do ibnd = 1, nbnd_sig
        do iw = 1, nwsigma
           resig_diag (iw,ibnd) = real( sigma_band_c(ibnd, ibnd, iw))
           dresig_diag (iw,ibnd) = resig_diag (iw,ibnd) + real(sigma_band_ex(ibnd,ibnd)) - real( vxc(ibnd,ibnd) )
           imsig_diag (iw,ibnd) = aimag ( sigma_band_c (ibnd, ibnd, iw) )
           a_diag (iw,ibnd) = one/pi * abs ( imsig_diag (iw,ibnd) ) / &
               ( abs( w_ryd(iw) - et(ibnd, ikq) - ( resig_diag (iw,ibnd) + sigma_band_ex(ibnd, ibnd) - vxc(ibnd,ibnd) ) )**2.d0 &
               + abs ( imsig_diag (iw,ibnd) )**2.d0 )
        enddo
        call qp_eigval ( nwsigma, w_ryd, dresig_diag(1,ibnd), et(ibnd,ikq), et_qp (ibnd), z(ibnd) )
     enddo

 ! Now take the trace (get rid of phase arbitrariness of the wfs)
 ! (alternative and more appropriate: calculate non-diagonal, elements of
 ! degenerate subspaces and diagonalize)
 ! count degenerate manifolds and degeneracy...
  nman = 1
  ndeg = 1

  do ibnd = 2, nbnd_sig
     if ( abs( et (ibnd, ikq) - et (ibnd-1, ikq)  ) .lt. 1.d-5 ) then
        ndeg (nman) = ndeg(nman) + 1
     else
        nman = nman + 1
     endif
  enddo

  write(6,'(" Manifolds")')
  write (stdout, *) nman, (ndeg (iman) ,iman=1,nman)
  write(6,*)
  
  ! ...and take the trace over the manifold
  
  ibnd = 0
  jbnd = 0

  do iman = 1, nman
    resig_diag_tr = 0.d0
    imsig_diag_tr = 0.d0
    a_diag_tr = 0.d0
    et_qp_tr = 0.d0
    z_tr = 0.d0
    vxc_tr = 0.d0
    sigma_ex_tr = 0.0d0

    do ideg = 1, ndeg(iman)
       ibnd = ibnd + 1
       resig_diag_tr = resig_diag_tr + resig_diag (:,ibnd)
       imsig_diag_tr = imsig_diag_tr + imsig_diag (:,ibnd)
       a_diag_tr = a_diag_tr + a_diag (:,ibnd)
       et_qp_tr = et_qp_tr + et_qp (ibnd)
       z_tr = z_tr + z (ibnd)
       vxc_tr = vxc_tr + real(vxc(ibnd,ibnd))
       sigma_ex_tr = sigma_ex_tr + real(sigma_band_ex(ibnd,ibnd))
    enddo

    do ideg = 1, ndeg(iman)
      jbnd = jbnd + 1
      resig_diag (:,jbnd) = resig_diag_tr / float( ndeg(iman) )
      imsig_diag (:,jbnd) = imsig_diag_tr / float( ndeg(iman) )
      a_diag (:,jbnd) = a_diag_tr / float( ndeg(iman) )
      et_qp (jbnd) = et_qp_tr / float( ndeg(iman) )
      z (jbnd) = z_tr / float( ndeg(iman) )
      vxc_diag (jbnd) = vxc_tr / float( ndeg(iman) )
      sigma_ex_diag(jbnd) = sigma_ex_tr/float(ndeg(iman))
    enddo
  enddo

  if(nbnd_sig.le.8) single_line=.true.
  if(nbnd_sig.gt.8) single_line=.false.

  write(stdout,*)
  write(stdout,'("REsigma")')
  do iw = 1, nwsigma
    if(single_line) then
       write(stdout,'(9f14.7)') wsigma(iw), (RYTOEV*resig_diag (iw,ibnd), ibnd=1,nbnd_sig)
    else
       write(stdout,'(9f14.7)', advance='no') wsigma(iw), (RYTOEV*resig_diag (iw,ibnd), ibnd=1,8)
    endif

    if(nbnd_sig.gt.8) then
    do ideg = 9, nbnd_sig, 8 
       if(ideg+7.lt.nbnd_sig) write(stdout,9005,advance='no') (RYTOEV*resig_diag (iw,ideg:ideg+7)) 
       if(ideg+7.ge.nbnd_sig) write(stdout,9005) (RYTOEV*resig_diag (iw,ideg:nbnd_sig)) 
    enddo
    endif
  enddo

  write(stdout,*)
  write(stdout,'("IMsigma")')
  do iw = 1, nwsigma
     if(single_line) then
        write(stdout,'(9f15.8)') wsigma(iw), (RYTOEV*imsig_diag (iw,ibnd), ibnd=1,nbnd_sig)
     else
        write(stdout,'(9f15.8)', advance='no') wsigma(iw), (RYTOEV*imsig_diag (iw,ibnd), ibnd=1,8)
     endif
     if(nbnd_sig.gt.8) then
     do ideg = 9, nbnd_sig, 8
        if(ideg+7.lt.nbnd_sig) write(stdout, 9005, advance='no') (RYTOEV*imsig_diag (iw,ibnd), ibnd=ideg,ideg+7)
        if(ideg+7.ge.nbnd_sig) write(stdout, 9005) (RYTOEV*imsig_diag (iw,ibnd), ibnd=ideg,nbnd_sig)
     enddo
     endif
  enddo

  write(stdout,*)
  write(stdout,'("ASpec")')
  do iw = 1, nwsigma
     if(single_line) then
        write(stdout,'(9f15.8)') wsigma(iw), (a_diag (iw,ibnd)/RYTOEV,ibnd=1,nbnd_sig)
     else
        write(stdout,'(9f15.8)',advance='no') wsigma(iw), (a_diag (iw,ibnd)/RYTOEV, ibnd=1,8)
     endif

     if(nbnd_sig.gt.8) then
     do ideg = 9, nbnd_sig, 8
        if(ideg+7.lt.nbnd_sig) write(stdout, 9005,advance='no') (a_diag (iw,ibnd)/RYTOEV, ibnd=ideg,ideg+7)
        if(ideg+7.ge.nbnd_sig) write(stdout, 9005) (a_diag (iw,ibnd)/RYTOEV, ibnd=ideg,nbnd_sig)
     enddo
     endif
  enddo
  write(stdout,*)

  !
  ! print output files according to user requirement
  !
  CALL pp_output_xml(output%pp_re_corr_iw, ikq, xk(:,ikq), wsigma, resig_diag * RYTOEV)
  CALL pp_output_xml(output%pp_im_corr_iw, ikq, xk(:,ikq), wsigma, imsig_diag * RYTOEV)
  CALL pp_output_xml(output%pp_spec_iw,    ikq, xk(:,ikq), wsigma, a_diag / RYTOEV)

  9005 format(8(1x,f14.7))
RETURN
end subroutine print_matel_im

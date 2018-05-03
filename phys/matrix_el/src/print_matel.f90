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
subroutine print_matel(ikq, vxc, sigma_band_ex, sigma_band_c, w_ryd, nwsigma) 

use cell_base,            only : at
use constants,            only : RYTOEV, pi
use gwcom,                only : output
use gwsigma,              only : nbnd_sig
use io_global,            only : stdout
use kinds,                only : DP
use klist,                only : xk
use pp_output_mod,        only : pp_output, pp_output_xml
use wvfct,                only : nbnd, et

implicit none

integer                   ::   nwsigma
complex(DP)               ::   sigma_band_c(nbnd_sig, nbnd_sig, nwsigma),&
                               sigma_band_ex(nbnd_sig, nbnd_sig), vxc(nbnd_sig,nbnd_sig)
complex(DP)               ::   czero
real(DP)                  ::   wsigma(nwsigma)
real(DP)                  ::   w_ryd(nwsigma), xkcryst(3)
real(DP)                  ::   one
real(DP)                  ::   resig_diag(nwsigma,nbnd_sig), imsig_diag(nwsigma,nbnd_sig),&
                               et_qp(nbnd_sig), a_diag(nwsigma,nbnd_sig)
real(DP)                  ::   dresig_diag(nwsigma,nbnd_sig), vxc_tr, vxc_diag(nbnd_sig),&
                               sigma_ex_tr, sigma_ex_diag(nbnd_sig)
real(DP)                  ::   resig_diag_tr(nwsigma), imsig_diag_tr(nwsigma), a_diag_tr(nwsigma),&
                               et_qp_tr, z_tr, z(nbnd_sig)
integer                   ::   iw, ibnd, jbnd
integer                   ::   iman, nman, ndeg(nbnd_sig), ideg, ikq
logical                   ::   single_line

     one   = 1.0d0 
     czero = (0.0d0, 0.0d0)
     nbnd  = nbnd_sig 
     wsigma(:) = w_ryd(:) * RYTOEV

     do ibnd = 1, nbnd_sig
        do iw = 1, nwsigma
           resig_diag (iw,ibnd) = real(sigma_band_c(ibnd, ibnd, iw))
           dresig_diag (iw,ibnd) = resig_diag (iw,ibnd) + real(sigma_band_ex(ibnd,ibnd)) - real( vxc(ibnd,ibnd) )
           imsig_diag (iw,ibnd) = aimag ( sigma_band_c (ibnd, ibnd, iw) )
           a_diag (iw,ibnd) = one/pi * abs ( imsig_diag (iw,ibnd) ) / &
               ( abs( w_ryd(iw) - et(ibnd, ikq) - (resig_diag (iw,ibnd) + sigma_band_ex(ibnd, ibnd) - vxc(ibnd,ibnd) ) )**2.d0 &
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

  xkcryst(:) = xk(:, ikq)
  call cryst_to_cart(1, xkcryst, at, -1)
  if(single_line) then
     write(stdout,'(/4x,"LDA eigenval (eV)", 8(1x,f7.2))')  et(1:nbnd_sig, ikq)*RYTOEV
  else
     !write(stdout,'(/4x,"LDA eigenval (eV)", 8(1x,f7.2))', advance='no')  et(1:8, ikq)*RYTOEV
     write(stdout,'(/4x,"LDA eigenval (eV)", 8(1x,f7.2))')  et(1:8, ikq)*RYTOEV
  endif
  if(nbnd_sig.gt.8) then
  do ideg = 9, nbnd_sig, 8 
     !if(ideg+7.lt.nbnd_sig) write(stdout,9000, advance='no')  et(ideg:ideg+7, ikq)*RYTOEV
     if(ideg+7.lt.nbnd_sig) write(stdout,9000)  et(ideg:ideg+7, ikq)*RYTOEV
     if(ideg+7.ge.nbnd_sig) write(stdout,9000)  et(ideg:nbnd_sig, ikq)*RYTOEV
  enddo
  endif

  write(stdout,'(/4x,"GWKpoint cart :", 3(1x,f8.4))') xk(:,ikq)
  write(stdout,'(/4x,"GWKpoint cryst:", 3(1x,f8.4))') xkcryst(:)

  if(single_line) then
     write(stdout, '(4x,"GW qp energy (eV)",8(1x,f7.2))')  et_qp(1:nbnd_sig)*RYTOEV
  else
     !write(stdout, '(4x,"GW qp energy (eV)",8(1x,f7.2))', advance='no')  et_qp(1:8)*RYTOEV
     write(stdout, '(/4x,"GW qp energy (eV)",8(1x,f7.2))')  et_qp(1:8)*RYTOEV
  endif

  if(nbnd_sig.gt.8) then
  do ideg = 9, nbnd_sig, 8 
     !if(ideg+7.lt.nbnd_sig) write(stdout,9000, advance='no') et_qp(ideg:ideg+7)*RYTOEV
     if(ideg+7.lt.nbnd_sig) write(stdout,9000) et_qp(ideg:ideg+7)*RYTOEV
     if(ideg+7.ge.nbnd_sig) write(stdout,9000)  et_qp(ideg:nbnd_sig)*RYTOEV
  enddo
  endif

  if(single_line) then
     write(stdout,'(4x,"Vxc expt val (eV)",8(1x,f7.2))')  vxc_diag(1:nbnd_sig)*RYTOEV
  else
     !write(stdout,'(4x,"Vxc expt val (eV)",8(1x,f7.2))', advance='no')  vxc_diag(1:8)*RYTOEV
     write(stdout,'(/4x,"Vxc expt val (eV)",8(1x,f7.2))')  vxc_diag(1:8)*RYTOEV
  endif

  if(nbnd_sig.gt.8) then
  do ideg = 9, nbnd_sig, 8 
     !if(ideg+7.lt.nbnd_sig) write(stdout,9000, advance='no')  vxc_diag(ideg:ideg+7)*RYTOEV
     if(ideg+7.lt.nbnd_sig) write(stdout,9000)  vxc_diag(ideg:ideg+7)*RYTOEV
     if(ideg+7.ge.nbnd_sig) write(stdout,9000)  vxc_diag(ideg:nbnd_sig)*RYTOEV
  enddo
  endif

  if(single_line) then
     write(stdout,'(4x,"Sigma_ex val (eV)",8(1x,f7.2))')  sigma_ex_diag(1:nbnd_sig)*RYTOEV
  else
     !write(stdout,'(4x,"Sigma_ex val (eV)",8(1x,f7.2))', advance='no')  sigma_ex_diag(1:8)*RYTOEV
     write(stdout,'(/4x,"Sigma_ex val (eV)",8(1x,f7.2))')  sigma_ex_diag(1:8)*RYTOEV
  endif

  if(nbnd_sig.gt.8) then
  do ideg = 9, nbnd_sig, 8 
     !if(ideg+7.lt.nbnd_sig) write(stdout,9000, advance='no')  sigma_ex_diag(ideg:ideg+7)*RYTOEV
     if(ideg+7.lt.nbnd_sig) write(stdout,9000)  sigma_ex_diag(ideg:ideg+7)*RYTOEV
     if(ideg+7.ge.nbnd_sig) write(stdout,9000)  sigma_ex_diag(ideg:nbnd_sig)*RYTOEV
  enddo
  endif

  if(single_line) then
     write(stdout,'(4x,"QP renorm",8x, 8(1x,f7.2))')  z(1:nbnd_sig)
  else
     !write(stdout,'(4x,"QP renorm",8x, 8(1x,f7.2))', advance='no')  z(1:8)
     write(stdout,'(/4x,"QP renorm",8x, 8(1x,f7.2))')  z(1:8)
  endif

  if(nbnd_sig.gt.8) then
  do ideg = 9, nbnd_sig, 8 
     !if(ideg+7.lt.nbnd_sig) write(stdout,9000, advance='no')  z(ideg:ideg+7)
     if(ideg+7.lt.nbnd_sig) write(stdout,9000)  z(ideg:ideg+7)
     if(ideg+7.ge.nbnd_sig) write(stdout,9000)  z(ideg:nbnd_sig)
  enddo
  endif

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
  CALL pp_output(output%pp_dft,      xk(:,ikq), et(:,ikq) * RYTOEV)
  CALL pp_output(output%pp_gw,       xk(:,ikq), et_qp * RYTOEV)
  CALL pp_output(output%pp_vxc,      xk(:,ikq), vxc_diag * RYTOEV)
  CALL pp_output(output%pp_exchange, xk(:,ikq), sigma_ex_diag * RYTOEV)
  CALL pp_output(output%pp_renorm,   xk(:,ikq), z)
  CALL pp_output_xml(output%pp_re_corr, ikq, xk(:,ikq), wsigma, resig_diag * RYTOEV)
  CALL pp_output_xml(output%pp_im_corr, ikq, xk(:,ikq), wsigma, imsig_diag * RYTOEV)
  CALL pp_output_xml(output%pp_spec,    ikq, xk(:,ikq), wsigma, a_diag / RYTOEV)

  9000 format(21x, 8(1x,f7.2))
  9005 format(8(1x,f14.7))
RETURN
end subroutine print_matel

!----------------------------------------------------------------
  SUBROUTINE  qp_eigval ( nw, w, sig, et, et_qp, z )
!----------------------------------------------------------------
!
  use constants,     only : RYTOEV
  use kinds,         only : DP

  IMPLICIT NONE

  integer :: nw, iw, iw1, iw2
  real(DP) :: w(nw), sig(nw), et, et_qp, dw, w1, w2, sig_et, sig1, sig2, z, sig_der, one
  
  one = 1.0d0
  dw = w(2)-w(1)

 if ((et.lt.w(1)+dw).or.(et.gt.w(nw)-dw)) then
 !call errore ('qp_eigval','original eigenvalues outside the frequency range of the self-energy',1)
  write(6,'("original eigenvalues outside the frequency range of the self-energy")')
  write(6,'(3f7.2)') RYTOEV*et, RYTOEV*(w(1)+dw), RYTOEV*(w(nw) - dw)
 return
 endif

  iw = 1
  iw1 = 1
  iw2 = 1
  do while ((iw.lt.nw).and.(w(iw).lt.et))
    iw = iw + 1
    iw1 = iw-1
    iw2 = iw
  enddo

  w1 = w(iw1)
  w2 = w(iw2)
  sig1 = sig(iw1)
  sig2 = sig(iw2)

  sig_et = sig1 + ( sig2 - sig1 ) * (et-w1) / (w2-w1)
  sig_der = ( sig2 - sig1 ) / ( w2 - w1 )
  z = one / ( one - sig_der)

! temporary - until I do not have Vxc
  et_qp = et + z * sig_et
  END SUBROUTINE qp_eigval
!----------------------------------------------------------------

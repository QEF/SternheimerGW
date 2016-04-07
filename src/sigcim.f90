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
subroutine sigma_c_im(ik0) 
!G TIMES W PRODUCT sigma_correlation_imaginary frequency.
  use kinds,         only : DP
  use kinds_gw,      only : i8b
  use io_global,     only : stdout, ionode_id, ionode, meta_ionode
  use io_files,      only : iunigk, prefix, tmp_dir
  use lsda_mod,      only : nspin
  use constants,     only : e2, fpi, RYTOEV, tpi, eps8, pi
  use disp,          only : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  use control_gw,    only : lgamma, eta, godbyneeds, padecont, cohsex, modielec, &
                            trunc_2d, tmp_dir_coul, high_io
  use klist,         only : wk, xk, nkstot, nks, lgauss
  use wvfct,         only : nbnd, npw, npwx, igk, g2kin
  use eqv,           only : evq
  use freq_gw,       only : fpol, fiu, nfs, nfsmax, &
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                            wgreen, wsigma, wsigmamin, wsigmamax, &
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw, &
                            w0pmw, wgtcoul
  use units_gw,      only : iuncoul, iungreen, iunsigma, lrsigma,&
                            lrcoul, lrgrn, iuwfc, lrwfc
  use qpoint,        only : xq, npwq, igkq, nksq, ikks, ikqs
  use gvect,         only : g, ngm, nl
  use cell_base,     only : tpiba2, tpiba, omega, alat, at,bg
  use symm_base,     only : nsym, s, time_reversal, t_rev, ftau, invs, nrot, invsym
  use lr_symm_base,  only : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  use wavefunctions_module, only : evc
  use control_flags,        only : noinv
  use ener,          only : ef
  use gwsigma,       only : sigma_c_st, gcutcorr
  use mp_global,     only : mp_global_end
  use mp_world,      only : nproc, mpime
  use mp_images,     only : nimage, my_image_id, intra_image_comm,   &
                            me_image, nproc_image, inter_image_comm
  use mp,            only : mp_sum, mp_barrier, mp_bcast
  use mp_pools,      only : inter_pool_comm

  implicit none

  complex(DP)         :: ci, czero
  complex(DP)         :: phase
  complex(DP)         :: cprefac, dz
!Sigma arrays
  complex (DP), allocatable :: sigma(:,:,:)
  complex (DP), allocatable :: sigma_g(:,:,:)
!Pade arrays
  complex(DP), allocatable :: z(:), u(:), a(:)
!W arrays 
  complex(DP), allocatable :: scrcoul_g (:,:,:)
  complex(DP), allocatable :: scrcoul_pade_g (:,:)
  complex(DP), allocatable :: scrcoul(:,:)
!G arrays:
  complex(DP), allocatable :: greenf_g(:,:,:), greenfr(:,:)
  complex(DP), allocatable  ::  eigv(:,:)
!v array
!complex(DP), allocatable ::  barcoul(:,:), barcoulr(:,:), barcoul_R(:,:)
  real(DP) :: qg2, qg, qxy, qz
  real(DP) :: w_ryd(nwcoul), w_rydsig(nwsigma)
  real(DP) :: xq_ibk(3), xq_ibz(3)
!q-vector of coulomb potential:
  real(DP) :: xq_coul(3)
  real(DP) :: rcut, spal
!CHECK FOR NAN's
  real(DP)    :: ar, ai
!For dirac delta fxn.
  real(DP)    :: dirac, x, support, zcut
  real(DP)    :: ehomo, elumo, mu
  real (DP)   :: xk1(3), aq(3), xk1_old(3)
  real(DP)    :: nsymm1
  real(DP)    :: wgt(nsym), xk_un(3,nsym)
  real(DP)            :: wgtcoulry(nwcoul)
  real(DP), parameter :: eps=1.e-5_dp
  real(DP), allocatable :: xqtr(:,:), wqtr(:)
!FREQUENCY GRIDS/COUNTERS
  integer, allocatable  :: gmapsym(:,:)
  integer  :: iwim, iw, ikq
  integer  :: iw0, iw0mw, iw0pw
  integer  :: iqstart, iqstop
!COUNTERS
  integer :: ig, igp, irr, icounter, ir, irp
  integer :: iqs, nkr
  integer :: iq, ipol, iqrec
  integer :: ikmq, ik0, ik, nkpool
  integer :: rec0, ios
  integer :: counter, ierr
  integer :: inversym, screening
!SYMMETRY
  integer     :: isym, jsym, isymop, nig0
!For G^NA
  integer     :: igkq_ig(npwx) 
  integer     :: igkq_tmp(npwx) 
  integer     :: ibnd
  integer     :: iw0start, iw0stop
  integer     :: nnr, nqstr
  integer(i8b):: unf_recl
!For running PWSCF need some variables 
  logical             :: pade_catch
  logical             :: found_q, trev
  logical             :: limq, found
  logical, external   :: eqvect
  logical             :: invq_k(nsym)
!File related:
  character(len=256)  :: tempfile, filename
!Star of q:
  real(DP) :: sxq(3,48), xqs(3,48)
  integer  :: imq, isq(48), nqstar
  integer  :: nsq(48), i, nsymrot, iq1

#define DIRECT_IO_FACTOR 8 
! iG(W-v)
   allocate ( scrcoul_g       (gcutcorr, gcutcorr, nfs)     )
   allocate ( scrcoul_pade_g  (gcutcorr, gcutcorr)          )
   allocate ( greenf_g        (gcutcorr, gcutcorr, 2*nwcoul))
   allocate ( sigma_g         (gcutcorr, gcutcorr, nwsigma))

!Writing these to disk
!These go on the big grid...
!Technically only need gmapsym up to gcutcorr or ngmgrn...
   allocate ( gmapsym  (ngm, nsym)   )
   allocate ( eigv     (ngm, nsym)   )
   allocate (z(nfs), a(nfs), u(nfs))

   nnr = sigma_c_st%dfftt%nnr
   wgtcoulry(:) = wgtcoul(:)/RYTOEV
   w_ryd(:)    = wcoul(:)/RYTOEV
   w_rydsig(:) = wsigma(:)/RYTOEV

   write(stdout,'(/4x,"Direct product GW for k0(",i3," ) = (", 3f12.7," )")') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
   write(stdout,'(4x, "nfs, ", i4, " nwsigma, ", i4)') nfs, nwsigma
   write(stdout,'(4x, "nrsco, ", i4, " nfs, ", i4)') sigma_c_st%dfftt%nnr, nfs
   write(stdout,'(4x, "nsym, nsymq, nsymbrav ", 3i4)') nsym, nsymq, nrot 
   write(stdout,'(4x, "gcutcorr", i4 )') gcutcorr
   allocate (wqtr(nq1*nq2*nq3))
   allocate (xqtr(3, nq1*nq2*nq3))
!Generate grid of kpoints with inversion symmetry.
   !    call kpoint_grid(nsym, .false., .false., s, t_rev,& 
   !                     bg, nq1*nq2*nq3, 0,0,0, nq1, nq2, nq3, nqstr, xqtr, wqtr)
   !else
   !    nqstr = nqs
   !endif
   write(stdout,'(4x, "num of q points in convolution: ", i4 )') nqs
   !do iq = 1, nqstr
   !   write(stdout, '(5x,i3, 4f14.9)') iq, xqtr(1,iq), xqtr(2,iq), xqtr(3,iq), wqtr(iq)
   !end do
   !do iq = 1, nqstr
   !   xq(:) = xqtr(:, iq)
   !   write(stdout, '(5x,i3, 4f14.9)') iq, xqtr(1,iq), xqtr(2,iq), xqtr(3,iq), wqtr(iq)
   !   trev = .false.
   !   call find_trev(xq, s, invs, iqtr,isym,trev)
   !   write(stdout,*) trev
   !   write(stdout, '(5x,i3, 4f14.9)') iqtr, x_q(1,iqtr), x_q(2,iqtr), x_q(3,iqtr), wq(iqtr)
   !enddo
   ci = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
!2D Truncation
   zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
   call start_clock('sigmac')
   call gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)
!Set appropriate weights for points in the brillouin zone.
!Weights of all the k-points are in odd positions in list.
!nksq is number of k points not including k+q.
!Every processor needs access to the files: _gw0si.coul1 and _gw0si.green1
   call mp_barrier(inter_image_comm)
#ifdef __PARA
!OPEN coulomb file (only written to by head node).
   filename = trim(prefix)//"."//"coul1"
   tempfile = trim(tmp_dir_coul) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
   open(iuncoul, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
#endif
  call para_img(nwsigma, iw0start, iw0stop)
  write(stdout, '(5x, "nwsigma ",i4, " iw0start ", i4, " iw0stop ", i4)') nwsigma, iw0start, iw0stop
  write(stdout,'("Starting Frequency Integration")')
  call get_homo_lumo (ehomo, elumo)
  if(.not.lgauss) then
    mu = ehomo + 0.5d0*(elumo-ehomo)
  else
    mu = ef
  endif
  call mp_barrier(inter_pool_comm)
  call mp_bcast(mu, ionode_id ,inter_pool_comm)
  call mp_barrier(inter_pool_comm)
  nsymm1 = 1.0d0/dble(nsym)
  call para_pool(nqs,iqstart,iqstop)
  do iq = iqstart, iqstop
     trev = .false.
     xq(:) = x_q(:,iq)
     CALL star_q(xq(1), at, bg, nsym, s, invs, nqstar, sxq, isq, nsq, imq, .false. )
   !  write( 1000+mpime, * )
   !  write( 1000+mpime, '(5x,a,i4)') 'Number of q in the star = ', nqstar
   !  write( 1000+mpime, '(5x,a)') 'List of q in the star:'
   !  write( 1000+mpime, '(7x,i4,i4,i4,3f14.9)') (iq1, nsq(iq1), isq(iq1), (sxq(i,iq1), i=1,3), iq1=1,nqstar)
   !  write( 1000+mpime, '(7x,i4,i4)') (iq1, isq(iq1), iq1=1,nsym)
     cprefac = wq(iq)*dcmplx(-1.0d0, 0.0d0)/tpi
     do iq1 = 1, nqstar
        scrcoul_g(:,:,:) = dcmplx(0.0d0, 0.0d0)
        if(.not.modielec) call davcio(scrcoul_g, lrcoul, iuncoul, iq, -1)
        nsymrot=0
        do isym=1,nsym
           if (isq(isym) == iq1) then
               nsymrot=nsymrot+1
               if (nsymrot == 1) isymop=isym
           endif
        enddo
        !do iw = 1, nfs
        !   scrcoul_pade_g(:,:) = scrcoul_g(:,:,iw)
        !   do ig = 1, gcutcorr
        !      do igp = 1, gcutcorr
!\eps^-1_{Sq}(\G,\G') = \eps^{-1}{q}(s^{-1}G,S^{-1}G')
                 !scrcoul_g(gmapsym(ig,invs(isymop)), gmapsym(igp,invs(isymop)), iw) = scrcoul_pade_g(ig,igp)
                 !scrcoul_g(gmapsym(ig,invs(isymop)), gmapsym(igp,invs(isymop)), iw) = scrcoul_pade_g(ig,igp)
        !         scrcoul_g(ig, igp, iw) = scrcoul_pade_g(gmapsym(ig,isymop),gmapsym(igp,isymop))
        !      enddo
        !   enddo
        !enddo
        if(nsymrot == 0) then
           call errore('dfile_star','no symmetry relates q at star(q)',iq)
        endif
       !call rotate(xq(1), aq, s, nsym, invs(isymop))
        xk1 = xk_kpoints(:,ik0) - sxq(:,iq1) 
        call coulpade(scrcoul_g(1,1,1), xq(:))
        call green_linsys_shift_im(greenf_g(1,1,1), xk1(1), 1, mu, 2*nwcoul)
        nig0    = 1
        if(iw0stop-iw0start+1.gt.0) THEN
           do iw0 = iw0start, iw0stop
              do iw = 1, nwcoul
                 dz =  dcmplx(float(nsq(iq1))*nsymm1*wgtcoulry(iw),0.0d0)*cprefac
                 call construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), abs(w_ryd(iw)-w_rydsig(iw0)) )
                 call sigprod(isymop, dz, scrcoul_pade_g(1,1), greenf_g(1,1,iw), sigma_g(1,1,iw0), gmapsym(1,1))
                 call construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), abs(w_rydsig(iw0)+w_ryd(iw)))
                 call sigprod(isymop, dz, scrcoul_pade_g(1,1), greenf_g(1,1,iw+nwcoul), sigma_g(1,1,iw0), gmapsym(1,1))
              enddo ! on frequency convolution over w'
           enddo ! on iw0
        endif
     enddo ! star_q
  enddo ! iq
  deallocate ( eigv           )
  deallocate ( gmapsym        )
  deallocate ( greenf_g       )
  deallocate ( scrcoul_pade_g )
  deallocate ( scrcoul_g      )
  deallocate ( z, a, u        )
#ifdef __PARA
  call mp_barrier(inter_pool_comm)
  call mp_sum(sigma_g, inter_pool_comm)
  call mp_barrier(inter_image_comm)
  call mp_sum(sigma_g, inter_image_comm)
#endif __PARA
!Now write Sigma in G space to file. 
  if (meta_ionode) THEN
      call davcio (sigma_g, lrsigma, iunsigma, ik0, 1)
      write(6,'(4x,"Sigma Written to File")')
      call stop_clock('sigmac')
  endif !ionode
  deallocate ( sigma_g  )
  call mp_barrier(inter_image_comm)
  return
end subroutine sigma_c_im

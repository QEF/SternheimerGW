  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
subroutine sigma_c_im_mod(ik0) 
!G TIMES W PRODUCT sigma_correlation_imaginary frequency.
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode_id, ionode, meta_ionode
  USE io_files,      ONLY : iunigk, prefix, tmp_dir
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec, trunc_2d, tmp_dir_coul
  USE klist,         ONLY : wk, xk, nkstot, nks, lgauss
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE eqv,           ONLY : evq, eprec
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, &
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                            wgreen, wsigma, wsigmamin, wsigmamax, &
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw, &
                            w0pmw, wgtcoul
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gvect,         ONLY : g, ngm, nl
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at,bg
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE modes,         ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : noinv
  USE ener,          ONLY : ef
  USE gwsigma,       ONLY : sigma_c_st, gcutcorr
  USE mp_global,     ONLY : mp_global_end
  USE mp_world,      ONLY : nproc, mpime
  USE mp_images,     ONLY : nimage, my_image_id, intra_image_comm,   &
                            me_image, nproc_image, inter_image_comm
  USE mp,            ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_pools,      ONLY : inter_pool_comm

  IMPLICIT NONE

  COMPLEX(DP)         :: ci, czero
  COMPLEX(DP)         :: phase
!Sigma arrays
  COMPLEX (DP), ALLOCATABLE :: sigma(:,:,:)
  COMPLEX (DP), ALLOCATABLE :: sigma_g(:,:,:)
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:,:), greenfr(:,:)
  COMPLEX(DP) :: cprefac, dz
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
!v array
!COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:), barcoul_R(:,:)
  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: w_ryd(nwcoul), w_rydsig(nwsigma)
  REAL(DP) :: xq_ibk(3), xq_ibz(3)
!q-vector of coulomb potential:
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: rcut, spal
!CHECK FOR NAN's
  REAL(DP)     :: ar, ai
!For dirac delta fxn.
  REAL(DP)     :: dirac, x, support, zcut
  REAL(DP) :: ehomo, elumo, mu
  REAL (DP)   :: xk1(3), aq(3), xk1_old(3)
  REAL(DP)    :: sxq(3,48), xqs(3,48)
  REAL(DP)    :: nsymm1
  REAL(DP)    :: wgt(nsym), xk_un(3,nsym)
  REAL(DP)            :: wgtcoulry(nwcoul)
  REAL(DP), parameter :: eps=1.e-5_dp
!FREQUENCY GRIDS/COUNTERS
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  INTEGER  :: iwim, iw, ikq
  INTEGER  :: iw0, iw0mw, iw0pw
  INTEGER  :: iqstart, iqstop
!COUNTERS
  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER :: iqs, nkr
  INTEGER :: iq, ipol, iqrec
  INTEGER :: ikmq, ik0, ik, nkpool
  INTEGER :: rec0, ios
  INTEGER :: counter, ierr
  INTEGER :: inversym, screening
!SYMMETRY
  INTEGER               :: isym, jsym, isymop, nig0
!For G^NA
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ss(3,3)
  INTEGER     :: ibnd
  INTEGER     :: iw0start, iw0stop
!Complete file name
  INTEGER     :: imq, isq(48), nqstar, nkpts
  INTEGER     :: i, ikstar
  INTEGER     :: ixk1, iqrec_old
  INTEGER     :: isym_k(nsym), nig0_k(nsym), iqrec_k(nsym)
  INTEGER     :: nnr
  INTEGER*8   :: unf_recl
!For running PWSCF need some variables 
  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limq, inv_q, found
  LOGICAL, EXTERNAL   :: eqvect
  LOGICAL             :: invq_k(nsym)
!File related:
  CHARACTER(len=256)  :: tempfile, filename

#define DIRECT_IO_FACTOR 8 
! iG(W-v)
   allocate ( greenf_g        (gcutcorr, gcutcorr, 2*nwcoul))
   allocate ( sigma_g         (gcutcorr, gcutcorr, nwsigma) )
   allocate ( scrcoul_g       (gcutcorr, gcutcorr, nfs)     )
   allocate ( scrcoul_pade_g  (gcutcorr, gcutcorr)          )

!These go on the big grid...
!Technically only need gmapsym up to gcutcorr or ngmgrn...
   allocate ( gmapsym  (ngm, nrot)   )
   allocate ( eigv     (ngm, nrot)   )
!This is a memory hog...
   allocate (z(nfs), a(nfs), u(nfs))

   nnr = sigma_c_st%dfftt%nnr
   wgtcoulry(:) = wgtcoul(:)/RYTOEV

   w_ryd(:) = wcoul(:)/RYTOEV
   w_rydsig(:) = wsigma(:)/RYTOEV

   write(stdout,'(/4x,"Direct product GW for k0(",i3," ) = (", 3f12.7," )")') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
   write(stdout,'(4x, "nfs, ", i4, " nwsigma, ", i4)') nfs, nwsigma
   write(stdout,'(4x, "nrsco, ", i4, " nfs, ", i4)') sigma_c_st%dfftt%nnr, nfs
   write(stdout,'(4x, "nsym, nsymq, nsymbrav ", 3i4)'), nsym, nsymq, nrot 
   write(stdout,'(4x, "gcutcorr", i4 )') gcutcorr
   ci = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
!2D Truncation
   zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
   CALL start_clock('sigmac')
   CALL gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
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
  CALL para_img(nwsigma, iw0start, iw0stop)
  write(stdout, '(5x, "nwsigma ",i4, " iw0start ", i4, " iw0stop ", i4)') nwsigma, iw0start, iw0stop
  write(stdout,'("Starting Frequency Integration")')
  CALL get_homo_lumo (ehomo, elumo)
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
  xk1_old(:) =  -400.0
!  write(1000+mpime,*),'iqstart, stop ', iqstart, iqstop
  do iq = iqstart, iqstop
     scrcoul_g(:,:,:) = dcmplx(0.0d0, 0.0d0)
     if(.not.modielec) CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, -1)
     cprefac = wq(iq)*dcmplx(-1.0d0, 0.0d0)/tpi
     xq(:) = x_q(:,iq)
     CALL coulpade(scrcoul_g(1,1,1), xq(1))
!write(1000+mpime,*), xq
     do isymop = 1, nsym
!write(1000+mpime,*), isymop
        CALL rotate(xq, aq, s, nsym, invs(isymop))
        xk1 = xk_kpoints(:,ik0) - aq(:)
!only recalculate as required
        if (.not.(abs(xk1(1)-xk1_old(1)).lt.eps .and.   &
                  abs((xk1(2)-xk1_old(2))).lt.eps .and. & 
                  abs((xk1(3)-xk1_old(3))).lt.eps)) THEN
                  CALL green_linsys_shift_im(greenf_g(1,1,1), xk1(1), 1, mu, 2*nwcoul)
        endif
        isym     = 1
        nig0     = 1
        inv_q   = .false.
        if(iw0stop-iw0start+1.gt.0) THEN
           do iw0 = iw0start, iw0stop
              do iw = 1, nwcoul
               call construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), abs(w_ryd(iw)-w_rydsig(iw0)))
               dz = dcmplx(nsymm1*wgtcoulry(iw),0.0d0)*cprefac
               call sigprod(isymop, dz, scrcoul_pade_g(1,1), greenf_g(1,1,iw), sigma_g(1,1,iw0), gmapsym(1,1))
               call construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), (w_rydsig(iw0)+w_ryd(iw)))
               call sigprod(isymop, dz, scrcoul_pade_g(1,1), greenf_g(1,1,iw+nwcoul), sigma_g(1,1,iw0), gmapsym(1,1))
              enddo ! on frequency convolution over w'
           enddo ! on iw0  
        endif
     enddo ! isymop
enddo!iq
deallocate ( eigv           )
deallocate ( gmapsym        )
deallocate ( greenf_g       )
deallocate ( scrcoul_pade_g )
deallocate ( scrcoul_g      )
deallocate ( z, a, u        )
#ifdef __PARA
 CALL mp_barrier(inter_pool_comm)
 CALL mp_sum(sigma_g, inter_pool_comm)
 CALL mp_barrier(inter_image_comm)
 CALL mp_sum(sigma_g, inter_image_comm)
#endif __PARA
 if (meta_ionode) THEN
!Now write Sigma in G space to file. 
   CALL davcio (sigma_g, lrsigma, iunsigma, ik0, 1)
   write(6,'(4x,"Sigma Written to File")')
   CALL stop_clock('sigmac')
 endif !ionode
 deallocate ( sigma_g  )
 CALL mp_barrier(inter_image_comm)
RETURN
end subroutine sigma_c_im_mod

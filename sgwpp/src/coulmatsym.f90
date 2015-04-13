! Copyright (C) 2004-2009 Andrea Benassi and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!------------------------------
SUBROUTINE coulmatsym()
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE cell_base,            ONLY : tpiba2, omega, at, alat, bg
  USE ener,                 ONLY : ef
  USE io_global,            ONLY : ionode, stdout
  USE start_k,              ONLY : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE klist,                ONLY : nks, nkstot, ngauss, degauss, xk, wk, nelec
  USE wavefunctions_module, ONLY : evc 
  USE wvfct,                ONLY : nbnd, et, npw, igk, npwx,  g2kin, ecutwfc
  USE lsda_mod,             ONLY : nspin
  USE ktetra,               ONLY : ntetra, tetra, ltetra
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE start_k,              ONLY : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE io_files,             ONLY : nwordwfc, iunwfc, diropn
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat
  USE symm_base,            ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE gvect,                ONLY : ngm, g, nl
  USE gvecs,                ONLY : nls, nlsm
  USE control_coulmat,      ONlY : degaussfs, nbndmin, debye_e, do_lind, ngcoul, do_diag
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_world,             ONLY : world_comm, mpime
  USE mp_pools,             ONLY : nproc_pool, me_pool, my_pool_id, inter_pool_comm, npool
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : inter_image_comm, intra_image_comm, &
                                   my_image_id, nimage, root_image

IMPLICIT NONE

  COMPLEX(DP) :: vcnknpkp
  REAL(DP)    :: enk, enpkp
  REAL(DP)    :: norm
  REAL(DP)    :: En, DOSofE(2), N0, phase
  REAL(DP)    :: degaussw0, w0g1, w0g2
  REAL(DP)    :: qg2, xq(3), xkp(3)
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  COMPLEX(DP) :: psink(dffts%nnr,nbnd), psinpkp(dffts%nnr,nbnd), psi_temp(dffts%nnr), fnknpkp(dffts%nnr)
  COMPLEX(DP) :: pwg0(dffts%nnr)
  INTEGER     :: ik, ikp, ibnd, jbnd, ig, igp, isymop, isym
  INTEGER     :: iq , iq1, i
  INTEGER     :: nkp, nkp_abs, ipool
!IO
  LOGICAL     :: exst
!to put in coul struct:
  REAL(DP)    :: mu, mustar, muloc
  REAL(DP)    :: munnp(2, 2)
  REAL(DP)    :: dosnnp(2)
  REAL(DP)    :: kcut
! SYMMMETRY
  LOGICAL     :: found_q, inv_q, minus_q
  INTEGER     :: iqrec, irotmq
  REAL (DP)   :: gi(3,48), gimq(3), aq(3)
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
  COMPLEX(DP), ALLOCATABLE :: vc(:,:)
  INTEGER :: ir, niG0
  INTEGER :: nqstart, nqstop
!For Star of q.
  REAL(DP) :: sxq(3,48), xqs(3,48)
  INTEGER  :: imq, isq(48), nqstar, nqs
  INTEGER  :: nsq(48)
  REAL(DP) :: xkp_loc(3), xk_loc(3)

  REAL(DP) :: ehomo, elomo, bandwidth
  REAL(DP) :: muk(nks)


  ALLOCATE ( gmapsym  (ngm, nsym)   )
  ALLOCATE ( eigv     (ngm, nsym)   )

  gmapsym(:,:) = 0
  eigv(:,:) = 0

  CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

  ALLOCATE(vc(ngcoul,ngcoul))
!  First calculate dosef: 
!  dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
!  N(Ef) in the equation for lambda is the DOS per spin
!  dosef = dosef / two
!  CALCULATE Dos at Fermi Level gaussian:
  if(ltetra) then
   WRITE( stdout,'(/5x,"USING TETRAHEDRAL INTEGRATION")')
   CALL dos_t(et,nspin,nbnd, nks,ntetra,tetra, ef, DOSofE)
  else
  WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
       &        "ngauss,degauss=",i4,f12.6/)') ngauss, degauss
   CALL dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, ef, DOSofE)
  endif
! tetra=.true.
! use tetrahedron method:
! CALL dos_t(et, nspin, nbnd, nks, ntetra, tetra, ef, DOSofE)
! want density of states per spin.
  N0 = DOSofE(1)/2.d0
  if (degaussfs.eq.0.0d0) then
     degaussw0 = 0.05
  else
     degaussw0 = degaussfs
  endif

  En      = ef
  mu         = 0.0d0
  munnp(:,:) = 0.0d0
  dosnnp(:) = 0.0d0
  !kcut = 0.33333
  kcut = 0.12

  !write(stdout, * )  nbnd
  !write(stdout, '(5X, 6f12.5)' ) munnp(1:nbnd,1:nbnd)

  call parallelize(nks, nqstart, nqstop)

  write(1000+mpime, *) nqstart, nqstop
  !write(1000+mpime, *) xk(:, 1:nks)
  !write(1000+mpime, *) wk(1:nks)
!Again only collinear case
  ibnd = NINT( nelec ) / 2
  ehomo = MAXVAL( et(ibnd,  1:nkstot) )
  elomo = MINVAL( et(:,  1:nkstot))
  bandwidth = ef - elomo
!print*, "Ef Bandwidth: ", bandwidth*rytoev
  do iq = nqstart, nqstop
     write(1000+mpime, *) 
     write(1000+mpime, '(5x, "iq: ", i4)') iq
     write(1000+mpime, '(5x, "xq: ", 3f9.5)') xk(:, iq)
!Only need sym ops in star of q the rest are just nq*Vnkn'k'
     CALL star_q(xk(:,iq), at, bg, nsym, s, invs, nqs, sxq, isq, nsq, imq, .false. )
  !  IF (verbosity) THEN
     WRITE( 1000+mpime, * )
     WRITE( 1000+mpime, '(5x,a,i4)') 'Number of q in the star = ', nqs
     WRITE( 1000+mpime, '(5x,a)') 'List of q in the star:'
     WRITE( 1000+mpime, '(7x,i4,i4,i4,3f14.9)') (iq1, nsq(iq1), isq(iq1), (sxq(i,iq1), i=1,3), iq1=1,nqs)
  !ENDIF
     muloc = 0.0d0
     do ik = 1, nks
!    if (mod(ik,2).eq.0) print*, "xk:  ", ik
!    Load coulomb interaction at xq
        vc = (0.0d0, 0.0d0)
        call load_coul(vc, iq)
!do S = nsymops:
!need to modify this:
!On a gen'l grid ktokpmq should also return the symmetry operation which rotates
!xkp to the brillouin zone. then we keep track of that symmetry else and use it
!to invert one of kpoints in the IBZ. Keeping with the idea that we only need
!IBZk and IBZq.
        psink(:,:)   = (0.0d0,0.0d0)
        CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL get_buffer (evc, nwordwfc, iunwfc, ik)
        psink(nls(igk(1:npw)), :) = evc(1:npw,:)
        do ibnd = nbndmin, nbnd
          CALL invfft ('Wave', psink(:,ibnd), dffts)
!MODIFYING FOR STAR Q
          do isymop = 1, nsym
!        x_k' = x_k-S^-1*xq
            xq(:) = xk(:,iq)
            CALL rotate(xq, aq, s, nsym, invs(isymop))
            xkp(:) = xk(:, ik)-aq(:)
!        Find symmop what gives us k'
            inv_q=.false.
            call find_qG_ibz(xkp, s, iqrec, isym, nig0, found_q, inv_q)
            if(iqrec.gt.nks) call errore('coulmatsym','COULD NOT MAP k+q to IBZ',1)
            ikp = iqrec
            psinpkp(:,:) = (0.0d0,0.0d0)
            CALL gk_sort (xk(1,iqrec), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
            CALL get_buffer(evc, nwordwfc, iunwfc, ikp)
            psinpkp(nls(gmapsym(igk(1:npw), isym)),:) = evc(1:npw,:)
            if(nig0.gt.1) then
             pwg0(:) = dcmplx(0.0d0, 0.0d0)
             pwg0(nls(nig0)) = dcmplx(1.0d0, 0.0d0)
             CALL invfft('Wave', pwg0(:), dffts)
            endif
            do jbnd = nbndmin, nbnd
!calculate f_{nk,npkp}(\G)
             psi_temp = (0.0d0, 0.0d0)
             psi_temp = psinpkp(:,jbnd)
             CALL invfft ('Wave', psi_temp(:), dffts)
             fnknpkp = (0.0d0,0.0d0)
!calc psi_{\k-\q+\G}
             if(nig0.gt.1) then
               do ir = 1, dffts%nnr  
                  psi_temp(ir) = psi_temp(ir)*pwg0(ir) 
               enddo
             endif              
             do ir = 1, dffts%nnr  
                fnknpkp(ir) = fnknpkp(ir) + conjg(psi_temp(ir))*psink(ir,ibnd)
             enddo
             CALL fwfft ('Wave', fnknpkp(:), dffts)
!Eigenvalues
             enk = (et(ibnd, ik) - ef)
             enpkp = (et(jbnd, ikp) - ef)

             w0g1 = w0gauss ( enk / degaussw0, 0) / degaussw0
             w0g2 = w0gauss ( enpkp / degaussw0, 0) / degaussw0

!Again we want per spin.
             vcnknpkp = 0.0d0
             if(.not.do_lind) then
               if(.not.do_diag) then
                 do ig = 1, ngcoul 
                  do igp = 1, ngcoul
                     phase = eigv(ig,isymop)*conjg(eigv(igp,isymop))
                     vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,igp)*fnknpkp(nls(igp))*phase
                  enddo
                 enddo
               else
                 do ig = 1, ngcoul 
                     phase = eigv(ig,isymop)*conjg(eigv(ig,isymop))
                     vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,ig)*fnknpkp(nls(ig))*phase
                 enddo
               endif
             else
               do ig = 1, ngcoul 
                  phase = eigv(ig,isymop)*conjg(eigv(igp,isymop))
                  vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,ig)*fnknpkp(nls(ig))*phase
               enddo
             endif
!mu = mu + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)
!wk weight includes 2*spin index... so for Coulomb we kill that...
!and to get per spin we need to kill it in the wk factor.
             mu = mu + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
             muk(ik) = muk(ik) + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
!Need to separate the bands in k space
!for MgB2 to disambiguate the bands I think we only
!need to check that the k point lies within a certain distance of
!the gamma point. I think we only need the first corner of the Gamma
!point!
             xk_loc(:)   = xk(:,ik)
             xkp_loc(:)  = xk(:,ikp)

             CALL cryst_to_cart(1, xk_loc(:), at, -1)
             CALL cryst_to_cart(1, xkp_loc(:), at, -1)

             if(sqrt(xk_loc(1)**2 + xk_loc(2)**2) .lt. kcut ) then
                if ((sqrt((xkp_loc(1))**2 + (xkp_loc(2)))**2).lt.kcut) then
                   munnp(1, 1) = munnp(1,1) + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
                   if(ik.eq.1) dosnnp(1) = dosnnp(1) + w0g1*(wk(ik)/2.0)
                else
                   munnp(1, 2) = munnp(1,2) + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
!                   dosnnp(1) = dosnnp(1) + w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
                endif
             else
                if ((sqrt((xkp_loc(1))**2 + (xkp_loc(2)))**2).lt.kcut) then
                   munnp(2, 1) = munnp(2,1) + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
!                  dosnnp(2,1) = dosnnp(2,1) + w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
                else
                   munnp(2, 2) = munnp(2,2) + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
                   if(ik.eq.1) dosnnp(2) = dosnnp(2) + w0g1*(wk(ik)/2.0)
                endif
             endif
             muloc = muloc + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)
          enddo!jbnd
         enddo!isym
        enddo!ibnd
       enddo!ik
       write(1000+mpime, '(5x,"\mu(iq) " 1f12.5)') muloc*wk(1)/2.0
       write(1000+mpime, '(5X, 6f12.5)' ) munnp(1:2,1:2)
     enddo!iq
     CALL mp_sum(mu, inter_image_comm)!reduce over q points
     CALL mp_sum(munnp, inter_image_comm)!reduce over q points
!Factors not included when we calculate V^{c}_{nkn'k'}.

     mu = mu/(omega*nsym)
     munnp = munnp/(omega*nsym)

   write(stdout,*) nk1, nk2, nk3
   write(stdout,*) omega
   write(stdout,*) nsym

   write(stdout, '(5X, "nbndmin ", i4)'), nbndmin
   write(stdout, '(5X, "Ef ", f12.7, " bandwidth ", f12.7)'), ef*rytoev, bandwidth*rytoev
   write(stdout, '(5X, "N(0) ", f12.7)'),  N0/rytoev
   write(stdout, '(5X, "debye temp Ry", f12.7)'), debye_e
   write(stdout, '(5X, "\mu", f12.7)'), mu
   write(stdout, '(5X, "\mu^{*} ", f12.7)'), mu/(1+mu*log((ef)/debye_e))
   write(stdout, '(5X, "\mu^{*}", f12.7)' ) mu/(1+mu*log((bandwidth)/debye_e))

   write(stdout, * ) 
   write(stdout, * ) munnp(:,:)
   write(stdout, * ) 
   write(stdout, * ) dosnnp(:)
   write(stdout, * ) 

   write(stdout, '(5X, 6f12.7)' ) munnp(1:2,1:2)

   write(stdout, '(5X, 6f12.7)' ) dosnnp(1:2)

  if (ionode) CALL davcio (muk, lrcoulmat, iuncoulmat, 1, 1)
  write (stdout, '(5X, "nbndmin ", i4)'), nbndmin

END SUBROUTINE coulmatsym

SUBROUTINE load_coul(vc, iq)
  USE kinds,           ONLY : DP
  USE klist,           ONLY : xk
  USE control_coulmat, ONLY : ngcoul, do_lind
  USE cell_base,            ONLY : tpiba2, omega, at, alat
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat, iuncoul, lrcoul
  USE gvect,                ONLY : g
  USE mp_world,         ONLY : world_comm, mpime

IMPLICIT NONE

  REAL(DP)     :: lind_eps
  INTEGER      :: iq, ig, igp
  COMPLEX(DP)  :: vc(ngcoul,ngcoul)
  LOGICAL      :: limq
  REAL(DP)     :: xq(3)
  REAL(DP)     :: qg2

  xq(:) = xk(:,iq)
  vc    = (0.0d0,0.0d0)
  if(.not.do_lind) then 
  !   print*, "reading xq:  ", xq, iq
     CALL davcio (vc, lrcoul, iuncoul, iq, -1)
     do ig = 1, ngcoul
        vc(ig,ig) = vc(ig,ig) + 1.0d0
     enddo 
  endif

  do ig = 1, ngcoul
     qg2 = (g(1,ig) + xq(1))**2 + (g(2,ig) + xq(2))**2 + (g(3,ig)+xq(3))**2
     limq = (qg2.lt.eps8) 
     if (limq) cycle
     if (.not.do_lind) then
        do igp = 1, ngcoul
           vc(ig,igp) = vc(ig,igp)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
        enddo
     else
           vc(ig,ig) = (1.0d0/lind_eps(qg2))*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
     endif
  enddo
END SUBROUTINE load_coul

REAL(DP) FUNCTION lind_eps(qg2)
  USE kinds,       ONLY : DP
  USE dielectric,  ONLY : qtf, kf
  IMPLICIT NONE
  REAL(DP) :: qg2, x
  x = qg2/(2.0d0*kf)
  lind_eps = 1.0d0 + (qtf)**2.0d0/(2.0d0*qg2**2)*(1.0d0+(1.0d0/(2.0d0*x))*(1.0d0-x**2)*log(abs((1.0d0+x)/(1.0d0-x))))
  RETURN
END FUNCTION

! Copyright (C) 2004-2009 Andrea Benassi and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!------------------------------
SUBROUTINE coulmatsym()
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE cell_base,            ONLY : tpiba2, omega, at, alat
  USE ener,                 ONLY : ef
  USE io_global,            ONLY : ionode, stdout
  USE start_k,              only : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE klist,                ONLY : nks, nkstot, ngauss, degauss, xk, wk
  USE wavefunctions_module, ONLY : evc 
  USE wvfct,                ONLY : nbnd, et, npw, igk, npwx,  g2kin, ecutwfc
  USE lsda_mod,   ONLY : nspin
  USE ktetra,     ONLY : ntetra, tetra, ltetra
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE start_k,              only : nks_start, xk_start, wk_start, &
                                   nk1, nk2, nk3, k1, k2, k3
  USE io_files,             ONLY : nwordwfc, iunwfc, diropn
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE gvect,         ONLY : ngm, g, nl
  USE gvecs,                ONLY : nls, nlsm
  USE control_coulmat,      ONlY : degaussfs, nbndmin, debye_e, do_lind, ngcoul

IMPLICIT NONE

  COMPLEX(DP) :: vcnknpkp
  REAL(DP)    :: enk, enpkp
  REAL(DP)    :: norm, nqs
  REAL(DP)    :: En, DOSofE(2), N0
  REAL(DP)    :: degaussw0, w0g1, w0g2
  REAL(DP)    :: qg2, xq(3), xkp(3)
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  COMPLEX(DP) :: psink(dffts%nnr,nbnd), psinpkp(dffts%nnr,nbnd), psi_temp(dffts%nnr), fnknpkp(dffts%nnr)
  INTEGER     :: ik, ikp, ibnd, jbnd, ig, igp, isym
  INTEGER     :: iq
  INTEGER     :: nkp, nkp_abs, ipool
!IO
  LOGICAL     :: exst
!to put in coul struct:
  REAL(DP)    :: mu, mustar
! SYMMMETRY
  LOGICAL     :: found_q, inv_q, minus_q
  INTEGER     :: iqrec, irotmq
  REAL (DP)   :: gi(3,48), gimq(3)
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
  COMPLEX(DP), ALLOCATABLE :: vc(:,:)
  INTEGER :: ir

  ALLOCATE ( gmapsym  (ngm, nrot)   )
  ALLOCATE ( eigv     (ngm, nrot)   )
  CALL gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
  ALLOCATE(vc(ngcoul,ngcoul))

 !First calculate dosef: 
 !dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
 !N(Ef) in the equation for lambda is the DOS per spin
 !dosef = dosef / two
  WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
       &        "ngauss,degauss=",i4,f12.6/)') ngauss, degauss
  !CALCULATE Dos at Fermi Level gaussian:
   CALL dos_g(et,nspin,nbnd, nks,wk,degauss,ngauss, ef, DOSofE)
  !tetra=.true.
  !use tetrahedron method:
  !CALL dos_t(et, nspin, nbnd, nks, ntetra, tetra, ef, DOSofE)
  !want density of states per spin.
  N0 = DOSofE(1)/2.d0
  if (degaussfs.eq.0.0d0) then
     degaussw0 = degauss
  else
     degaussw0 = degaussfs
  endif
  En      = ef
! nqs     = dble(nk1*nk2*nk3)
  nqs     = nks
  mu = 0.0d0
!Loop over IBZ on kpoints:
  print*, xk(:,1:nks)
  print*, wk(1:nks)

  do iq = 1, nks
     do ik = 1, nks
!    if (mod(ik,2).eq.0) print*, "xk:  ", ik
!Load coulomb interaction at xq
        vc = (0.0d0, 0.0d0)
        call load_coul(vc, iq)
!do S = nsymops:
!need to modify this:
!        CALL ktokpmq (xk(1,ik), xq(1), -1, ipool, nkp, nkp_abs)
!On a gen'l grid ktokpmq should also return the symmetry operation which rotates
!xkp to the brillouin zone. then we keep track of that symmetry else and use it
!to invert one of kpoints in the IBZ. Keeping with the idea that we only need
!IBZk and IBZq.
        xq(:) = xk(:,iq)
        xkp(:) = xk(:,ik)-xq(:)
!       print*, 'xkp'
!       print*, xkp
!Find symmop what gives us k'
        inv_q=.false.
        call find_qG_ibz(xkp, s, iqrec, isym, found_q, inv_q)
        print *,'iqrec', iqrec, isym
!       print*, xk(:, iqrec)
        ikp = iqrec
        psink(:,:)   = (0.0d0,0.0d0)
        psinpkp(:,:) = (0.0d0,0.0d0)
        CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
        psink(nls(igk(1:npw)), :) = evc(1:npw,:)
        CALL gk_sort (xkp(1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, 2*nwordwfc, iunwfc, ikp, -1 )
!       psinpkp(nls(igk(1:npw)),:) = evc(1:npw,:)
!SYMM
        psinpkp(nls(gmapsym(igk(1:npw), isym)),:) = evc(1:npw,:)

        do ibnd = nbndmin, nbnd
           CALL invfft ('Wave', psink(:,ibnd), dffts)
           do jbnd = nbndmin, nbnd
!calculate f_{nk,npkp}(\G)
              psi_temp = (0.0d0, 0.0d0)
              psi_temp = psinpkp(:,jbnd)
              CALL invfft ('Wave', psi_temp(:), dffts)
              fnknpkp = (0.0d0,0.0d0)
              do ir = 1, dffts%nnr  
                 fnknpkp(ir) = fnknpkp(ir) +  conjg(psi_temp(ir))*psink(ir,ibnd)
              enddo
              CALL fwfft ('Wave', fnknpkp(:), dffts)
!Eigenvalues
              enk = (et(ibnd, ik) - ef)
              enpkp = (et(jbnd, ikp) - ef)
              !print*, enk, enpkp
              w0g1 = w0gauss ( enk / degaussw0, 0) / degaussw0
              w0g2 = w0gauss ( enpkp / degaussw0, 0) / degaussw0
              !print*, w0g1, w0g2
!Again we want per spin.
              vcnknpkp = 0.0d0
              if(.not.do_lind) then
                do ig = 1, ngcoul 
                   do igp = 1, ngcoul
                      vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,igp)*fnknpkp(nls(ig))
                   enddo
                enddo
              else
                do ig = 1, ngcoul 
                   vcnknpkp = vcnknpkp + conjg(fnknpkp(nls(ig)))*vc(ig,ig)*fnknpkp(nls(ig))
                enddo
              endif
!mu = mu + (1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)
!wk weight includes 2*spin index... so for Coulomb we kill that...
!and to get per spin we need to kill it.
              mu = mu +(1.0d0/N0)*vcnknpkp*w0g1*w0g2*(wk(ik)/2.0)*(wk(iq)/2.0)
          enddo!jbnd
       enddo!ibnd
     enddo!iq
   enddo!ik
!CALL mp_sum(mu, inter_pool_comm)
!Factors not included when we calculate V^{c}_{nkn'k'}.
!mu = mu/omega/(2*nk1*nk2*nk3)
  mu = mu/omega
  print*, nk1,nk2,nk3
  print*,omega
  print*, "Ef", ef*rytoev, "N(0)", N0/rytoev
  print*, "debye temp Ry", debye_e
  print*, "\mu", mu
  print*, "\mu^{*}", mu/(1+mu*log((ef)/debye_e))
END SUBROUTINE coulmatsym

SUBROUTINE load_coul(vc, iq)
  USE kinds,           ONLY : DP
  USE klist,           ONLY : xk
  USE control_coulmat, ONLY : ngcoul, do_lind
  USE cell_base,            ONLY : tpiba2, omega, at, alat
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat, iuncoul, lrcoul
  USE gvect,                ONLY : g
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
     print*, "reading xq:  ", xq, iq
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
  x = qg2/(2*kf)
  lind_eps = 1.0d0 + (qtf)**2.0/(2.0*qg2**2)*(1.0d0+(1.0d0/(2.0*x))*(1-x**2)*log(abs((1+x)/(1-x))))
  RETURN
END FUNCTION

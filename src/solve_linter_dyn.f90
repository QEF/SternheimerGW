  !
  !-----------------------------------------------------------------------
  subroutine solve_linter_dyn ( dvbare, dvscf, xxq, et, vr, w, maxscf, maxbcgsolve, alpha_mix, tprec, convt)
  !-----------------------------------------------------------------------
  !
  ! this works for one perturbation at a time
  !
  ! Fully converge dpsi for each one of the first "maxscf" dvscf iterations,    
  ! then perform only "maxbcgsolve" steps for dpsi at each dvscf iteration.
  ! This should be similar to the "direct minimization method" in Paratec
  ! setting maxscf = 10**10 we have the conventional minimization
  ! (full relaxation at every scf step)
  !
  !-----------------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
#ifdef __PARA
  use para, only : mypool
#endif
  implicit none
  !
  integer :: maxscf, maxbcgsolve
  logical :: tprec
  ! switches on/off preconditioning
  !
  complex(kind=DP) :: dvbare(nr)
  ! the perturbation in real space
  complex(kind=DP) :: dvscf(nr), vr(nr)
  real(DP) :: et(nbnd_occ, nks), w, alpha_mix
  !
  complex(kind=DP) :: vr_dyn (nr)
  ! local potential plus the dynamical part w + i * eta
  !
  integer :: ik, ikk, ikq, iter, ibnd, jbnd, ios, ig, ir
  complex(kind=DP) :: evc (ngm, nbnd_occ), evq (ngm, nbnd_occ)
  real(kind=DP) :: g2kin(ngm), dr2, wgt, qg2, xxq(3), x
  !
  complex(kind=DP) :: dpsi(ngm,nbnd_occ), dvpsi(ngm,nbnd_occ), &
                      dvpsi0(ngm,nbnd_occ), dvscfin(nr), hpsi(ngm), &
                      hpsi2(ngm,nbnd_occ), ps2(nbnd_occ,nbnd_occ), &
                      drhoscf(nr), dvscfout(nr)
  complex(kind=DP) :: dpsip(ngm,nbnd_occ), dpsim(ngm,nbnd_occ)
  complex(kind=DP) :: aux(nr), aux1(nr), aux2(nr)
  complex(kind=DP) :: ps(nbnd_occ), auxg(ngm)
  complex(kind=DP) :: ZDOTC, cw
  real(DP) :: eprec(nbnd_occ), h_diag(ngm, nbnd_occ), anorm, meandvb
  logical :: conv_root, convt
  integer :: lter
  external ch_psi_all_eta

  !
  !  loop over the iterations
  !
  iter = 0
  convt = .false.
  do while (iter.lt.nmax_iter .and. .not.convt)
     !
     iter = iter + 1
     lter = 0
     drhoscf = czero
     !
     do ik = 1, nksq
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! reads unperturbed wavefuctions u(k) and u(k+q)
        !
        read ( iunwfc, rec = ikk, iostat = ios) evc
        read ( iunwfc, rec = ikq, iostat = ios) evq
        !
        ! compute the kinetic energy for k+q
        !
        do ig = 1, ngm
          g2kin(ig) = ( (xk(1,ikq) + g(1,ig))**2.d0 + &
                        (xk(2,ikq) + g(2,ig))**2.d0 + &
                        (xk(3,ikq) + g(3,ig))**2.d0 ) * tpiba2
        enddo
        !
        if (iter.eq.1) then
          !
          !  dpsi and dvscfin are set to zero
          !
          dpsip = czero
          dpsim = czero
          dvscfin = czero
          !
          ! dvbare*psi is calculated for this k point and all bands...
          ! (we compute the product in real space)
          !
          dvpsi = czero
          !
          do ibnd = 1, nbnd_occ
            call dv_psi ( dvbare, evc (1:ngm, ibnd), dvpsi (1:ngm,ibnd) )
          enddo
          !
          ! writes dvpsi for this k-point on iunit iubar
          !
          write ( iubar, rec = ik, iostat = ios) dvpsi
          !
        else
          !
          ! read  dvbare*psi for this k-point on iunit iubar
          !
          read ( iubar, rec = ik, iostat = ios) dvpsi
          !
          ! dvpsi =  dvbare*psi + dvscfin*psi 
          !
          do ibnd = 1, nbnd_occ
            call dv_psi ( dvscfin, evc (1:ngm, ibnd), dvpsi (1:ngm,ibnd) )
          enddo
          !
          !  read dpsi for this k-point from iudwf
          !  this saves us about 10% of the total CG iterations on 1 scf cycle
          !  (although we cannot use it if we want to go for multishift solvers)
          !
          read ( iudwfp, rec = ik, iostat = ios) dpsip
          read ( iudwfm, rec = ik, iostat = ios) dpsim
          !
          ! initial guess set to zero - this seems to work better
          ! than the reading from previous scf cycle
          !
   !      dpsip = czero
   !      dpsim = czero
          !
        endif
        !
        !  - ( 1 - P_occ^{k+q} ) * dvpsi
        !
        call emptyproj ( evq, dvpsi)       
        !
        ! iterative solution of the linear system 
        ! (H-et)*dpsi=   - ( 1 - P_occ^{k+q} ) * (dvbare+dvscf)*psi
        !              [                  dvpsi                     ]
        !
        if (tprec) then
          !
          ! I use the preconditioner of Teter, Payne, Allan [PRB 40, 12255 (1988)]
          !
          do ibnd = 1, nbnd_occ
             do ig = 1, ngm
                 auxg (ig) = g2kin (ig) * evq (ig, ibnd)
             enddo
             ! the misterious factor 1.35 seems to be a way to match 
             ! approxiately the TPA preconditioner and the simple one 1/x. FG
      !      eprec (ibnd) = 1.35d0 * ZDOTC (ngm, evq (1, ibnd), 1, auxg, 1)
             eprec (ibnd) = ZDOTC (ngm, evq (1, ibnd), 1, auxg, 1)
          enddo
          do ibnd = 1, nbnd_occ
             do ig = 1, ngm
                x = g2kin(ig)/eprec (ibnd)
                h_diag(ig,ibnd) = (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0) &
                                / (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0+16.d0*x**4.d0)
                !
                ! original preconditioning used in PH 
                !
                ! h_diag(ig,ibnd) = 1.d0/(max(1.d0,x)-et(ibnd,ikk)+w) ! preconditioner as in VdW code
                ! 
             enddo
          enddo
          !
        else
          ! NO preconditioning
          !
          h_diag = 1.d0
          !
        endif
        !
        ! --- now obtain dpsi from cBiCG
        !
        ! [ The following works for any complex frequency: in order to obtain
        !   points on the real axis, just set cw = dcmplx (w, 0.d0) 
        !
        cw = dcmplx (0.d0, w)
     !  cw = dcmplx ( w, 0.00734981) ! 0.1 eV imaginary component
        !
        ! complex frequency +w
        !
        call bcgsolve_all   (ch_psi_all_eta, et(:,ikk), dvpsi, dpsip, h_diag, &
             ngm, ngm, tr_cgsolve, ik, lter, conv_root, anorm, nbnd_occ, &
             g2kin, vr, evq,  cw )
        !
        ! complex frequency -w
        !
        call bcgsolve_all   (ch_psi_all_eta, et(:,ikk), dvpsi, dpsim, h_diag, &
             ngm, ngm, tr_cgsolve, ik, lter, conv_root, anorm, nbnd_occ, &
             g2kin, vr, evq, -cw )
        !
!       if (.not.conv_root) &
!          write( 6, '(4x,"ik",i4," linter: one or more roots not converged ",e10.3)') &
!          ik , anorm
!       write(6,'("cgsolve_all:",2x,3i5,3x,e9.3)') iter, ik, lter, anorm
!       !
!       !
!       ! DEBUG: calculate (H-et+alpha*Pv)*dpsi-dvpsi, this should be zero 
!       ! if dpsi is the correct solution - O.K. this is checked
!       !
!       call ZGEMM ('C', 'N', nbnd_occ , nbnd_occ, ngm, (1.d0, 0.d0) , evq, &
!         ngm, dpsi, ngm, (0.d0, 0.d0) , ps2, nbnd_occ)
!       call ZGEMM ('N', 'N', ngm, nbnd_occ, nbnd_occ, dcmplx(alpha_pv,0.d0), evq, &
!         ngm, ps2, nbnd_occ, czero, hpsi2, ngm)
!       do ibnd = 1, nbnd_occ
!         call h_psi_c ( dpsi(:,ibnd), hpsi, g2kin, vr_dyn)
!         hpsi = hpsi - et(ibnd,ikk) * dpsi(:,ibnd) - dvpsi(:,ibnd)
!         hpsi = hpsi +  hpsi2 (:, ibnd)
!         write(6,*) '--------------------'
!         do ig = 1, ngm
!           write(6,'(3i5,3(2x,2f15.5))') &
!             ik, ibnd, ig, 100*dpsi(ig,ibnd), 100*dvpsi(ig,ibnd), 100*hpsi(ig)
!         enddo
!       enddo
        !
        ! writes dpsi for this k point on iunit iudwf
        !
        write ( iudwfp, rec = ik, iostat = ios) dpsip
        write ( iudwfm, rec = ik, iostat = ios) dpsim
        !
        ! contribution to drhoscf from this kpoint
        !
        wgt = wk(ikk) / omega
!
! My interpretation of the factor 1/omega (which in principle
! does not exist in my notes):
! The fft is performed in fftw.f90 as
! G->r: sum_G c(G) exp(iGr)
! r->G: sum_r psi(r) exp(-iGr) / Nr 
! while what we really need is the following:
! psi(r) = 1/sqrt(Omega) sum_G c(G) exp(iGr)
! c(G) = 1/sqrt(Omega) sum_r psi(r) exp(-iGr) / Nr
! Therefore in order to calculate the induced charge density
! in real space we need to multiply by 1/sqrt(Omega) for evc.
! At the same time, the linear system which is solved has
! -(1-P_occ)dV psi on the rhs, and psi is the "wrong" one,
! in the sense that it misses the 1/sqrt(Omega) factor.
! As a consequence, the dpsi misses the 1/sqrt(Omega) factor
! and we need to bring it back. That is why we multiply drho by 1/Omega
! NOTE: this will apply to the transformation of W and G as well.
!
        !
        do ibnd = 1, nbnd_occ
          !
          ! psi^\star * ( dpsip + dpsim )
          !
          aux1 = czero
          do ig = 1, ngm
            aux1 ( nl ( ig ) ) = evc (ig, ibnd)
          enddo
          call cfft3 ( aux1, nr1, nr2, nr3,  1)
          !
          aux2 = czero
          do ig = 1, ngm
            aux2 ( nl ( ig ) ) = dpsip (ig, ibnd) + dpsim (ig, ibnd)
          enddo
          call cfft3 ( aux2, nr1, nr2, nr3,  1)
          !
          do ir = 1, nr
            drhoscf (ir) = drhoscf (ir) + wgt * conjg (aux1 (ir) ) * aux2 (ir)
          enddo
          !
        enddo
        !
     enddo 
     !
     ! here we have drhoscf for this iteration
     ! compute the corresponding Hartree potential (RPA)
     !
     dvscfout = czero
     call cfft3 ( drhoscf, nr1, nr2, nr3, -1)
     !
     ! here we enforce zero average variation of the charge density 
     ! if the bare perturbation does not have a constant term
     ! (otherwise the numerical error, coupled with a small denominator
     ! in the coulomb term, gives rise to a spurious dvscf response)
     !
     meandvb = sqrt ( (sum(dreal(dvbare)))**2.d0 + (sum(aimag(dvbare)))**2.d0 ) / float(nr)
     if (meandvb.lt.1.d-8) drhoscf ( nl(1) ) = 0.d0
     !
     do ig = 1, ngm
       qg2 = (g(1,ig)+xxq(1))**2 + (g(2,ig)+xxq(2))**2 + (g(3,ig)+xxq(3))**2
       if (qg2 > 1.d-8) &
         dvscfout ( nl(ig) ) =  e2 * fpi * drhoscf ( nl(ig) ) / (tpiba2 * qg2)
     enddo
     !
     call cfft3 ( dvscfout, nr1, nr2, nr3,  1)
     !
     ! we mix with the old potential
     !
     ! modif broyden mixing, complex potential
     !
     call mix_potential_c ( nr, dvscfout, dvscfin, alpha_mix, dr2, tr2_ph, iter, nmix_ph, convt)
     !
     ! convergence criterion
     !
     convt = dr2.lt.tr2_ph
     !
!    write(6,'(4x, "scf iteration ",i3,": dr2 = ",e8.2)') iter, dr2
     !
     ! the factor 2 below is because we actually run +iw and -iw
     ! so the average number for each one of +iw and -iw is lter/nksq/2 
     !
     write(6,'(4x, "scf iteration ",i3,": dr2 = ",e8.2,3x,"average CG iter ",f5.1)') &
       iter, dr2, float(lter)/float(nksq)/float(2)
#ifdef __PARA
!    write(1000+mypool,'(4x, "scf iteration ",i3,": dr2 = ",e8.2,3x,"average CG iter ",f5.1)') &
!      iter, dr2, float(lter)/float(nksq)/float(2)
#endif
     !
  enddo
  !
  ! at this point dvscfin is the converged Hartree screening.
  ! the screened coulomb interaction corresponds to dv_bare + dv_hartree (RPA)
  !
  dvscf = dvscfin + dvbare
  !
  return
  end subroutine solve_linter_dyn
  !
  !-----------------------------------------------------------------------
  subroutine dv_psi (dv, psi, dvpsi)
  !-----------------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
  implicit none
  !
  complex(kind=DP) :: dv(nr), dvpsi(ngm)
  complex(kind=DP) :: dvscf(nr)
  integer :: ibnd, ig, ir
  complex(kind=DP) :: psi (ngm)
  complex(kind=DP) :: aux(nr)
  !
  aux = czero
  do ig = 1, ngm
    aux ( nl ( ig ) ) = psi (ig)
  enddo
  call cfft3 ( aux, nr1, nr2, nr3,  1)
  do ir = 1, nr
    aux (ir) = aux(ir) * dv (ir)
  enddo
  ! back to G-space (fft order of G-vectors)
  call cfft3 ( aux, nr1, nr2, nr3, -1)
  ! switch to magnitude-order of G-vectors
  do ig = 1, ngm
    dvpsi(ig) = dvpsi(ig) + aux( nl(ig) )
  enddo
  !
  end subroutine dv_psi
  !-----------------------------------------------------------------------
  subroutine emptyproj ( evq, dvpsi )
  !-----------------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
  implicit none
  !
  complex(kind=DP) :: evq (ngm, nbnd_occ), dvpsi(ngm,nbnd_occ) 
  complex(kind=DP) :: ps(nbnd_occ), auxg(ngm)
  integer :: ibnd, jbnd
  complex(kind=DP) :: ZDOTC
  !
  do ibnd = 1, nbnd_occ
    auxg = czero
    do jbnd = 1, nbnd_occ
      ps(jbnd) = - ZDOTC(ngm, evq(:,jbnd), 1, dvpsi(:,ibnd), 1)
      call ZAXPY (ngm, ps (jbnd), evq (:, jbnd), 1, auxg, 1)
    enddo
    call DAXPY (2 * ngm, one, auxg, 1, dvpsi (:, ibnd), 1)
  enddo
  call DSCAL (2 * ngm * nbnd_occ, - 1.d0, dvpsi, 1)
  !
  end subroutine emptyproj

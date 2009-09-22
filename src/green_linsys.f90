  !
  !----------------------------------------------------------------
  subroutine green_linsys ( vr, g2kin, xxk, nw, w, green, igstart, igstop)
  !----------------------------------------------------------------
  ! 
  use parameters
  use constants
  use gspace
  use kspace
!#ifdef __PARA
!  USE para
!  USE mp_global,  ONLY : nproc, mpime, nproc_pool, my_pool_id, me_pool
!  USE mp, ONLY:  mp_barrier
!#endif
  implicit none
  !
  real(dbl) :: xxk(3)
  integer :: ig, igp, nw
  integer :: igstart, igstop
  integer :: lter, n1
  real(DP) :: h_diag(ngm), anorm, rdummy
  logical :: conv_root, tprec
  integer :: ig1, ig2, idummy, iw, ierr, ibnd, ios, recl, unf_recl, ikf, fold(nq)
  real(dbl) :: x, delta, kplusg(3), g2kin(ngm), eval(nbnd_occ)
  real(dbl) :: eval_all(nbnd), w(nw), w_ryd(nw)
  complex(dbl) :: vr(nr), psi(ngm, nbnd_occ), rhs(ngm), gr(ngm), cw, gr_A(ngm), gr_N(ngm)
  complex(dbl) :: psi_all(ngm,nbnd), gr_exp(nw)
  complex(dbl) :: green(nrs,nrs,nw), cdummy(ngm, nbnd_occ)
  logical :: convt, found
  complex(kind=DP) ::  auxg(ngm)
  complex(kind=DP) :: ZDOTC
  real(DP) :: eprec(nbnd_occ), dirac
  external ch_psi_all_eta

  call start_clock('green_linsys')
 
!  !
!  do iw=1,nw
!    w(iw) = float(iw-1)/float(nw-1)*50-25
!  enddo

  w_ryd = w/ryd2ev

  lter = 0
  n1 = 1
  rdummy = 0.d0
  cdummy = (0.d0,0.d0)
  idummy = 0

  !
  !  generate occupied wavefunctions for the k0-q point
  !
! write(6,'(/4x,a)') repeat('-',67)
! write(6,'(4x,"k0-q point: occupied eigenvalues (eV)")')
! write(6,'(4x,a/)') repeat('-',67)
  !
  ! the k-dependent kinetic energy in Ry
  !
  do ig = 1, ngm
    kplusg = xxk + g(:,ig)
    g2kin ( ig ) = tpiba2 * dot_product ( kplusg, kplusg )
  enddo
  !
  call eigenstates2 ( xxk, vr, g2kin, psi, eval )
  !
  write ( stdout, '(4x,"k0-q = (",3f12.7," )",10(3x,f7.3))') xxk, eval*ryd2ev
  !
  tprec = .true.
  if (tprec) then
    !
    ! I use the preconditioner of Teter, Payne, Allan [PRB 40, 12255 (1988)]
    ! I set the kinetic energy reference to the highest occupied state
    !
    ibnd = nbnd_occ
    do ig = 1, ngm
        auxg (ig) = g2kin (ig) * psi (ig, ibnd)
    enddo
    eprec (ibnd) = ZDOTC (ngm, psi (1, ibnd), 1, auxg, 1)
    do ig = 1, ngm
       x = g2kin(ig)/eprec (ibnd)
       h_diag(ig) = (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0) &
                   / (27.d0+18.d0*x+12.d0*x*x+8.d0*x**3.d0+16.d0*x**4.d0)
    enddo
    !
  else
    ! NO preconditioning
    !
    h_diag = 1.d0
    !
  endif
  !
  !  loop on bare perturbations ig and fixed k point
  !
  !  [note: no barriers inside this loop as the number of perturbation
  !  is different across processors]
  !
  do iw = 1, nw
    !
!   write(6,'(4x,3x,"iw = ",i5," of ",i5)') iw,nw
    !
    do ig = igstart, igstop 
      !
!     write(6,'(4x,"ig = ",i5)') ig
      !
!     write(6,'(4x,"Green linsys: k =",3f7.3,"  G =",3f7.3,"  w(eV) =",3f7.3)') &
!        xxk,g(:,ig), w(iw)
      !
      cw = dcmplx(w_ryd(iw),eta)
      !
      rhs = czero
      !
      ! rhs must be -1: the spin factor does not enter since in <nk|Sigma|nk>
      ! only states of the same spin couple [cf HL86]
      !
!     rhs ( ig ) = (-2.d0,0.d0)
      rhs ( ig ) = -cone
      !
      lter = 0
      gr_A = czero
      call bcgsolve_all (ch_psi_all_eta, rdummy, rhs, gr_A, h_diag, &
         ngm, ngm, tr_cgsolve, idummy, lter, conv_root, anorm, n1, g2kin, vr, cdummy, cw)
      !
!     write(6,'(4x,"bcgsolve_all iterations:",i5)') lter
      !
      gr_N = czero
      do ibnd = 1, nbnd_occ
         x = w_ryd(iw) - eval(ibnd)
         dirac = eta / pi / (x**2.d0+eta**2.d0)
         !
         ! no spin factor (see above)
!        gr_N = gr_N + 2.d0 * twopi * ci * conjg ( psi(ig,ibnd) ) * psi(1:ngm,ibnd) * dirac
         gr_N = gr_N + twopi * ci * conjg ( psi(ig,ibnd) ) * psi(1:ngm,ibnd) * dirac
      enddo
      !
      gr = gr_A + gr_N
      !
!     write(400+mypool,*) w(iw),real(gr(4)),aimag(gr(4))
      !
      ! keep only the G-vectors 1:ngms for the Green's function
      !
      do igp = 1, ngms
        green (ig,igp,iw) = gr (igp) 
      enddo
      !
      ! green (-,igp,-) is in SIZE order of G-vectors
      !
    enddo
    !
  enddo
  !
!  ! direct calculation - debug only
!  ! The Sternheimer method and the direct calculation match perfectly
!  !
!  ig1 = 2
!  ig2 = 4
!  call  eigenstates_all ( vr, g2kin, psi_all, eval_all )
!  !
!  gr_exp = czero
!  do ibnd = 1, nbnd
!    if (ibnd.le.nbnd_occ) then
!      delta =  1.d0
!    else
!      delta = -1.d0
!    endif
!    do iw = 1, nw
!     gr_exp(iw) = gr_exp(iw) + psi_all(ig1,ibnd)*conjg(psi_all(ig2,ibnd)) &
!                / ( w_ryd(iw) - eval_all(ibnd) - ci * delta * eta)
!    enddo
!  enddo
!  gr_exp = 2.d0 * gr_exp ! spin
!  do iw = 1, nw
!    write(300+mypool,'(3f15.10)') w(iw), gr_exp(iw)
!  enddo
  !
  call stop_clock('green_linsys')
  ! 
  return
  end subroutine green_linsys
  !----------------------------------------------------------------
  !

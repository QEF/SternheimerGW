  !
  !----------------------------------------------------------------
  subroutine sigma_matel (  ik0, vr, xk0, nw, w )
  !----------------------------------------------------------------
  ! 
  use parameters
  use constants
  use gspace
  use kspace
#ifdef __PARA
  USE para
  USE mp_global,  ONLY : nproc, mpime, nproc_pool, my_pool_id, me_pool
  USE mp, ONLY:  mp_barrier
#endif
  implicit none
  !
  integer :: ig, igp, nw, iw, ibnd, jbnd, ios, ipol, ik0, ir
  real(dbl) :: xk0(3), kplusg(3), g2kin(ngm), et(nbnd), w(nw), w_ryd(nw), v_xc(nr)
  complex(dbl) :: vr(nr), evc(ngm,nbnd), sigma(ngms,ngms,nw), aux(ngms), vpsi(ngm), psic(nr)
  complex(kind=DP) :: ZDOTC, sigma_band(nbnd_sig,nbnd_sig,nw), vxc(nbnd_sig,nbnd_sig)
  real(dbl) :: resig_diag(nw,nbnd_sig), imsig_diag(nw,nbnd_sig), et_qp(nbnd_sig), a_diag(nw,nbnd_sig)
  real(dbl) :: resig_diag_tr(nw), imsig_diag_tr(nw), a_diag_tr(nw), et_qp_tr, z_tr, z(nbnd_sig)
  integer :: iman, nman, ndeg(nbnd_sig), ideg
  !
  call start_clock ('sigma_matel')
  !
  w_ryd = w/ryd2ev
  !
  write(stdout,'(/4x,"k0(",i3," ) = (",3f7.3," )")') ik0, (xk0 (ipol) , ipol = 1, 3)
  !
  ! NOTE - I calculate the eigenstates of -xk0
  ! in order to have c_k(-G) = [c_-k(G)]*
  ! Because of my convention on the FFTs on G,G' in the paper,
  ! below we mix G and -G in the sandwitches. The easiest way
  ! to perform the calculation is to use teh eigenvectors for -xk0  
  ! and take their cc to obtain c(-G) for xk0
  !
  do ig = 1, ngm
    ! note the -xk0 for the reason above!
    kplusg = -xk0 + g(:,ig)
    g2kin ( ig ) = tpiba2 * dot_product ( kplusg, kplusg )
  enddo
  !
  call  eigenstates_all ( vr, g2kin, evc, et )
  !
  write(stdout,'(4x,"LDA eigenval (eV)",8(1x,f7.3))') et(1:nbnd_sig)*ryd2ev
  !
#ifdef __PARA
  ! only proc 0 reads from file and does the product
  ! (need some sort of parallelization here)
  if (me.eq.1.and.mypool.eq.1) then
#endif
  !
  ! MATRIX ELEMENTS OF THE XC POTENTIAL
  !
  open(unit=110,file='vxc.dat')
  rewind(110)
  do ir = 1, nr
    read(110,*) v_xc(ir)
  enddo 
  close(110)
  !
  do jbnd = 1, nbnd_sig
    !
    psic = czero
    do ig = 1, ngm
      psic ( nl ( ig ) ) = evc(ig, jbnd)
    enddo
    call cfft3 ( psic, nr1, nr2, nr3,  1)
    do ir = 1, nr
      psic (ir) = psic(ir) * v_xc (ir)
    enddo
    call cfft3 ( psic, nr1, nr2, nr3, -1)
    do ig = 1, ngm
      vpsi(ig) = psic( nl(ig) )
    enddo
    !
    do ibnd = 1, nbnd_sig
       vxc(ibnd,jbnd) = ZDOTC (ngm, evc (1, ibnd), 1, vpsi, 1)
    enddo
    !
  enddo
  !
  ! MATRIX ELEMENTS OF THE SELF-ENERGY
  !
  read ( iunsigma, rec = ik0, iostat = ios) sigma
!@
!sigma = sigma / 2.d0 ! HAL silicon.sigma had a wrong spin factor 2
!@
  !
  ! following the convention in the paper, thsi should be
  ! <i|Sigma|j> = sum_G,G' u_ik^*(-G) <G|Sigma|G'> u_jk(-G')
  !             = sum_G,G' u_i,-k(G) <G|Sigma|G'> [u_j,-k(G')]*
  !
  do ibnd = 1, nbnd_sig 
   do jbnd = 1, nbnd_sig
    do iw = 1, nw
      !
      sigma_band (ibnd, jbnd, iw) = czero
      !
      do ig = 1, ngms
        aux = sigma (ig, 1:ngms, iw) 
        sigma_band (ibnd, jbnd, iw) = sigma_band (ibnd, jbnd, iw) + &
            evc (ig, ibnd) * ZDOTC (ngms, evc (1:ngms, jbnd), 1, aux, 1)
      enddo 
      !
    enddo
   enddo
  enddo
  !
  ! Now calculate the expectation value of the self-energy
  ! using the diagonal matrix elements
  ! NOTE: we cannot calculate the QP correction until we
  ! have the expt value of the Vxc. I could get Vxc from the
  ! density and the LDA expression.
  !
  do ibnd = 1, nbnd_sig
    !
    do iw = 1, nw
      resig_diag (iw,ibnd) = real( sigma_band (ibnd, ibnd, iw) - vxc(ibnd,ibnd) )
      imsig_diag (iw,ibnd) = aimag ( sigma_band (ibnd, ibnd, iw) )
      a_diag (iw,ibnd) = one/pi * abs ( imsig_diag (iw,ibnd) ) / &
         ( abs ( w_ryd(iw) - et(ibnd) - resig_diag (iw,ibnd) )**2.d0 &
          + abs ( imsig_diag (iw,ibnd) )**2.d0 ) 
    enddo
    !
    call qp_eigval ( nw, w_ryd, resig_diag(1,ibnd), et(ibnd), et_qp (ibnd), z(ibnd) )  
    !
  enddo
  !
  ! Now take the trace (get rid of phase arbitrariness of the wfs)
  ! (alternative and more approrpiate: calculate the nondiag on the
  ! deg subspaces and diagonalize)
  !
  ! count degenerate manifolds and degeneracy...
  !
  nman = 1
  ndeg = 1
  do ibnd = 2, nbnd_sig
    if ( abs( et (ibnd) - et (ibnd-1)  ) .lt. 1.d-5 ) then
      ndeg (nman) = ndeg(nman) + 1
    else
      nman = nman + 1
    endif
  enddo
! write (stdout, *) nman, (ndeg (iman) ,iman=1,nman)
  !
  ! ...and take the trace over the manifold
  !
  ibnd = 0
  jbnd = 0
  do iman = 1, nman
    !
    resig_diag_tr = 0.d0
    imsig_diag_tr = 0.d0
    a_diag_tr = 0.d0
    et_qp_tr = 0.d0
    z_tr = 0.d0
    !
    do ideg = 1, ndeg(iman)
      ibnd = ibnd + 1
      resig_diag_tr = resig_diag_tr + resig_diag (:,ibnd)
      imsig_diag_tr = imsig_diag_tr + imsig_diag (:,ibnd)
      a_diag_tr = a_diag_tr + a_diag (:,ibnd)
      et_qp_tr = et_qp_tr + et_qp (ibnd)
      z_tr = z_tr + z (ibnd)
    enddo
    !
    do ideg = 1, ndeg(iman)
      jbnd = jbnd + 1 
      resig_diag (:,jbnd) = resig_diag_tr / float( ndeg(iman) )
      imsig_diag (:,jbnd) = imsig_diag_tr / float( ndeg(iman) )
      a_diag (:,jbnd) = a_diag_tr / float( ndeg(iman) )
      et_qp (jbnd) = et_qp_tr / float( ndeg(iman) )
      z (jbnd) = z_tr / float( ndeg(iman) )
    enddo
    !
  enddo
  !
  write(stdout,'(4x,"GW  expt val (eV)",8(1x,f7.3))') et_qp(1:nbnd_sig)*ryd2ev
  write(stdout,'(4x,"GW  renorm       ",8(1x,f7.3)/)') z(1:nbnd_sig)
  !
  do iw = 1, nw
    write(stdout,'(9f15.8)') w(iw), (ryd2ev*resig_diag (iw,ibnd), ibnd=1,nbnd_sig)
  enddo
  write(stdout,*)
  do iw = 1, nw
    write(stdout,'(9f15.8)') w(iw), (ryd2ev*imsig_diag (iw,ibnd), ibnd=1,nbnd_sig)
  enddo
  write(stdout,*)
  do iw = 1, nw
    write(stdout,'(9f15.8)') w(iw), (a_diag (iw,ibnd)/ryd2ev, ibnd=1,nbnd_sig)
  enddo
  !
#ifdef __PARA
  endif
#endif
  !
  call stop_clock ('sigma_matel')
  ! 
  return
  end subroutine sigma_matel
  !----------------------------------------------------------------
  subroutine qp_eigval ( nw, w, sig, et, et_qp, z )  
  !----------------------------------------------------------------
  !
  use parameters
  use constants
  implicit none
  integer :: nw, iw, iw1, iw2
  real(DP) :: w(nw), sig(nw), et, et_qp, dw, w1, w2, sig_et, sig1, sig2, z, sig_der
  !
  dw = w(2)-w(1)
  !
  if ((et.lt.w(1)+dw).or.(et.gt.w(nw)-dw)) &
    call errore ('qp_eigval','original eigenvalues outside the frequency range of the self-energy',1)
  iw = 1
  do while ((iw.lt.nw).and.(w(iw).lt.et))
    iw = iw + 1
    iw1 = iw-1
    iw2 = iw
  enddo
  w1 = w(iw1)
  w2 = w(iw2)
  sig1 = sig(iw1)
  sig2 = sig(iw2)
  !
  sig_et = sig1 + ( sig2 - sig1 ) * (et-w1) / (w2-w1)
  !
  sig_der = ( sig2 - sig1 ) / ( w2 - w1 )
  z = one / ( one - sig_der)
  !
  ! temporary - until I do not have Vxc
  !
  et_qp = et + z * sig_et
  !
  end subroutine qp_eigval
  !----------------------------------------------------------------
  !







